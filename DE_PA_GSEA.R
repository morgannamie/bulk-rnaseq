#!/usr/bin/env Rscript

################################################################################
#Load Libraries
################################################################################
suppressPackageStartupMessages({
  library(argparse)           
  library(DESeq2)             
  library(clusterProfiler)   
  library(org.Hs.eg.db)     
  library(ReactomePA)       
  library(enrichplot)       
  library(ggplot2)          
  library(ashr)          
})

################################################################################
#Parse Command-Line Arguments
################################################################################
#Set up argument parsing to allow flexible command-line usage of the script
parser <- ArgumentParser("Simplified DESeq2 → Pathway → GSEA pipeline")
parser$add_argument("--counts",required=TRUE, help="Counts CSV (genes × samples)")
parser$add_argument("--coldata",required=TRUE, help="Metadata CSV (rows = samples); must include 'condition'")
parser$add_argument("--control_name", required=TRUE, help="Reference condition name")
parser$add_argument("--p_adj",type="double", default=0.05, help="Adjusted p-value cutoff for significance")
parser$add_argument("--lfc", type="double", default=1.5,  help="Fold-change cutoff (linear)")
parser$add_argument("--outdir",default="bulk_rna_results", help="Output directory path")

args <- parser$parse_args()
dir.create(args$outdir, showWarnings=FALSE)  #Make sure output directory exists

################################################################################
#Read Data & Differential Expression
################################################################################
#Load input data: count matrix and sample metadata
counts  <- read.csv(args$counts, row.names=1, check.names=FALSE)
coldata <- read.csv(args$coldata, row.names=1)

#Validate input: ensure 'condition' column exists and control level is present
if (!"condition" %in% colnames(coldata)) stop("Metadata must have a 'condition' column")
if (!args$control_name %in% coldata$condition) stop("control_name not found in condition levels")

#Create DESeq2 object using ~condition design
dds <- DESeqDataSetFromMatrix(countData=counts, colData=coldata, design=~condition)

#Filter out low-count genes (counts <10 across all samples) to improve power
dds <- dds[rowSums(counts(dds)) >= 10, ]

#Run DESeq2 pipeline (normalization, dispersion estimation, testing)
dds <- DESeq(dds)

#Save the full DESeq2 object
saveRDS(dds, file.path(args$outdir, "dds_full.rds"))

#Export normalized counts for downstream visualization
norm_counts <- counts(dds, normalized=TRUE)
write.csv(as.data.frame(norm_counts), file.path(args$outdir, "normalized_counts.csv"), row.names=TRUE)

conds <- setdiff(levels(dds$condition), args$control_name)
all_res <- lapply(conds, function(lvl) {
  res <- results(dds, contrast=c("condition", lvl, args$control_name), alpha=args$p_adj)
  res <- lfcShrink(dds, contrast=c("condition", lvl, args$control_name), res=res, type="ashr")  #Shrink log2FC to reduce noise
  df  <- as.data.frame(res[order(res$padj), ])  #Sort by adjusted p-value
  write.csv(df, file.path(args$outdir, paste0("DE_results_", lvl, "_vs_", args$control_name, ".csv")), row.names=TRUE)
  df
})

#Combine results across all comparisons
res_df <- do.call(rbind, all_res)

#Filter for statistically significant and biologically relevant genes
sig <- subset(res_df, padj < args$p_adj & abs(log2FoldChange) >= log2(args$lfc))
sig_df <- sig
sig_df$ENSEMBL <- rownames(sig_df)  #Add ENSEMBL ID as a column
sig_df <- sig_df[, c("ENSEMBL", names(res_df))]  #Reorder columns

#Save significant gene list
write.csv(sig_df, file.path(args$outdir, "significant_genes.csv"), row.names=FALSE)

################################################################################
#Pathway Analysis (ORA: over-representation analysis)
################################################################################
#Convert ENSEMBL to ENTREZID for enrichment functions
sig_entrez <- bitr(unique(rownames(sig)), fromType="ENSEMBL", toType="ENTREZID", OrgDb=org.Hs.eg.db)

#Run GO, KEGG, and Reactome enrichment on significant gene set
ego <- enrichGO(gene=sig_entrez$ENTREZID, OrgDb=org.Hs.eg.db, ont="BP", pAdjustMethod="BH", pvalueCutoff=0.05, qvalueCutoff=0.1, readable=TRUE)
ekegg <- enrichKEGG(gene=sig_entrez$ENTREZID, organism="hsa", pAdjustMethod="BH", pvalueCutoff=0.05)
reactome_enrich <- enrichPathway(gene=sig_entrez$ENTREZID, pAdjustMethod="BH", pvalueCutoff=0.05, readable=TRUE)

#Export enrichment results
write.csv(as.data.frame(ego), file.path(args$outdir, "GO_BP.csv"), row.names=FALSE)
write.csv(as.data.frame(ekegg), file.path(args$outdir, "KEGG.csv"), row.names=FALSE)
write.csv(as.data.frame(reactome_enrich), file.path(args$outdir, "Reactome.csv"), row.names=FALSE)

#Generate and save dotplots of enrichment results
plot_ego <- dotplot(ego) + ggtitle("GO BP Enrichment")
plot_ekegg <- dotplot(ekegg) + ggtitle("KEGG Enrichment")
plot_reactome_enrich <- dotplot(reactome_enrich) + ggtitle("Reactome Enrichment")
ggsave(file.path(args$outdir, "go_pathway_enrichment.png"), plot=plot_ego, width=12, height=8, dpi=300)
ggsave(file.path(args$outdir, "kegg_pathway_enrichment.png"), plot=plot_ekegg, width=12, height=8, dpi=300)
ggsave(file.path(args$outdir, "reactome_pathway_enrichment.png"), plot=plot_reactome_enrich, width=12, height=8, dpi=300)

################################################################################
#GSEA Pre‑ranked with Deduplication
################################################################################
#For GSEA, rank *all* genes by log2 fold change, not just significant ones

#Convert all ENSEMBL to ENTREZID
all_eg <- bitr(rownames(res_df), fromType="ENSEMBL", toType="ENTREZID", OrgDb=org.Hs.eg.db)

#Create a ranked list of genes by log2FC (keep only one entry per ENTREZID)
df <- data.frame(ENSEMBL=rownames(res_df), log2FoldChange=res_df$log2FoldChange, stringsAsFactors=FALSE)
merged <- merge(df, all_eg, by="ENSEMBL")
merged <- merged[order(-abs(merged$log2FoldChange)), ]       #Sort by effect size magnitude
merged_unique <- merged[!duplicated(merged$ENTREZID), ]      #Keep top-ranking gene per ENTREZID

#Format named vector for GSEA input
ranks <- merged_unique$log2FoldChange
names(ranks) <- merged_unique$ENTREZID
ranks <- sort(ranks, decreasing=TRUE)

#Run GSEA for GO, KEGG, and Reactome
gsea_go <- gseGO(geneList=ranks, OrgDb=org.Hs.eg.db, ont="BP", pAdjustMethod="BH", pvalueCutoff=0.05, verbose=FALSE)
gsea_kegg <- gseKEGG(geneList=ranks, organism="hsa", pAdjustMethod="BH", pvalueCutoff=0.05, verbose=FALSE)
gsea_reactome <- gsePathway(geneList=ranks, pAdjustMethod="BH", pvalueCutoff=0.05)

#Save GSEA results
write.csv(as.data.frame(gsea_go),ffile.path(args$outdir, "GSEA_GO.csv"), row.names=FALSE)
write.csv(as.data.frame(gsea_kegg),ffile.path(args$outdir, "GSEA_KEGG.csv"), row.names=FALSE)
write.csv(as.data.frame(gsea_reactome), file.path(args$outdir, "GSEA_Reactome.csv"), row.names=FALSE)

#Plot GSEA dotplots
plot_gsea_go <- dotplot(gsea_go) + ggtitle("GSEA GO BP")
plot_gsea_kegg <- dotplot(gsea_kegg) + ggtitle("GSEA KEGG")
plot_gsea_reactome <- dotplot(gsea_reactome) + ggtitle("GSEA Reactome")
ggsave(file.path(args$outdir, "gsea_go_pathway_enrichment.png"), plot=plot_gsea_go, width=12, height=8, dpi=300)
ggsave(file.path(args$outdir, "gsea_kegg_pathway_enrichment.png"), plot=plot_gsea_kegg, width=12, height=8, dpi=300)
ggsave(file.path(args$outdir, "gsea_reactome_pathway_enrichment.png"), plot=plot_gsea_reactome, width=12, height=8, dpi=300)

message("Pipeline complete! Outputs in: ", args$outdir)
