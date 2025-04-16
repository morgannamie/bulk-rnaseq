# bulk-rnaseq

## Overview
An analysis of bulk RNAseq data utilizing a publicly available prostate cancer dataset available at https://www.ncbi.nlm.nih.gov/sra?term=SRP060235.

## From FASTQs to Counts
This pipeline uses count data generated from raw RNA-seq FASTQ files. The process includes:
- **Quality Control:** Raw reads are checked for quality using tools like FastQC
- **Trimming:** Adapters and low-quality bases are removed using tools such as Trimmomatic or Cutadapt
- **Alignment:** Cleaned reads are mapped to a reference genome using aligners like STAR or HISAT2, or pseudo-aligners like Salmon or Kallisto
- **Post-Alignment Processing:** Aligned reads are sorted, indexed, and quality-assessed
- **Quantification:** Reads overlapping gene features are counted (using tools like featureCounts or htseq-count) to create a gene-by-sample count matrix

The analysis covers:

- differential expression
- geneset enrichment
- pathway analysis

# DESeq2 → Pathway → GSEA Pipeline

This pipeline provides an end-to-end workflow for analyzing bulk RNA-seq data. It covers differential expression analysis, pathway enrichment (GO, KEGG, Reactome), and Gene Set Enrichment Analysis (GSEA) with a focus on deduplicating gene identifiers. The script is written in R and uses several essential libraries to generate statistical results and visualization plots.

## Overview

The pipeline performs the following steps:

1. **Load Required Libraries**  
   Loads libraries for argument parsing, differential expression analysis (DESeq2), enrichment analysis (clusterProfiler, ReactomePA, enrichplot), and plotting (ggplot2).

2. **Command-Line Argument Parsing**  
   Uses the `argparse` package to accept user inputs including count data file, metadata file, control condition, significance thresholds, fold-change cutoff, and output directory.

3. **Data Reading and Pre-processing**  
   - Reads a counts matrix (genes × samples) and a metadata file (samples as rows, must include a `condition` column).  
   - Validates that the specified control condition exists.  
   - Filters out low-count genes to improve statistical power.

4. **Differential Expression Analysis (DESeq2)**  
   - Constructs a DESeq2 object with a `~condition` design.  
   - Normalizes the data, estimates dispersion, tests for differential expression, and applies adaptive shrinkage (using the `ashr` method) to reduce noise.
   - Exports the full DESeq2 object, normalized counts, and differential expression results for each condition vs. the control.

5. **Pathway Enrichment Analysis (ORA)**  
   - Converts ENSEMBL IDs to ENTREZ IDs for gene enrichment functions.  
   - Performs enrichment analysis for GO Biological Process, KEGG, and Reactome pathways on significant genes.  
   - Generates and exports dotplot visualizations for each enrichment analysis.

6. **Gene Set Enrichment Analysis (GSEA)**  
   - Ranks all genes by their log2 fold change and deduplicates based on ENTREZ IDs.  
   - Runs GSEA on GO, KEGG, and Reactome gene sets.  
   - Saves the GSEA results and generates dotplot visualizations.

7. **Output**  
   All outputs, including R objects, CSV files, and PNG plots, are saved in the specified output directory.

8. **Downstream Analysis**

  The DESeq2 object (`dds`) produced by this script is used to generate an interactive HTML report. This report includes run details and key differential expression  plots such as MA, volcano, PCA, and heatmap visualizations. This can be seen in the https://github.com/morgannamie/Differential-Expression-Reports repo


## Prerequisites

Make sure the following R packages are installed:

- **argparse**
- **DESeq2**
- **clusterProfiler**
- **org.Hs.eg.db**
- **ReactomePA**
- **enrichplot**
- **ggplot2**
- **ashr**

You can install them as follows:

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("DESeq2", "clusterProfiler", "org.Hs.eg.db", "ReactomePA", "enrichplot"))
install.packages(c("argparse", "ggplot2", "ashr"))
