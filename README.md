# Data Visualization using ComplexHeatmap package

## RNA-Seq Data Analysis for Mouse Samples

This project focuses on the analysis of RNA-seq data obtained from mouse samples. The main objectives are to process and normalize the gene expression data, select genes of interest based on predefined lists, and generate heatmaps to visualize the expression patterns across different experimental conditions.

## Overview

The analysis begins by loading the necessary libraries and setting up the working environment. The gene expression data is read from a CSV file and converted into a DGEList object using the edgeR package. The samples are grouped based on their experimental conditions (WT, BPCre, BCre, PCre, and ACC).

After filtering out genes with low expression levels, the data is normalized using the TMM (Trimmed Mean of M-values) method. Dispersions are estimated, and normalization factors are calculated.

## Gene Selection

Three gene lists of interest are loaded from text files:

1. `IFNG.GS.txt`: A list of genes related to the Interferon-gamma (IFNG) gene set.
2. `ISG_RS_genes.csv`: A list of Interferon-stimulated genes (ISGs) and related genes.
3. `list_720_epifactors.csv`: A list of 720 epifactors (epigenetic factors).

The Ensembl IDs in the gene expression data are mapped to gene symbols using BioMart. The gene lists are then filtered to include only the genes present in the expression data.

## Heatmap Generation

The gene expression data is scaled by row and transposed for heatmap visualization. Sample information, including the experimental condition for each sample, is loaded from a separate file.

Heatmaps are generated using the ComplexHeatmap package for each gene list, comparing the expression patterns between the WT (wild-type) condition and other experimental conditions (BPCre, BCre, and PCre).

The heatmaps are color-coded based on the z-score values, with blue representing lower expression and red representing higher expression. Clustering is performed on both rows (genes) and columns (samples) using Spearman correlation distance and Ward's linkage method.

The resulting heatmaps are exported as PNG images with filenames indicating the gene list and experimental conditions compared.

## Results

The generated heatmaps provide a visual representation of the gene expression patterns for the selected gene lists across different experimental conditions. These heatmaps can be used for further analysis and interpretation of the RNA-seq data.

## Usage

To reproduce the analysis, ensure that you have the required input files (gene expression data, gene lists, and sample information) in the appropriate locations specified in the code. Additionally, make sure that all the necessary R packages are installed.

## Dependencies

The following R packages are required for this analysis:

- ComplexHeatmap
- RColorBrewer
- circlize
- colorRamps
- biomaRt
- edgeR
