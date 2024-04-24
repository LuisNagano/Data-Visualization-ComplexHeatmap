# Load necessary libraries
library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)
library(colorRamps)
library(biomaRt)
library(edgeR)

# Define base path for files
setwd("C:/Users/bruep/OneDrive/Área de Trabalho/Kleiton/Heatmaps/Mouse")
base_path <- "C:/Users/bruep/OneDrive/Área de Trabalho/Kleiton/Heatmaps/Mouse"

# Load data
file_path <- file.path(base_path, "rawcounts_RNAseq_Mouse.csv")
counts <- read.csv(file_path, row.names = 1)

# Create DGEList object
y <- DGEList(counts = counts)

# Define groups for each sample based on column names
group_vector <- sub("(\\D+).*", "\\1", colnames(y))
group_vector <- gsub("WT", "WT", group_vector) # Wild Type
group_vector <- gsub("BPCre", "BPCre", group_vector) # BPCre Group
group_vector <- gsub("BCre", "BCre", group_vector) # BCre Group
group_vector <- gsub("PCre", "PCre", group_vector) # PCre Group
group_vector <- gsub("ACC", "ACC", group_vector) # ACC Group

# Update y object with group information
y$samples$group <- factor(group_vector)

# Filter genes with low expression, now considering groups
keep <- filterByExpr(y, group = y$samples$group)
y <- y[keep, , keep.lib.sizes = FALSE]

# Calculate normalization factors using TMM method
y <- calcNormFactors(y)

# Estimate dispersions
y <- estimateDisp(y)

# Set colors for each group
colors <- c("red", "blue", "green", "purple", "orange")
names(colors) <- levels(y$samples$group)

# Set plot size and create plot with more space on the right for legend
png(file.path(base_path, "MDS_Plot_Normalized_Samples.png"), width = 800, height = 600)
par(mar = c(5.1, 4.9, 4.1, 12.1)) # Adjust margins (bottom, left, top, right)
plotMDS(y, col = colors[y$samples$group], pch = 20, cex.lab = 1.5, cex.axis = 1.5, cex = 3, main = "MDS Plot of Normalized Samples", cex.main = 2)
# Draw legend outside the main plot, more to the right
legend("topright", inset = c(-0.2, 0.3), legend = levels(y$samples$group), fill = colors, pch = 20, cex = 1.5, pt.cex = 3, bty = "n", xpd = TRUE)
dev.off()

# Print normalization factors
print(y$samples$norm.factors)

# Compute CPMs for further analysis
y <- cpm(y, log = TRUE)

# Function to scale data by row
scale_rows <- function(x) {
  m <- apply(x, 1, mean, na.rm = TRUE)
  s <- apply(x, 1, sd, na.rm = TRUE)
  return((x - m) / s)
}

# Load gene lists
genes_of_interest_1 <- tolower(read.table(file.path(base_path, "IFNG.GS.txt"), header = FALSE, stringsAsFactors = FALSE)$V1)
genes_of_interest_2 <- tolower(read.csv(file.path(base_path, "ISG_RS_genes.csv"), header = FALSE, stringsAsFactors = FALSE)$V1)
genes_of_interest_3 <- tolower(read.csv(file.path(base_path, "list_720_epifactors.csv"), header = FALSE, stringsAsFactors = FALSE)$V1)

# Load gene expression data
gene_expression_data <- y
gene_expression_data <- gene_expression_data[, order(colnames(gene_expression_data))] # order columns alphabetically
colnames(gene_expression_data) <- tolower(colnames(gene_expression_data)) # Convert to lowercase

# Map Ensembl IDs to Gene Symbols
ensembl <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
ensembl_ids <- rownames(gene_expression_data)
gene_symbols <- getBM(attributes = c("ensembl_gene_id", "external_gene_name"), 
                      filters = "ensembl_gene_id", 
                      values = ensembl_ids, 
                      mart = ensembl)

gene_symbols$external_gene_name <- tolower(gene_symbols$external_gene_name)

# Replace Ensembl IDs with Gene Symbols in gene lists of interest
genes_of_interest_1 <- gene_symbols[gene_symbols$external_gene_name %in% genes_of_interest_1, "ensembl_gene_id"]
genes_of_interest_2 <- gene_symbols[gene_symbols$external_gene_name %in% genes_of_interest_2, "ensembl_gene_id"]
genes_of_interest_3 <- gene_symbols[gene_symbols$external_gene_name %in% genes_of_interest_3, "ensembl_gene_id"]

# Filter data to include only genes of interest
filtered_data_1 <- gene_expression_data[rownames(gene_expression_data) %in% genes_of_interest_1, ]
filtered_data_2 <- gene_expression_data[rownames(gene_expression_data) %in% genes_of_interest_2, ]
filtered_data_3 <- gene_expression_data[rownames(gene_expression_data) %in% genes_of_interest_3, ]
filtered_data_3_clean <- filtered_data_3[apply(filtered_data_3, 1, function(row) length(unique(row)) > 1), ]

# Transpose and scale data
tab_heatmap_1 <- filtered_data_1
scale_tab_1 <- scale_rows(as.matrix(tab_heatmap_1))
tab_heatmap_2 <- filtered_data_2
scale_tab_2 <- scale_rows(as.matrix(tab_heatmap_2))
tab_heatmap_3 <- filtered_data_3_clean
scale_tab_3 <- scale_rows(as.matrix(tab_heatmap_3))

# Load sample conditions
sample_info <- read.delim(paste0(base_path, "/phenotype.txt"), row.names = 1)
sampleCondition <- sample_info$Phenotype

# Order gene expression matrix based on sample condition
ordered_indices <- order(sampleCondition)
gene_expression_data <- gene_expression_data[, ordered_indices]
sampleCondition <- sampleCondition[ordered_indices]

# Update conditions to specific names
sampleCondition <- gsub("wt", "WT", sampleCondition)
sampleCondition <- gsub("bpcre", "BPCre", sampleCondition)
sampleCondition <- gsub("bcre", "BCre", sampleCondition)
sampleCondition <- gsub("pcre", "PCre", sampleCondition)
sampleCondition <- gsub("acc", "ACC", sampleCondition)

# Adjust color mapping to include all found subgroups
colors_for_subgroups <- c("WT" = "red", 
                          "BPCre" = "blue", 
                          "BCre" = "yellow", 
                          "PCre" = "black",
                          "ACC" = "orange")

# Create dataframe for heatmap annotation
df1 <- data.frame(subgroups = sampleCondition, row.names = colnames(tab_heatmap_1))
ha1 <- HeatmapAnnotation(df = df1, col = list(subgroups = colors_for_subgroups))

df2 <- data.frame(subgroups = sampleCondition, row.names = colnames(tab_heatmap_2))
ha2 <- HeatmapAnnotation(df = df2, col = list(subgroups = colors_for_subgroups))

df3 <- data.frame(subgroups = sampleCondition, row.names = colnames(tab_heatmap_3))
ha3 <- HeatmapAnnotation(df = df3, col = list(subgroups = colors_for_subgroups))

# Define colors and other options for the heatmap
breaks <- seq(-2, 2, by = 0.1)

# Generate heatmaps for genes_of_interest_1
ht1_1 <- Heatmap(scale_tab_1[, sampleCondition %in% c("WT", "BPCre")],
                 name = "zscore", column_title = "WT vs BPCre - IFNG.GS genes",
                 row_names_gp = gpar(fontsize = 5),
                 top_annotation = ha1[sampleCondition %in% c("WT", "BPCre"), ],
                 show_row_names = F, show_column_names = F,
                 cluster_rows = T, cluster_columns = F,
                 col = colorRamp2(breaks, colorRampPalette(rev(brewer.pal(n = 10, name = "RdBu")))(41)),
                 show_column_dend = T, show_row_dend = T, 
                 clustering_distance_columns = "spearman", clustering_method_columns = "ward.D2", 
                 clustering_distance_rows = "spearman", clustering_method_rows = "ward.D2")

ht1_2 <- Heatmap(scale_tab_1[, sampleCondition %in% c("WT", "BCre")],
                 name = "zscore", column_title = "WT vs BCre - IFNG.GS genes",
                 row_names_gp = gpar(fontsize = 5),
                 top_annotation = ha1[sampleCondition %in% c("WT", "BCre"), ],
                 show_row_names = F, show_column_names = F,
                 cluster_rows = T, cluster_columns = F,
                 col = colorRamp2(breaks, colorRampPalette(rev(brewer.pal(n = 10, name = "RdBu")))(41)),
                 show_column_dend = T, show_row_dend = T, 
                 clustering_distance_columns = "spearman", clustering_method_columns = "ward.D2", 
                 clustering_distance_rows = "spearman", clustering_method_rows = "ward.D2")

ht1_3 <- Heatmap(scale_tab_1[, sampleCondition %in% c("WT", "PCre")],
                 name = "zscore", column_title = "WT vs PCre - IFNG.GS genes", 
                 row_names_gp = gpar(fontsize = 5),
                 top_annotation = ha1[sampleCondition %in% c("WT", "PCre"), ],
                 show_row_names = F, show_column_names = F,
                 cluster_rows = T, cluster_columns = F,
                 col = colorRamp2(breaks, colorRampPalette(rev(brewer.pal(n = 10, name = "RdBu")))(41)),
                 show_column_dend = T, show_row_dend = T, 
                 clustering_distance_columns = "spearman", clustering_method_columns = "ward.D2", 
                 clustering_distance_rows = "spearman", clustering_method_rows = "ward.D2")

# Generate heatmaps for genes_of_interest_2
ht2_1 <- Heatmap(scale_tab_2[, sampleCondition %in% c("WT", "BPCre")],
                 name = "zscore", column_title = "WT vs BPCre - ISG_RS genes",
                 row_names_gp = gpar(fontsize = 5),
                 top_annotation = ha2[sampleCondition %in% c("WT", "BPCre"), ],
                 show_row_names = F, show_column_names = F,
                 cluster_rows = T, cluster_columns = F,
                 col = colorRamp2(breaks, colorRampPalette(rev(brewer.pal(n = 10, name = "RdBu")))(41)),
                 show_column_dend = T, show_row_dend = T, 
                 clustering_distance_columns = "spearman", clustering_method_columns = "ward.D2", 
                 clustering_distance_rows = "spearman", clustering_method_rows = "ward.D2")

ht2_2 <- Heatmap(scale_tab_2[, sampleCondition %in% c("WT", "BCre")],
                 name = "zscore", column_title = "WT vs BCre - ISG_RS genes",
                 row_names_gp = gpar(fontsize = 5),
                 top_annotation = ha2[sampleCondition %in% c("WT", "BCre"), ],
                 show_row_names = F, show_column_names = F,
                 cluster_rows = T, cluster_columns = F,
                 col = colorRamp2(breaks, colorRampPalette(rev(brewer.pal(n = 10, name = "RdBu")))(41)),
                 show_column_dend = T, show_row_dend = T, 
                 clustering_distance_columns = "spearman", clustering_method_columns = "ward.D2", 
                 clustering_distance_rows = "spearman", clustering_method_rows = "ward.D2")

ht2_3 <- Heatmap(scale_tab_2[, sampleCondition %in% c("WT", "PCre")],
                 name = "zscore", column_title = "WT vs PCre - ISG_RS genes", 
                 row_names_gp = gpar(fontsize = 5),
                 top_annotation = ha2[sampleCondition %in% c("WT", "PCre"), ],
                 show_row_names = F, show_column_names = F,
                 cluster_rows = T, cluster_columns = F,
                 col = colorRamp2(breaks, colorRampPalette(rev(brewer.pal(n = 10, name = "RdBu")))(41)),
                 show_column_dend = T, show_row_dend = T, 
                 clustering_distance_columns = "spearman", clustering_method_columns = "ward.D2", 
                 clustering_distance_rows = "spearman", clustering_method_rows = "ward.D2")

# Generate heatmaps for genes_of_interest_3
ht3_1 <- Heatmap(scale_tab_3[, sampleCondition %in% c("WT", "BPCre")],
                 name = "zscore", column_title = "WT vs BPCre - list_720_epifactors",
                 row_names_gp = gpar(fontsize = 5),
                 top_annotation = ha3[sampleCondition %in% c("WT", "BPCre"), ],
                 show_row_names = F, show_column_names = F,
                 cluster_rows = T, cluster_columns = F,
                 col = colorRamp2(breaks, colorRampPalette(rev(brewer.pal(n = 10, name = "RdBu")))(41)),
                 show_column_dend = T, show_row_dend = T, 
                 clustering_distance_columns = "spearman", clustering_method_columns = "ward.D2", 
                 clustering_distance_rows = "spearman", clustering_method_rows = "ward.D2")

ht3_2 <- Heatmap(scale_tab_3[, sampleCondition %in% c("WT", "BCre")],
                 name = "zscore", column_title = "WT vs BCre - list_720_epifactors",
                 row_names_gp = gpar(fontsize = 5),
                 top_annotation = ha3[sampleCondition %in% c("WT", "BCre"), ],
                 show_row_names = F, show_column_names = F,
                 cluster_rows = T, cluster_columns = F,
                 col = colorRamp2(breaks, colorRampPalette(rev(brewer.pal(n = 10, name = "RdBu")))(41)),
                 show_column_dend = T, show_row_dend = T, 
                 clustering_distance_columns = "spearman", clustering_method_columns = "ward.D2", 
                 clustering_distance_rows = "spearman", clustering_method_rows = "ward.D2")

ht3_3 <- Heatmap(scale_tab_3[, sampleCondition %in% c("WT", "PCre")],
                 name = "zscore", column_title = "WT vs PCre - list_720_epifactors", 
                 row_names_gp = gpar(fontsize = 5),
                 top_annotation = ha3[sampleCondition %in% c("WT", "PCre"), ],
                 show_row_names = F, show_column_names = F,
                 cluster_rows = T, cluster_columns = F,
                 col = colorRamp2(breaks, colorRampPalette(rev(brewer.pal(n = 10, name = "RdBu")))(41)),
                 show_column_dend = T, show_row_dend = T, 
                 clustering_distance_columns = "spearman", clustering_method_columns = "ward.D2", 
                 clustering_distance_rows = "spearman", clustering_method_rows = "ward.D2")

# Export the heatmaps as PNG
png(filename = paste0(base_path, "/ht1_1.png"), width = 9, height = 7, units = "in", res = 500)
print(ht1_1)
dev.off()

png(filename = paste0(base_path, "/ht1_2.png"), width = 9, height = 7, units = "in", res = 500)
print(ht1_2)
dev.off()

png(filename = paste0(base_path, "/ht1_3.png"), width = 9, height = 7, units = "in", res = 500)
print(ht1_3)
dev.off()

png(filename = paste0(base_path, "/ht2_1.png"), width = 9, height = 7, units = "in", res = 500)
print(ht2_1)
dev.off()

png(filename = paste0(base_path, "/ht2_2.png"), width = 9, height = 7, units = "in", res = 500)
print(ht2_2)
dev.off()

png(filename = paste0(base_path, "/ht2_3.png"), width = 9, height = 7, units = "in", res = 500)
print(ht2_3)
dev.off()

png(filename = paste0(base_path, "/ht3_1.png"), width = 9, height = 7, units = "in", res = 500)
print(ht3_1)
dev.off()

png(filename = paste0(base_path, "/ht3_2.png"), width = 9, height = 7, units = "in", res = 500)
print(ht3_2)
dev.off()

png(filename = paste0(base_path, "/ht3_3.png"), width = 9, height = 7, units = "in", res = 500)
print(ht3_3)
dev.off()

