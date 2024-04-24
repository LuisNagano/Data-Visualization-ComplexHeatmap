# Load necessary libraries
library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)
library(colorRamps)

# Function to scale data by row
scale_rows <- function(x) {
  m <- apply(x, 1, mean, na.rm = TRUE)
  s <- apply(x, 1, sd, na.rm = TRUE)
  return((x - m) / s)
}

# Define the base path for files
base_path <- "C:/Users/bruep/OneDrive/Área de Trabalho/Kleiton/Heatmaps"

# Load gene list 1
genes_of_interest_1 <- read.table(paste0(base_path, "/IFNG.GS.txt"), header = FALSE, stringsAsFactors = FALSE)
genes_of_interest_1 <- genes_of_interest_1$V1

# Load gene list 2
genes_of_interest_2 <- read.csv(paste0(base_path, "/ISG_RS_genes.csv"), header = FALSE, stringsAsFactors = FALSE)
genes_of_interest_2 <- genes_of_interest_2$V1

# Load gene list 3
genes_of_interest_3 <- read.csv(paste0(base_path, "/list_720_epifactors.csv"), header = FALSE, stringsAsFactors = FALSE)
genes_of_interest_3 <- genes_of_interest_3$V1

# Load gene expression data
gene_expression_data <- read.delim(paste0(base_path, "/vsd_TCGA_GTEX_Adrenal_COCs_samples.txt"), row.names = 1)

# Filter data to include only genes of interest 1
filtered_data_1 <- gene_expression_data[rownames(gene_expression_data) %in% genes_of_interest_1, ]

# Filter data to include only genes of interest 2
filtered_data_2 <- gene_expression_data[rownames(gene_expression_data) %in% genes_of_interest_2, ]

# Filter data to include only genes of interest 3 and remove rows with equal values or no variance
filtered_data_3 <- gene_expression_data[rownames(gene_expression_data) %in% genes_of_interest_3, ]
filtered_data_3_clean <- filtered_data_3[apply(filtered_data_3, 1, function(row) length(unique(row)) > 1), ]

# Transpose and scale the data
tab_heatmap_1 <- filtered_data_1
scale_tab_1 <- scale_rows(as.matrix(tab_heatmap_1))

tab_heatmap_2 <- filtered_data_2
scale_tab_2 <- scale_rows(as.matrix(tab_heatmap_2))

tab_heatmap_3 <- filtered_data_3_clean
scale_tab_3 <- scale_rows(as.matrix(tab_heatmap_3))

# Load sample conditions
sample_info <- read.delim(paste0(base_path, "/phenotype.txt"), row.names = 1)
sampleCondition <- sample_info$Phenotype
sampleCondition <- gsub("Normal ", "Normal", sampleCondition) # Remove extra space

# Update conditions to specific names
sampleCondition <- gsub("COC1", "ACC-COC1", sampleCondition)
sampleCondition <- gsub("COC2", "ACC-COC2", sampleCondition)
sampleCondition <- gsub("COC3", "ACC-COC3", sampleCondition)
sampleCondition <- gsub("Normal", "Normal Adrenocortical Gland", sampleCondition)

# Adjust color mapping to include all found subgroups
colors_for_subgroups <- c("ACC-COC1" = "red", 
                          "ACC-COC2" = "blue", 
                          "ACC-COC3" = "yellow", 
                          "Normal Adrenocortical Gland" = "black")

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
ht1_1 <- Heatmap(scale_tab_1[, sampleCondition %in% c("Normal Adrenocortical Gland", "ACC-COC1")],
                 name = "zscore", column_title = "Normal vs ACC-COC1 - IFNG.GS genes",
                 row_names_gp = gpar(fontsize = 5),
                 top_annotation = ha1[sampleCondition %in% c("Normal Adrenocortical Gland", "ACC-COC1"), ],
                 show_row_names = F, show_column_names = F,
                 cluster_rows = T, cluster_columns = F,
                 col = colorRamp2(breaks, colorRampPalette(rev(brewer.pal(n = 10, name = "RdBu")))(41)),
                 show_column_dend = T, show_row_dend = T, 
                 clustering_distance_columns = "spearman", clustering_method_columns = "ward.D2", 
                 clustering_distance_rows = "spearman", clustering_method_rows = "ward.D2")

ht1_2 <- Heatmap(scale_tab_1[, sampleCondition %in% c("Normal Adrenocortical Gland", "ACC-COC2")],
                 name = "zscore", column_title = "Normal vs ACC-COC2 - IFNG.GS genes",
                 row_names_gp = gpar(fontsize = 5),
                 top_annotation = ha1[sampleCondition %in% c("Normal Adrenocortical Gland", "ACC-COC2"), ],
                 show_row_names = F, show_column_names = F,
                 cluster_rows = T, cluster_columns = F,
                 col = colorRamp2(breaks, colorRampPalette(rev(brewer.pal(n = 10, name = "RdBu")))(41)),
                 show_column_dend = T, show_row_dend = T, 
                 clustering_distance_columns = "spearman", clustering_method_columns = "ward.D2", 
                 clustering_distance_rows = "spearman", clustering_method_rows = "ward.D2")

ht1_3 <- Heatmap(scale_tab_1[, sampleCondition %in% c("Normal Adrenocortical Gland", "ACC-COC3")],
                 name = "zscore", column_title = "Normal vs ACC-COC3 - IFNG.GS genes", 
                 row_names_gp = gpar(fontsize = 5),
                 top_annotation = ha1[sampleCondition %in% c("Normal Adrenocortical Gland", "ACC-COC3"), ],
                 show_row_names = F, show_column_names = F,
                 cluster_rows = T, cluster_columns = F,
                 col = colorRamp2(breaks, colorRampPalette(rev(brewer.pal(n = 10, name = "RdBu")))(41)),
                 show_column_dend = T, show_row_dend = T, 
                 clustering_distance_columns = "spearman", clustering_method_columns = "ward.D2", 
                 clustering_distance_rows = "spearman", clustering_method_rows = "ward.D2")

ht1_4 <- Heatmap(scale_tab_1[, sampleCondition %in% c("ACC-COC1", "ACC-COC2", "ACC-COC3")],
                 name = "zscore", column_title = "COC1 vs COC2 vs COC3 - IFNG.GS genes",
                 row_names_gp = gpar(fontsize = 5),
                 top_annotation = ha1[sampleCondition %in% c("ACC-COC1", "ACC-COC2", "ACC-COC3"), ],
                 show_row_names = F, show_column_names = F,
                 cluster_rows = T, cluster_columns = F,
                 col = colorRamp2(breaks, colorRampPalette(rev(brewer.pal(n = 10, name = "RdBu")))(41)),
                 show_column_dend = T, show_row_dend = T, 
                 clustering_distance_columns = "spearman", clustering_method_columns = "ward.D2", 
                 clustering_distance_rows = "spearman", clustering_method_rows = "ward.D2")

ht1_5 <- Heatmap(scale_tab_1[, sampleCondition %in% c("ACC-COC1", "ACC-COC2", "ACC-COC3", "Normal Adrenocortical Gland")],
                 name = "zscore", column_title = "COC1 vs COC2 vs COC3 vs Normal - IFNG.GS genes",
                 row_names_gp = gpar(fontsize = 5),
                 top_annotation = ha1[sampleCondition %in% c("ACC-COC1", "ACC-COC2", "ACC-COC3", "Normal Adrenocortical Gland"), ],
                 show_row_names = F, show_column_names = F,
                 cluster_rows = T, cluster_columns = F,
                 col = colorRamp2(breaks, colorRampPalette(rev(brewer.pal(n = 10, name = "RdBu")))(41)),
                 show_column_dend = T, show_row_dend = T, 
                 clustering_distance_columns = "spearman", clustering_method_columns = "ward.D2", 
                 clustering_distance_rows = "spearman", clustering_method_rows = "ward.D2")

# Generate heatmaps for genes_of_interest_2
ht2_1 <- Heatmap(scale_tab_2[, sampleCondition %in% c("Normal Adrenocortical Gland", "ACC-COC1")],
                 name = "zscore", column_title = "Normal vs ACC-COC1 - ISG_RS genes",
                 row_names_gp = gpar(fontsize = 5),
                 top_annotation = ha2[sampleCondition %in% c("Normal Adrenocortical Gland", "ACC-COC1"), ],
                 show_row_names = F, show_column_names = F,
                 cluster_rows = T, cluster_columns = F,
                 col = colorRamp2(breaks, colorRampPalette(rev(brewer.pal(n = 10, name = "RdBu")))(41)),
                 show_column_dend = T, show_row_dend = T, 
                 clustering_distance_columns = "spearman", clustering_method_columns = "ward.D2", 
                 clustering_distance_rows = "spearman", clustering_method_rows = "ward.D2")

ht2_2 <- Heatmap(scale_tab_2[, sampleCondition %in% c("Normal Adrenocortical Gland", "ACC-COC2")],
                 name = "zscore", column_title = "Normal vs ACC-COC2 - ISG_RS genes",
                 row_names_gp = gpar(fontsize = 5),
                 top_annotation = ha2[sampleCondition %in% c("Normal Adrenocortical Gland", "ACC-COC2"), ],
                 show_row_names = F, show_column_names = F,
                 cluster_rows = T, cluster_columns = F,
                 col = colorRamp2(breaks, colorRampPalette(rev(brewer.pal(n = 10, name = "RdBu")))(41)),
                 show_column_dend = T, show_row_dend = T, 
                 clustering_distance_columns = "spearman", clustering_method_columns = "ward.D2", 
                 clustering_distance_rows = "spearman", clustering_method_rows = "ward.D2")

ht2_3 <- Heatmap(scale_tab_2[, sampleCondition %in% c("Normal Adrenocortical Gland", "ACC-COC3")],
                 name = "zscore", column_title = "Normal vs ACC-COC3 - ISG_RS genes", 
                 row_names_gp = gpar(fontsize = 5),
                 top_annotation = ha2[sampleCondition %in% c("Normal Adrenocortical Gland", "ACC-COC3"), ],
                 show_row_names = F, show_column_names = F,
                 cluster_rows = T, cluster_columns = F,
                 col = colorRamp2(breaks, colorRampPalette(rev(brewer.pal(n = 10, name = "RdBu")))(41)),
                 show_column_dend = T, show_row_dend = T, 
                 clustering_distance_columns = "spearman", clustering_method_columns = "ward.D2", 
                 clustering_distance_rows = "spearman", clustering_method_rows = "ward.D2")

ht2_4 <- Heatmap(scale_tab_2[, sampleCondition %in% c("ACC-COC1", "ACC-COC2", "ACC-COC3")],
                 name = "zscore", column_title = "COC1 vs COC2 vs COC3 - ISG_RS genes",
                 row_names_gp = gpar(fontsize = 5),
                 top_annotation = ha2[sampleCondition %in% c("ACC-COC1", "ACC-COC2", "ACC-COC3"), ],
                 show_row_names = F, show_column_names = F,
                 cluster_rows = T, cluster_columns = F,
                 col = colorRamp2(breaks, colorRampPalette(rev(brewer.pal(n = 10, name = "RdBu")))(41)),
                 show_column_dend = T, show_row_dend = T, 
                 clustering_distance_columns = "spearman", clustering_method_columns = "ward.D2", 
                 clustering_distance_rows = "spearman", clustering_method_rows = "ward.D2")

ht2_5 <- Heatmap(scale_tab_2[, sampleCondition %in% c("ACC-COC1", "ACC-COC2", "ACC-COC3", "Normal Adrenocortical Gland")],
                 name = "zscore", column_title = "COC1 vs COC2 vs COC3 vs Normal - ISG_RS genes",
                 row_names_gp = gpar(fontsize = 5),
                 top_annotation = ha2[sampleCondition %in% c("ACC-COC1", "ACC-COC2", "ACC-COC3", "Normal Adrenocortical Gland"), ],
                 show_row_names = F, show_column_names = F,
                 cluster_rows = T, cluster_columns = F,
                 col = colorRamp2(breaks, colorRampPalette(rev(brewer.pal(n = 10, name = "RdBu")))(41)),
                 show_column_dend = T, show_row_dend = T, 
                 clustering_distance_columns = "spearman", clustering_method_columns = "ward.D2", 
                 clustering_distance_rows = "spearman", clustering_method_rows = "ward.D2")

# Generate heatmaps for genes_of_interest_3
ht3_1 <- Heatmap(scale_tab_3[, sampleCondition %in% c("Normal Adrenocortical Gland", "ACC-COC1")],
                 name = "zscore", column_title = "Normal vs ACC-COC1 - list_720_epifactors",
                 row_names_gp = gpar(fontsize = 5),
                 top_annotation = ha3[sampleCondition %in% c("Normal Adrenocortical Gland", "ACC-COC1"), ],
                 show_row_names = F, show_column_names = F,
                 cluster_rows = T, cluster_columns = F,
                 col = colorRamp2(breaks, colorRampPalette(rev(brewer.pal(n = 10, name = "RdBu")))(41)),
                 show_column_dend = T, show_row_dend = T, 
                 clustering_distance_columns = "spearman", clustering_method_columns = "ward.D2", 
                 clustering_distance_rows = "spearman", clustering_method_rows = "ward.D2")

ht3_2 <- Heatmap(scale_tab_3[, sampleCondition %in% c("Normal Adrenocortical Gland", "ACC-COC2")],
                 name = "zscore", column_title = "Normal vs ACC-COC2 - list_720_epifactors",
                 row_names_gp = gpar(fontsize = 5),
                 top_annotation = ha3[sampleCondition %in% c("Normal Adrenocortical Gland", "ACC-COC2"), ],
                 show_row_names = F, show_column_names = F,
                 cluster_rows = T, cluster_columns = F,
                 col = colorRamp2(breaks, colorRampPalette(rev(brewer.pal(n = 10, name = "RdBu")))(41)),
                 show_column_dend = T, show_row_dend = T, 
                 clustering_distance_columns = "spearman", clustering_method_columns = "ward.D2", 
                 clustering_distance_rows = "spearman", clustering_method_rows = "ward.D2")

ht3_3 <- Heatmap(scale_tab_3[, sampleCondition %in% c("Normal Adrenocortical Gland", "ACC-COC3")],
                 name = "zscore", column_title = "Normal vs ACC-COC3 - list_720_epifactors", 
                 row_names_gp = gpar(fontsize = 5),
                 top_annotation = ha3[sampleCondition %in% c("Normal Adrenocortical Gland", "ACC-COC3"), ],
                 show_row_names = F, show_column_names = F,
                 cluster_rows = T, cluster_columns = F,
                 col = colorRamp2(breaks, colorRampPalette(rev(brewer.pal(n = 10, name = "RdBu")))(41)),
                 show_column_dend = T, show_row_dend = T, 
                 clustering_distance_columns = "spearman", clustering_method_columns = "ward.D2", 
                 clustering_distance_rows = "spearman", clustering_method_rows = "ward.D2")

# Remove columns with zero variance
scale_tab_3_filtered <- scale_tab_3[, apply(scale_tab_3, 2, var, na.rm = TRUE) > 0]

ht3_4 <- Heatmap(scale_tab_3_filtered[, sampleCondition %in% c("ACC-COC1", "ACC-COC2", "ACC-COC3")],
                 name = "zscore", column_title = "COC1 vs COC2 vs COC3 - list_720_epifactors",
                 row_names_gp = gpar(fontsize = 5),
                 top_annotation = ha3[sampleCondition %in% c("ACC-COC1", "ACC-COC2", "ACC-COC3"), ],
                 show_row_names = F, show_column_names = F,
                 cluster_rows = T, cluster_columns = F,
                 col = colorRamp2(breaks, colorRampPalette(rev(brewer.pal(n = 10, name = "RdBu")))(41)),
                 show_column_dend = T, show_row_dend = T, 
                 clustering_distance_columns = "euclidean", clustering_method_columns = "ward.D2", 
                 clustering_distance_rows = "euclidean", clustering_method_rows = "ward.D2")

ht3_5 <- Heatmap(scale_tab_3_filtered[, sampleCondition %in% c("ACC-COC1", "ACC-COC2", "ACC-COC3", "Normal Adrenocortical Gland")],
                 name = "zscore", column_title = "COC1 vs COC2 vs COC3 vs Normal - list_720_epifactors",
                 row_names_gp = gpar(fontsize = 5),
                 top_annotation = ha3[sampleCondition %in% c("ACC-COC1", "ACC-COC2", "ACC-COC3", "Normal Adrenocortical Gland"), ],
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

png(filename = paste0(base_path, "/ht1_4.png"), width = 9, height = 7, units = "in", res = 500)
print(ht1_4)
dev.off()

png(filename = paste0(base_path, "/ht1_5.png"), width = 9, height = 7, units = "in", res = 500)
print(ht1_5)
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

png(filename = paste0(base_path, "/ht2_4.png"), width = 9, height = 7, units = "in", res = 500)
print(ht2_4)
dev.off()

png(filename = paste0(base_path, "/ht2_5.png"), width = 9, height = 7, units = "in", res = 500)
print(ht2_5)
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

png(filename = paste0(base_path, "/ht3_4.png"), width = 9, height = 7, units = "in", res = 500)
print(ht3_4)
dev.off()

png(filename = paste0(base_path, "/ht3_5.png"), width = 9, height = 7, units = "in", res = 500)
print(ht3_5)
dev.off()


