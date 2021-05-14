# Analysis of scRNA-seq human melanoma cell line data from Tirosh et al., 2016 https://science.sciencemag.org/content/352/6282/189/tab-pdf
# Data from GEO: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE72056

library(Seurat)
library(tidyverse)
library(scales)
library(pals)
library(cowplot)
library(fgsea)
library(msigdbr)
library(data.table)
library(mvhspatialplots)


# Load counts matrix and create Seurat object:
counts <- read.delim("/Users/hunterm/Dropbox/10x/Tirosh_data/GSE72056_melanoma_single_cell_revised_v2.txt")
counts_cells_only <- counts[4:nrow(counts),] 
rownames(counts_cells_only) <- NULL

# Remove MARCH1 and MARCH2 because for some reason they're duplicated in the counts matrix
counts_cells_only <- counts_cells_only[(counts_cells_only$Cell != c("MARCH1", "MARCH2")),] 
rownames(counts_cells_only) <- counts_cells_only$Cell
counts_cells_only <- counts_cells_only %>% .[,2:ncol(.)]

# Create Seurat object
Tirosh_data <- CreateSeuratObject(counts = counts_cells_only,
                                  min.cells = 3,
                                  min.features = 200,
                                  project = "Tirosh")


# Add metadata from authors to Seurat object
metadata <- counts[1:3,] %>% t() %>% data.frame() %>% rownames_to_column(var = "cell")
colnames(metadata) <- c("cell", "tumor", "is_malignant", "non_malignant_cell_type")

# Is tumor?
tumor_metadata <- AddMetadata(Tirosh_data,
                              metadata = metadata %>% dplyr::select(cell, tumor) %>% .[2:nrow(.),] %>% deframe(),
                              col.name = "tumor")

# Is malignant?
malignant_metadata <- metadata %>% dplyr::select(cell, is_malignant) %>% .[2:nrow(.),] 
malignant_metadata$is_malignant <- malignant_metadata$is_malignant %>% as.character()
malignant_metadata$is_malignant <- gsub(pattern = " 1", replacement = "no", x = malignant_metadata$is_malignant) %>%
  gsub(pattern = " 2", replacement = "yes", x = .) %>%
  gsub(pattern = " 0", replacement = "unresolved", x = .)

Tirosh_data <- AddMetaData(Tirosh_data,
                           metadata = malignant_metadata %>% deframe(),
                           col.name = "is_malignant")

# Non malignant cell types
celltype_metadata <- metadata %>% dplyr::select(cell, non_malignant_cell_type) %>% .[2:nrow(.),]
celltype_metadata$non_malignant_cell_type <- gsub(pattern = " 1", replacement = "T cell", x = celltype_metadata$non_malignant_cell_type) %>%
  gsub(pattern = " 2", replacement = "B cell", x = .) %>%
  gsub(pattern = " 3", replacement = "macrophage", x = .) %>%
  gsub(pattern = " 4", replacement = "endothelial", x = .) %>%
  gsub(pattern = " 5", replacement = "CAF", x = .) %>%
  gsub(pattern = " 6", replacement = "NK cell", x = .) %>%
  gsub(pattern = " 0", replacement = "tumor", x = .) # assuming anything that isn't any of the other cell types = tumor

Tirosh_data <- AddMetaData(Tirosh_data,
                           metadata = celltype_metadata %>% deframe(),
                           col.name = "cell_type")

# Normalization, dimensionality reduction and clustering
Tirosh_data <- SCTransform(Tirosh_data, verbose = F)
Tirosh_data <- FindVariableFeatures(object = Tirosh_data, selection.method = "vst", verbose = F)
Tirosh_data <- ScaleData(object = Tirosh_data, features = rownames(Tirosh_data), verbose = F)
Tirosh_data <- RunPCA(object = Tirosh_data, features = VariableFeatures(Tirosh_data), verbose = F)
ElbowPlot(Tirosh_data) # choosing 15 PCs for downstream analyses
Tirosh_data <- RunUMAP(Tirosh_data, dims = 1:15, verbose = F)
Tirosh_data <- FindNeighbors(Tirosh_data, dims = 1:15, verbose = F)
Tirosh_data <- FindClusters(Tirosh_data, verbose = T)












