# Import, normalization and clustering of snRNA-seq data

library(Seurat)
library(tidyverse)
library(scales)
library(pals)
library(cowplot)
library(fgsea)
library(msigdbr)
library(data.table)
library(mvhspatialplots)

# Create Seurat object
NS_counts <- Read10X_h5("/Users/hunterm/Dropbox/10x/nuc_seq/cellranger/outs/filtered_feature_bc_matrix.h5")
NS_new <- CreateSeuratObject(counts = NS_counts$`Gene Expression`,
                             project = "nuc-seq",
                             min.cells = 3,
                             min.features = 200)

# Add hashing data
NS_new[['HTO']] <- CreateAssayObject(counts = NS_counts$`Antibody Capture`[, colnames(NS_new)]) 

# Some QC metrics
NS_new[["percent.mt"]] <- PercentageFeatureSet(NS_new, pattern = "^mt-")
ribo_genes <- grep(pattern = "^rps|^rpl", x = rownames(NS_new), value = T)
NS_new[["percent.ribosomal"]] <- PercentageFeatureSet(NS_new, features = ribo_genes)

# Analysis of HTO data
NS_new <- SCTransform(NS_new, assay = "HTO", verbose = F)
NS_new <- HTODemux(NS_new, assay = "HTO", positive.quantile = 0.99)

# Run doubletFinder to check for potential doublets
library(DoubletFinder)
sweep.res.list_NS <- paramSweep_v3(NS_new, PCs = 1:15, sct = T)
sweep.stats_NS <- summarizeSweep(sweep.res.list_NS, GT = F)
bcmvn_NS <- find.pK(sweep.stats_NS)
homotypic.prop <- modelHomotypic(NS_new@meta.data$cell_type)
nExp_poi <- round(0.075*nrow(NS_new@meta.data))
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
NS_new <- doubletFinder_v3(NS_new, PCs = 1:15, pN = 0.25, pK = 0.09, nExp = nExp_poi, reuse.pANN = F, sct = T)
NS_new <- doubletFinder_v3(NS_new, PCs = 1:15, pN = 0.25, pK = 0.09, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.09_862", sct = T)

# There are a lot of doublets in the snRNA-seq dataset, so will filter them out before continuing with downstream analyses.
Idents(NS_new) <- "DF.classifications_0.25_0.09_747"
NS_new <- subset(NS_new, idents = "Singlet")

# Normalization and clustering of RNA data
# Need to run this after HTO demultimplexing, otherwise it will overwrite the SCT slot.
NS_new <- SCTransform(object = NS_new, verbose = F, assay = "RNA")
NS_new <- FindVariableFeatures(object = NS_new, selection.method = "vst")
NS_new <- ScaleData(object = NS_new, features = rownames(NS_new))
NS_new <- RunPCA(object = NS_new, features = VariableFeatures(NS_new))
ElbowPlot(NS_new)
NS_new <- FindNeighbors(NS_new, dims = 1:15, verbose = F)
NS_new <- FindClusters(NS_new, resolution = 0.8)
NS_new <- RunUMAP(NS_new, dims = 1:15, verbose = F)
NS_new <- FindClusters(NS_new, resolution = 0.2, verbose = F)
NS_new <- FindClusters(NS_new, resolution = 0.25, verbose = F)
NS_new <- FindClusters(NS_new, resolution = 0.3, verbose = F)

# Add cluster assignments
NS_metadata <- NS_new[[]] %>% rownames_to_column(var = "cell") %>% dplyr::select(cell, SCT_snn_res.0.2)
for (ii in 1:nrow(NS_metadata)) {
  cluster <- NS_metadata[ii,2] %>% as.character() %>% as.numeric()
  if (cluster %in% c(0,10)) {
    NS_metadata$cell_type[ii] <- "liver"
  } else if (cluster == 5) {
    NS_metadata$cell_type[ii] <- "macrophages"
  } else if (cluster %in% c(1,6)) {
    NS_metadata$cell_type[ii] <- "tumor"
  } else if (cluster == 4) {
    NS_metadata$cell_type[ii] <- "interface"
  } else if (cluster == 3) {
    NS_metadata$cell_type[ii] <- "intestinal"
  } else if (cluster == 2) {
    NS_metadata$cell_type[ii] <- "keratinocytes"
  } else if (cluster == 7) {
    NS_metadata$cell_type[ii] <- "fibroblasts"
  } else if (cluster == 9) {
    NS_metadata$cell_type[ii] <- "muscle"
  } else if (cluster == 11) {
    NS_metadata$cell_type[ii] <- "endothelial"
  } else if (cluster == 13) {
    NS_metadata$cell_type[ii] <- "unknown"
  } else if (cluster == 12) {
    NS_metadata$cell_type[ii] <- "erythrocytes"
  } else if (cluster == 14) {
    NS_metadata$cell_type[ii] <- "nervous system"
  } else if (cluster == 8) {
    NS_metadata$cell_type[ii] <- "NK cells"
  }
}

NS_new <- AddMetaData(NS_new,
                      metadata = NS_metadata %>% dplyr::select(cell, cell_type) %>% deframe(),
                      col.name = "cell_type")

## Subset interface cluster only
Idents(NS_new) <- "cell_type"
NS_new_interface <- subset(NS_new, idents = "interface")
NS_new_interface <- FindVariableFeatures(NS_new_interface, selection.method = "vst", verbose = F)
NS_new_interface <- ScaleData(object = NS_new_interface, features = rownames(NS_new_interface), verbose = F)
NS_new_interface <- RunPCA(NS_new_interface, features = VariableFeatures(NS_new_interface))
ElbowPlot(NS_new_interface)
NS_new_interface <- FindNeighbors(NS_new_interface, dims = 1:10, verbose = F)
NS_new_interface <- FindClusters(NS_new_interface, resolution = 0.8, verbose = F)
NS_new_interface <- RunUMAP(NS_new_interface, dims = 1:10, verbose = F)

# Subclustering of interface nuclei
NS_new_interface <- FindClusters(NS_new_interface, resolution = 0.3, verbose = F)
NS_interface_metadata <- NS_new_interface[[]] %>% rownames_to_column(var = "cell") %>% dplyr::select(cell, SCT_snn_res.0.3)
for (ii in 1:nrow(NS_interface_metadata)) {
  cluster <- NS_interface_metadata[ii,2] %>% as.character() %>% as.numeric()
  if (cluster == 0) {
    NS_interface_metadata$cell_type[ii] <- "interface (tumor-like)"
  } else if (cluster == 1) {
    NS_interface_metadata$cell_type[ii] <- "interface (muscle-like)"
  } else if (cluster %in% c(2,4)) {
    NS_interface_metadata$cell_type[ii] <- "interface (intestinal-like)"
  } else if (cluster == 5) {
    NS_interface_metadata$cell_type[ii] <- "interface (liver-like)"
  } else if (cluster == 3) {
    NS_interface_metadata$cell_type[ii] <- "interface (immune-like)"
  } 
}
NS_new_interface <- AddMetaData(NS_new_interface,
                                metadata = NS_interface_metadata %>% dplyr::select(cell, cell_type) %>% deframe(),
                                col.name = "cell_subtype")


