## Processing of Visium data and creation of integrated spatial object.

library(Seurat)
library(tidyverse)
library(scales)
library(pals)
library(cowplot)
library(fgsea)
library(msigdbr)
library(data.table)
library(mvhspatialplots)



## Processing of individual samples.
A1_new <- Load10X_Spatial(data.dir = '/Users/hunterm/Dropbox/MH_ST/visium_data_raw/MVH02_B1_BRAF_EGFP/')
B1_new <- Load10X_Spatial(data.dir = '/Users/hunterm/Dropbox/MH_ST/visium_data_raw/MVH03_C1_BRAF_EGFP/')
C1_new <- Load10X_Spatial(data.dir = '/Users/hunterm/Dropbox/MH_ST/visium_data_raw/MVH01_A1_BRAF_EGFP/')

A1_new$orig.ident <- "A1"
B1_new$orig.ident <- "B1"
C1_new$orig.ident <- "C1"

## Run normalization, dimensionality reduction and clustering.
# A1 sample:
A1_new <- SCTransform(A1_new, assay = "Spatial", verbose = F)
A1_new <- RunPCA(A1_new, assay = "SCT", verbose = F)
A1_new <- FindNeighbors(A1_new, reduction = "pca", dims = 1:6, verbose = F)
A1_new <- FindClusters(A1_new, verbose = F)
A1_new <- RunUMAP(A1_new, dims = 1:6, verbose = F)

# B1 sample:
B1_new <- SCTransform(B1_new, assay = "Spatial", verbose = F)
B1_new <- RunPCA(B1_new, assay = "SCT", verbose = F)
B1_new <- FindNeighbors(B1_new, reduction = "pca", dims = 1:6, verbose = F)
B1_new <- FindClusters(B1_new, verbose = F)
B1_new <- RunUMAP(B1_new, dims = 1:6, verbose = F)

# C1 sample:
C1_new <- SCTransform(C1_new, assay = "Spatial", verbose = F)
C1_new <- RunPCA(C1_new, assay = "SCT", verbose = F)
C1_new <- FindNeighbors(C1_new, reduction = "pca", dims = 1:6, verbose = F)
C1_new <- FindClusters(C1_new, verbose = F)
C1_new <- RunUMAP(C1_new, dims = 1:6, verbose = F)


## Integrate samples A, B and C.

# Increase memory to facilitate integration:
options(future.globals.maxSize = 8000 * 1024^2)



# Run SCTransform integration workflow: https://satijalab.org/seurat/archive/v3.0/integration.html

ABC.list <- list(A1_new, B1_new, C1_new)
genes.common <- Reduce(intersect, list(rownames(A1_new), rownames(B1_new), rownames(C1_new)))
list.features <- SelectIntegrationFeatures(object.list = ABC.list,
                                           nfeatures = 3000,
                                           assay = c("SCT", "SCT", "SCT"))
ABC.list <- PrepSCTIntegration(object.list = ABC.list,
                               anchor.features = list.features,
                               assay = "SCT",
                               verbose = F)
ABC.anchors <- FindIntegrationAnchors(object.list = ABC.list,
                                      normalization.method = "SCT",
                                      anchor.features = list.features,
                                      verbose = F)                       
ABC <- IntegrateData(anchorset = ABC.anchors,
                     features.to.integrate = genes.common,
                     normalization.method = "SCT", 
                     verbose = F)

# Rerun dimensionality reduction and clustering on integrated object.
ABC <- RunPCA(ABC, verbose = F)
ABC <- FindNeighbors(ABC, reduction = "pca", dims = 1:30)
ABC <- FindClusters(ABC, resolution = 0.8, verbose = F)
ABC <- RunUMAP(ABC, reduction = "pca", dims = 1:30)


# Add tissue type assignments to metadata.
tumor.clust <- c(1,3,6)
interface.clust <- c(11,16,17)
muscle.clust <- c(0,2)

Idents(ABC) <- "integrated_snn_res.0.8"
tumor.cells <- WhichCells(ABC, idents = tumor.clust)
tumor <- rep("tumor", length(tumor.cells))
tumor <- cbind(tumor.cells, tumor)
interface.cells <- WhichCells(ABC, idents = interface.clust)
interface <- rep("interface", length(interface.cells))
interface <- cbind(interface.cells, interface)
muscle.cells <- WhichCells(ABC, idents = muscle.clust)
muscle <- rep("muscle", length(muscle.cells))
muscle <- cbind(muscle.cells, muscle)
other.cells <- WhichCells(ABC, idents = c(1,3,6,11,16,17,0,2), invert = T)
other <- rep("other", length(other.cells))
other <- cbind(other.cells, other)

# combine all
all.metadata <- data.frame(rbind(tumor, muscle, interface, other)) %>% deframe()

# add tissue type metadata from old object to new object
ABC <- AddMetaData(object = ABC,
                       metadata = all.metadata,
                       col.name = "tissue.type")
