## Import and filtering of scRNA-seq reactions

library(Seurat)
library(tidyverse)
library(pals)
library(scales)


## Import counts matrix and make Seurat object for both reactions. 
# Note: calling these samples E1 and F1 since the ST samples are ABCD

# E1: (this is pretty stringent filtering but the same parameters work well for F1)
E1_counts <- Read10X(data.dir = "/Users/hunterm/Dropbox/10x/Rxn1-Fish2/raw_feature_bc_matrix/")
E1 <- CreateSeuratObject(counts = E1_counts,
                         min.features = 75,
                         min.cells = 1,
                         project = "E1")

E1[["percent.mt"]] <- PercentageFeatureSet(E1, pattern = "^mt-")
ribo_genes <- grep(pattern = "^rps|^rpl", x = rownames(E1), value = T)
E1[["percent.ribosomal"]] <- PercentageFeatureSet(E1, features = ribo_genes)
E1 <- subset(E1, subset = percent.mt <= 20)
E1 <- SCTransform(object = E1, verbose = F)
E1 <- FindVariableFeatures(object = E1, selection.method = "vst")
E1 <- ScaleData(object = E1, features = rownames(E1))
E1 <- RunPCA(object = E1, features = VariableFeatures(E1))
ElbowPlot(E1)
E1 <- FindNeighbors(E1, dims = 1:15, verbose = F)
E1 <- FindClusters(E1, resolution = 0.2)
E1 <- RunUMAP(E1, dims = 1:15, verbose = F)

# F1
F1_counts <- Read10X(data.dir = "/Users/hunterm/Dropbox/10x/Rxn2-Fish3/raw_feature_bc_matrix/")
F1 <- CreateSeuratObject(counts = F1_counts,
                         min.features = 75,
                         min.cells = 1,
                         project = "F1")
F1[["percent.mt"]] <- PercentageFeatureSet(F1, pattern = "^mt-")
ribo.genes <- grep(pattern = "^rps|^rpl|^mrpl|^mrps", x = rownames(F1), value = T)
F1[["percent.ribosomal"]] <- PercentageFeatureSet(F1, pattern = "^rps|^rpl|^mrpl|^mrps")
F1 <- subset(F1, subset = percent.mt <= 20)
F1 <- SCTransform(object = F1, verbose = F)
F1 <- FindVariableFeatures(object = F1, selection.method = "vst")
F1 <- ScaleData(object = F1, features = rownames(F1))
F1 <- RunPCA(object = F1, features = VariableFeatures(F1))
ElbowPlot(F1)
F1 <- FindNeighbors(F1, dims = 1:15, verbose = F)
F1 <- FindClusters(F1, resolution = 0.7)
F1 <- RunUMAP(F1, dims = 1:15, verbose = F)


# Remove cluster of tumor cells with very low number of UMIs.
# cluster 1 both E1 and F1
Idents(E1) <- "seurat_clusters"
Idents(F1) <- "seurat_clusters"

E1.filt <- subset(E1, idents = 1, invert = T)
F1.filt <- subset(F1, idents = 1, invert = T)

# Rerun dimensionality reductions
E1.filt <- FindVariableFeatures(E1.filt, selection.method = "vst", nFeatures = 2000)
E1.filt <- ScaleData(E1.filt, features = rownames(E1.filt))
E1.filt <- RunPCA(E1.filt, features = VariableFeatures(E1.filt))

F1.filt <- FindVariableFeatures(F1.filt, selection.method = "vst", nFeatures = 2000)
F1.filt <- ScaleData(F1.filt, features = rownames(F1.filt))
F1.filt <- RunPCA(F1.filt, features = VariableFeatures(F1.filt))


## Integrate the two datasets.
EF_list <- list(E1.filt, F1.filt)
EF_genes <- intersect(rownames(E1.filt), rownames(F1.filt))
EF_features <- SelectIntegrationFeatures(object.list = EF_list,
                                         nfeatures = 3000)
EF_list <- PrepSCTIntegration(object.list = EF_list,
                              anchor.features = EF_features,
                              verbose = F)
EF_anchors <- FindIntegrationAnchors(object.list = EF_list,
                                     normalization.method = "SCT",
                                     anchor.features = EF_features,
                                     verbose = F)
EF.filt <- IntegrateData(anchorset = EF_anchors,
                         features.to.integrate = EF_genes,
                         normalization.method = "SCT", 
                         verbose = FALSE)


EF.filt <- ScaleData(EF.filt, features = rownames(EF.filt), verbose = F)
EF.filt <- RunPCA(EF.filt, features = VariableFeatures(EF.filt))
EF.filt <- FindClusters(EF.filt, resolution = 0.25)
EF.filt <- RunUMAP(EF.filt, reduction = "pca", dims = 1:15, verbose = F)

# Add cluster annotations to metadata.
# cell_type
Idents(EF.filt) <- "integrated_snn_res.0.25"
tumor <- WhichCells(EF.filt, idents = c(0,2)) %>% data.frame() %>% mutate(integrated_cell_type = "tumor")
interface <- WhichCells(EF.filt, idents = c(3,11)) %>% data.frame() %>% mutate(integrated_cell_type = "interface")
ker <- WhichCells(EF.filt, idents = c(6,10)) %>% data.frame() %>% mutate(integrated_cell_type = "keratinocytes")
ery <- WhichCells(EF.filt, idents = c(1,5)) %>% data.frame() %>% mutate(integrated_cell_type = "erythrocytes")
unknown <- WhichCells(EF.filt, idents = 9) %>% data.frame() %>% mutate(integrated_cell_type = "unknown")
macro <- WhichCells(EF.filt, idents = 4) %>% data.frame() %>% mutate(integrated_cell_type = "macrophages")
Tcells <- WhichCells(EF.filt, idents = 8) %>% data.frame() %>% mutate(integrated_cell_type = "T cells")
neut <- WhichCells(EF.filt, idents = 7) %>% data.frame() %>% mutate(integrated_cell_type = "neutrophils")
celltypes <- rbind(tumor, interface, ker, ery, unknown, macro, Tcells, neut) %>% deframe()
rm(tumor, interface, ker, ery, unknown, macro, Tcells, neut)

EF.filt <- AddMetaData(EF.filt,
                       metadata = celltypes,
                       col.name = "integrated_cell_type")

# cell_subtype
Idents(EF.filt) <- "integrated_snn_res.0.25"
tumor <- WhichCells(EF.filt, idents = c(0,2) %>% data.frame() %>% mutate(integrated_cell_subtype = "tumor")
                    cilia <- WhichCells(EF.filt, idents = 11) %>% data.frame() %>% mutate(integrated_cell_subtype = "interface (muscle)")
                    div <- WhichCells(EF.filt, idents = 3) %>% data.frame() %>% mutate(integrated_cell_subtype = "interface (tumor)")
                    ker <- WhichCells(EF.filt, idents = c(6,10)) %>% data.frame() %>% mutate(integrated_cell_subtype = "keratinocytes")
                    ery <- WhichCells(EF.filt, idents = c(1,5)) %>% data.frame() %>% mutate(integrated_cell_subtype = "erythrocytes")
                    unknown <- WhichCells(EF.filt, idents = 9) %>% data.frame() %>% mutate(integrated_cell_subtype = "unknown")
                    macro <- WhichCells(EF.filt, idents = 4) %>% data.frame() %>% mutate(integrated_cell_subtype = "macrophages")
                    Tcells <- WhichCells(EF.filt, idents = 8) %>% data.frame() %>% mutate(integrated_cell_subtype = "T cells")
                    neut <- WhichCells(EF.filt, idents = 7) %>% data.frame() %>% mutate(integrated_cell_subtype = "neutrophils")
                    cellsubtypes <- rbind(tumor, cilia, div, ker, ery, unknown, macro, Tcells, neut) %>% deframe()
                    rm(tumor, cilia, div, ker, ery, unknown, macro, Tcells, neut)
                    
                    EF.filt <- AddMetaData(EF.filt,
                                           metadata = cellsubtypes,
                                           col.name = "integrated_cell_subtype")
                    
                    