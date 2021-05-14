# GSEA of interface clusters

library(Seurat)
library(tidyverse)
library(scales)
library(pals)
library(cowplot)
library(fgsea)
library(msigdbr)
library(data.table)
library(mvhspatialplots)

# Load pathways for GSEA
human.genes <- msigdbr(species = "Homo sapiens")
genesets.interest <- filter(human.genes, gs_cat == "H" | gs_subcat == "CP:KEGG" | gs_subcat == "CP:REACTOME" | gs_subcat == "BP")
pathways.interest <- genesets.interest %>% split(x = .$gene_symbol, f = .$gs_name)
# GO biological pathways
genesets.GOBP <- filter(human.genes, gs_subcat == "GO:BP")
pathways.GOBP <- genesets.GOBP %>% split(x = .$gene_symbol, f = .$gs_name)
# GO CC
genesets.GOCC <- filter(human.genes, gs_subcat == "GO:CC")
pathways.GOCC <- genesets.GOCC %>% split(x = .$gene_symbol, f = .$gs_name)


## GSEA of scRNA-seq interface cluster
# Find markers for interface cluster
Idents(EF.filt) <- "integrated_cell_type"
EF_interface_markers <- FindMarkers(EF.filt,
                                    ident.1 = "interface") %>% 
  rownames_to_column(var = "gene") %>% 
  arrange(-avg_logFC)

# Find markers for interface (muscle-like) cluster
Idents(EF.filt) <- "integrated_cell_subtype"
EF_muscle_interface_markers <- FindMarkers(EF.filt,
                                    ident.1 = "interface (muscle-like)") %>% 
  rownames_to_column(var = "gene") %>% 
  arrange(-avg_logFC)

# Find markers for interface (tumor-like) cluster
Idents(EF.filt) <- "integrated_cell_subtype"
EF_tumor_interface_markers <- FindMarkers(EF.filt,
                                    ident.1 = "interface (tumor-like)") %>% 
  rownames_to_column(var = "gene") %>% 
  arrange(-avg_logFC)

# Run GSEA with GO cellular component pathway set
EF_interface_GOCC <- fgsea(pathways = pathways.GOCC,
                           stats = EF_interface_markers %>% convert_FindMarkers_human(.),
                           nperm = 10000) %>%
  arrange(desc(NES)) %>%
  mutate(log10pval = -log10(pval)) %>%
  select(pathway, NES, pval, log10pval, everything())

EF_muscle_interface_GOCC <- fgsea(pathways = pathways.GOCC,
                           stats = EF_muscle_interface_markers %>% convert_FindMarkers_human(.),
                           nperm = 10000) %>%
  arrange(desc(NES)) %>%
  mutate(log10pval = -log10(pval)) %>%
  select(pathway, NES, pval, log10pval, everything())

EF_tumor_interface_GOCC <- fgsea(pathways = pathways.GOCC,
                           stats = EF_tumor_interface_markers %>% convert_FindMarkers_human(.),
                           nperm = 10000) %>%
  arrange(desc(NES)) %>%
  mutate(log10pval = -log10(pval)) %>%
  select(pathway, NES, pval, log10pval, everything())

# Run GSEA with GO biological processes pathway set
EF_interface_GOBP <- fgsea(pathways = pathways.GOBP,
                           stats = EF_interface_markers %>% convert_FindMarkers_human(.),
                           nperm = 10000) %>%
  arrange(desc(NES)) %>%
  mutate(log10pval = -log10(pval)) %>%
  select(pathway, NES, pval, log10pval, everything())

EF_muscle_interface_GOBP <- fgsea(pathways = pathways.GOBP,
                                  stats = EF_muscle_interface_markers %>% convert_FindMarkers_human(.),
                                  nperm = 10000) %>%
  arrange(desc(NES)) %>%
  mutate(log10pval = -log10(pval)) %>%
  select(pathway, NES, pval, log10pval, everything())

EF_tumor_interface_GOBP <- fgsea(pathways = pathways.GOBP,
                                 stats = EF_tumor_interface_markers %>% convert_FindMarkers_human(.),
                                 nperm = 10000) %>%
  arrange(desc(NES)) %>%
  mutate(log10pval = -log10(pval)) %>%
  select(pathway, NES, pval, log10pval, everything())


## GSEA of SRT interface clusters
# Find markers for interface cluster
Idents(ABC) <- "tissue.type"
ABC_interface_markers <- FindMarkers(ABC,
                                    ident.1 = "interface") %>% 
  rownames_to_column(var = "gene") %>% 
  arrange(-avg_logFC)

# Run GSEA with GO cellular component pathway set
ABC_interface_GOCC <- fgsea(pathways = pathways.GOCC,
                           stats = ABC_interface_markers %>% convert_FindMarkers_human(.),
                           nperm = 10000) %>%
  arrange(desc(NES)) %>%
  mutate(log10pval = -log10(pval)) %>%
  select(pathway, NES, pval, log10pval, everything())

# Run GSEA with GO biological processes pathway set
ABC_interface_GOBP <- fgsea(pathways = pathways.GOBP,
                           stats = ABC_interface_markers %>% convert_FindMarkers_human(.),
                           nperm = 10000) %>%
  arrange(desc(NES)) %>%
  mutate(log10pval = -log10(pval)) %>%
  select(pathway, NES, pval, log10pval, everything())

