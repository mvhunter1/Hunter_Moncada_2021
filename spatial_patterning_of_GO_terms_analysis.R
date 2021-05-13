# Spatial patterning of GO terms analysis

library(Seurat)
library(tidyverse)
library(msigdbr)
library(rdist)


# Load integrated Visium object
load("/Users/hunterm/Dropbox/MH_ST/Miranda_R/visium_Seurat_objs_fish/ABC_integrated_allgenes_new.R")
ST.int <- ABC

# Load metadata of the refined cluster annotations
ABC.meta <- read.csv('/Users/hunterm/Dropbox/MH_ST/Miranda_R/visium_Seurat_objs_fish/ABC_metadata_updated.csv')
ABC.meta$barcodes <- sapply(ABC.meta$barcodes, sub, pattern = '.1', replacement = '-1') 
rownames(ABC.meta) <- ABC.meta$barcodes

Idents(ST.int) <- ABC.meta$tissue.type # Rename tissue types

# Retrieve GO annotations for D. rerio
zf.genes <- msigdbr(species = "Danio rerio")
pathways.GOBP <- filter(zf.genes, gs_subcat == "BP") %>% split(x = .$gene_symbol, f = .$gs_name)


# Set variables to determine the sample and tissue region to perform analyses on
sample <- 'B1' # Change this to any of the following to switch between datasets: A1, B1, C1
spot.idx <- ST.int$orig.ident == sample

regions <- unique(ABC.meta[spot.idx,'tissue.type'])
region.of.interest <- c('interface', 'muscle', 'other')  # Change to any of 'regions' entries to switch

Idents(ST.int) <- "orig.ident"
SCT <- subset(ST.int, idents = sample)
roi.idx <- Idents(SCT)==region.of.interest
SCT$new.tissue.types <- Idents(SCT)

Idents(SCT) <- "tissue.type"
SCT <- subset(SCT, idents = region.of.interest)
exp <- as.data.frame(SCT@assays$SCT@data)


# Initialize data.frame to store results of analysis
GO.pval <- data.frame(matrix(1,nrow = length(names(pathways.GOBP)), ncol = 1),
                      row.names = names(pathways.GOBP))
colnames(GO.pval) <- 'p.value'

GO.expression <- data.frame(matrix(0,nrow = length(names(pathways.GOBP)), ncol = dim(exp)[2]),
                            row.names = names(pathways.GOBP))
colnames(GO.expression) <- names(exp)


##  Run analysis:
# Take a given GO term, compute average expression
for (GO in 1:length(names(pathways.GOBP))) {
  
  GO.genes <- pathways.GOBP[[GO]]
  GO.genes <- GO.genes[GO.genes %in% row.names(exp)] # Filter GO genes that aren't in the expression matrix. 
  
  GO.avg.exp <- apply(exp[GO.genes,], 2, mean)
  
  # Find spots that highly express this GO term (mean + 1.5 * std in this case)
  avg <- mean(GO.avg.exp)
  std <- sd(GO.avg.exp)
  exp.thresh <- avg + (1.5*std)
  
  GO.high.spots <- which(as.data.frame(GO.avg.exp) >= exp.thresh)
  GO.high.spots.bc <- colnames(SCT)[GO.high.spots]
  
  # Are there at least X number of spots?
  num.spots <- 5
  if (length(GO.high.spots.bc) < num.spots){
    next
  }
  
  # Compute the distance between these spots' XY coordinates
  coords <- rbind(ABC.meta[GO.high.spots.bc,'X'],
                  ABC.meta[GO.high.spots.bc,'Y'])
  distances <- rdist(t(coords), metric = 'euclidean') # Distance computation
  
  # Compute distances between equal number of random spots
  set.seed(1) # set seed for random number generator
  
  spot.idx <- ST.int$orig.ident == sample
  
  # Get list of barcodes for tumor region of sample
  bc<- intersect(names(spot.idx)[spot.idx],
                 ABC.meta$barcodes[ABC.meta$tissue.type == region.of.interest])
  all.coord <- cbind(as.data.frame(ABC.meta[bc,'X']),
                     as.data.frame(ABC.meta[bc,'Y']))
  
  if (dim(all.coord)[1] < length(GO.high.spots.bc)) {
    next
  }
  
  rand.idx <- sample(x = 1:dim(all.coord)[1],
                     size = length(GO.high.spots.bc))
  rand.coord <- all.coord[rand.idx,]
  
  rand.dist <- rdist(rand.coord, metric = 'euclidean') # Distance computation for random spots
  
  # Compare distance
  pval <- wilcox.test(distances, rand.dist, alternative = 'greater')
  
  GO.pval[names(pathways.GOBP[GO]), 'p.value']  <- pval$p.value
  GO.expression[names(pathways.GOBP[GO]),] <- apply(exp[GO.genes,], 2, mean)
  
  print(paste("GO term", GO, "out of", length(names(pathways.GOBP))))   
}

GO.pval <- GO.pval %>% arrange(p.value)