library(Seurat)
library(Signac)
#devtools::install_github("GreenleafLab/ArchR", ref="master", repos = BiocManager::repositories())
library(ArchR)
library(GenomicRanges)
library(EnsDb.Mmusculus.v79)
library(BSgenome.Mmusculus.UCSC.mm10)
library(DropletUtils)

library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(gplots)
library(gridExtra)

#remotes::install_github('chris-mcginnis-ucsf/DoubletFinder')
library(DoubletFinder)
library(scDblFinder)

library(TCseq)

library(Matrix)
library(patchwork)
library(parallel)

set.seed(1234)

# get gene annotations for mm10
annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
seqlevelsStyle(annotation) <- "UCSC"

palette <- c( "#FA9FB5", "#CC79A7",
              "#8b167c", "#ff595e",
              "#ffcc00", "#f89f00", 
              "#2C558F", "#4EB3D3",
              "#76c8b1", "#0c7488"
)

