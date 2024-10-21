library(Seurat)
library(Signac)
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

library(Matrix)
library(patchwork)
library(parallel)

set.seed(1234)

# get gene annotations for mm10
annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
seqlevelsStyle(annotation) <- "UCSC"

