library(Seurat)
library(Signac)
#devtools::install_github("GreenleafLab/ArchR", ref="master", repos = BiocManager::repositories())
library(ArchR)
library(GenomicRanges)
library(EnsDb.Mmusculus.v79)
library(BSgenome.Mmusculus.UCSC.mm10)
library(DropletUtils)
library(JASPAR2020)
library(TFBSTools)

library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(gplots)
library(gridExtra)

#remotes::install_github('chris-mcginnis-ucsf/DoubletFinder')
library(DoubletFinder)
library(scDblFinder)

library(TCseq)
library(gprofiler2)
library(UpSetR)

library(Matrix)
library(patchwork)
library(parallel)

library(rvest)
library(jsonlite)


set.seed(1234)

# get gene annotations for mm10
annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
seqlevelsStyle(annotation) <- "UCSC"

# chromvar
pfm_jaspar <- getMatrixSet(
  x = JASPAR2020,
  opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE)
)

palette <- c( "#FA9FB5", "#CC79A7",
              "#8b167c", "#ff595e",
              "#ffcc00", "#f89f00", 
              "#2C558F", "#4EB3D3",
              "#76c8b1", "#0c7488"
)

# for modisco report
human_to_mouse_dict <- list(
  "Bha15" = "Bach2",
  "Zn143" = "Zfp143",
  "Tyy1" = "Tcf7l1",
  "Znf76" = "Zfp76",
  "Zbt7a" = "Zbtb7a",
  "Znf524" = "Zfp524",
  "Prd16" = "Prdm16",
  "Zn528" = "Zfp528",
  "Itf2" = "Tcf4",
  "Gcr" = "Nr3c1",
  "Po3f1" = "Pou3f1",
  "Cot2" = "Map3k8",
  "Zn317" = "Zfp317",
  "Prd14" = "Prdm14",
  "Rorg" = "Rorc",
  "Zn264" = "Zfp264",
  "Coe1" = "Ebf1",
  "Hxa9" = "Hoxa9",
  "Zn436" = "Zfp436",
  "Zn322" = "Zfp322",
  "Kaiso" = "Zbtb33",
  "Znf143" = "Zfp143",
  "Znf8" = "Zfp8",
  "Nfac2" = "Nfatc2",
  "Zn667" = "Zfp667",
  "Po3f2" = "Pou3f2",
  "Znf384" = "Zfp384",
  "Htf4" = "Tcf4",
  "Po5f1" = "Pou5f1",
  "Zn449" = "Zfp449",
  "Zn431" = "Zfp431",
  "Tha" = "Thra",
  "Ap2b" = "Tfap2b",
  "Tha11" = "Thra1",
  "Andr" = "Ar",
  "Znf423" = "Zfp423",
  "Zn320" = "Zfp320",
  "Ap2c" = "Tfap2c",
  "Zn350" = "Zfp350",
  "Zn341" = "Zfp341",
  "Zbt18" = "Zbtb18",
  "Hxc9" = "Hoxc9"
)

replace_human_genes_with_mouse <- function(df, human_to_mouse_dict) {
  for (human_gene in names(human_to_mouse_dict)) {
    df[] <- lapply(df, function(x) ifelse(x == human_gene, human_to_mouse_dict[[human_gene]], x))
  }
  return(df)
}