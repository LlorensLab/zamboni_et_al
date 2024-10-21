# filtering
source("../source.R")

object <- readRDS("unmerged_obj_multiome.rds")
dirs <- list.dirs("data/", recursive = F, full.names = F)

# remove low quality nuclei
# remove nuclei that are above or below the 3IQR for each cluster
is_outlier <- function(feature, tail = "both"){
  # tail = c("both", "upper", "lower")
  upper = (feature > quantile(feature, .75) + 3*IQR(feature, na.rm = T))
  lower = (feature < quantile(feature, .25) - 3*IQR(feature, na.rm = T))
  
  if(tail == "upper") outliers = upper
  if(tail == "lower") outliers = lower
  if(tail == "both")  outliers = upper | lower
  return(outliers)
}

for(i in seq_along(object)){
  # remove cells with missing values in meta
  na_cells <- object[[i]]@meta.data[!complete.cases(object[[i]]@meta.data),] %>% rownames()
  object[[i]] <- subset(object[[i]], cells = na_cells, invert = TRUE)
  
  # identify groupwise outliers for each metric
  object[[i]]@meta.data <- object[[i]]@meta.data %>%
    group_by(seurat_clusters) %>%
    mutate(is_outlier_frag = is_outlier(nCount_peaks, tail = "both"),
           is_outlier_umi = is_outlier(nCount_RNA_raw, tail = "both"),
           is_outlier_genes = is_outlier(nFeature_RNA_raw, tail = "both"),
           is_outlier_mt = is_outlier(percent.mt, tail = "upper"),
           is_outlier_tss = is_outlier(TSS.enrichment, tail = "lower"),
           is_outlier_nucleosome = is_outlier(nucleosome_signal, tail = "upper")
    ) %>%
    as.data.frame() 
  
  rownames(object[[i]]@meta.data) <- colnames(object[[i]])
  
  outliers <- grep("is_outlier", colnames(object[[i]]@meta.data))
  object[[i]]$is_outlier = rowSums(object[[i]]@meta.data[outliers] == TRUE) > 0
  
}

object <- lapply(X = object, FUN = function(x) {
  x <- subset(x, subset = 
                is_outlier == "FALSE"
  )
}) # removes about 10% of cells

#subset low quality clusters manually
object$ALL_U_2 <- subset(object$ALL_U_2, idents = c("11"), invert = T)
object$DEPL_3dpi_6 <- subset(object$DEPL_3dpi_6, idents = "8", invert = T)
object$FT_28dpi_4 <- subset(object$FT_28dpi_4, idents = "5", invert = T)


object <- lapply(X = object, FUN = function(x) {
  #sig_pc <- length(which(x@reductions[["pca"]]@stdev > 2))
  x <- x %>%
    FindMultiModalNeighbors(
      reduction.list = list("pca", "lsi"), 
      dims.list = list(1:15, 2:15),
      modality.weight.name = "RNA.weight") %>%
    RunUMAP(
      nn.name = "weighted.nn",
      reduction.name = "wnn.umap", 
      reduction.key = "wnnUMAP_") %>%
    FindClusters(
      graph.name = "wsnn", 
      algorithm = 3, 
      resolution = 0.3)
  
  x <- RunUMAP(x, reduction = 'pca', dims = 1:15, assay = 'RNA', 
               reduction.name = 'rna.umap', reduction.key = 'rnaUMAP_')
  x <- RunUMAP(x, reduction = 'lsi', dims = 2:15, assay = 'macs2_peaks', 
               reduction.name = 'lsi.umap', reduction.key = 'lsiUMAP_')
  
  x <- FindClusters(x,
                    graph.name = "wsnn", 
                    algorithm = 3, 
                    resolution = 0.2)
  
})

saveRDS(object, "data/unmerged_multiome_cb_clean.rds")


#doublet
# RNA doublets
#best pK values obtained from paramsweep
best_pK <- c(0.005, 0.17, 0.18, 0.04, 0.01, 0.005, 
                    0.005, 0.01, 0.02, 0.005, 0.005, 0.02, 
                    0.03, 0.005, 0.01, 0.01, 0.005, 0.01, 0.14,
                    0.07, 0.08)

for(i in seq_along(object)){
  x <- object[[i]]
  best_pK <- best_pK[i]
  
  homotypic.prop <- modelHomotypic(x$seurat_clusters)           ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
  doublet_rate <- (0.0008*nrow(x@meta.data)-0.032)/100 # based on 10x linear estimate
  nExp_poi <- round(doublet_rate*nrow(x@meta.data))  
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
  
  sig_pc <- 15
  x <- doubletFinder_v3(x, 
                        PCs = 1:sig_pc, 
                        pN = 0.25, 
                        pK = as.numeric(as.character(best_pK)), 
                        nExp = nExp_poi, 
                        reuse.pANN = FALSE, 
                        sct = FALSE)
  x <- doubletFinder_v3(x, 
                        PCs = 1:sig_pc, 
                        pN = 0.25, 
                        pK = as.numeric(as.character(best_pK)), 
                        nExp = nExp_poi.adj, 
                        reuse.pANN = colnames(x@meta.data)[grep(colnames(x@meta.data), 
                                                                pattern = "pANN_")], 
                        sct = FALSE)
  
  df_class <- grep(colnames(x@meta.data), 
                   pattern = "DF.classifications_")
  
  x$DF_doublet <- ifelse(x[[colnames(x@meta.data)[df_class[1]]]] == "Doublet" &
                           x[[colnames(x@meta.data)[df_class[2]]]] == "Doublet",
                         "Doublet", "Singlet")

  object[[i]] <- x
}

# ATAC doublets
for(i in seq_along(object)){
  x <- object[[i]]
  doublet_rate <- (0.0008*nrow(x@meta.data)-0.032)/100 # based on 10x linear estimate
  artificialDoublets <- round(doublet_rate*nrow(x@meta.data)) 
  sig_pc <- 15
  scdbl_obj <- scDblFinder(x@assays$cellranger_peaks@counts,
                           clusters = x$seurat_clusters,
                           nfeatures = 1000, #aggregated features
                           dims = 1:sig_pc,
                           artificialDoublets = artificialDoublets,
                           aggregateFeatures = T,
                           verbose = T,
                           processing="normFeatures")
  
  scdbl_out <- as.data.frame(scdbl_obj@colData@listData)
  rownames(scdbl_out) <- scdbl_obj@colData@rownames
  x <- AddMetaData(x, scdbl_out)
  
  object[[i]] <- x
}

# remove unrelated cells (from another experiment)
object[[18]] <- subset(object[[18]], idents = c("0", "1", "3"))
object[[19]] <- subset(object[[19]], idents = c("1"))
object[[20]] <- subset(object[[20]], idents = c("0", "1"))
object[[21]] <- subset(object[[21]], idents = c("0"))

saveRDS(object, "unmerged_obj_multiome_clean.rds")

