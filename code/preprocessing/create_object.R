# make seurat objects
source("../source.R")

#load feature matrices
mtx_cluster <- readRDS("featurematrix_archr_clusterpeaks_multiome.rds")

path <- "data/"
dirs <- list.dirs(path, recursive = F, full.names = F)
inputFiles <- paste0(path, dirs, "/outs/atac_fragments.tsv.gz")
names(inputFiles) <- dirs

#create peaks object
#filter nuclei based on low count
unmerged_obj <- list()
for(i in seq_along(dirs)){
  assay <- CreateChromatinAssay(
    counts = mtx_cluster[[i]],
    fragments = inputFiles[i], 
    annotation = annotation
  )
  
  x <- CreateSeuratObject(assay, assay = "peaks")
  x$sample <- dirs[i]
  
  print(paste0("All cells: ", ncol(x)))
  x <- subset(x, subset = nCount_peaks > 300)
  print(paste0("Number of cells after filter (nCount_peaks > 300: ", ncol(x)))
  unmerged_obj[[i]] <- x
}

# add GEX data (with and without ambient RNA removal)
for(i in dirs){
  multiome <- unmerged_obj[[i]]
  
  #load gene expression matrix
  counts_gex <- Read10X(paste0(path, i, "/outs/filtered_feature_bc_matrix")) %>%
    .$`Gene Expression`
  #colnames(counts_gex) <- paste0(i, "_", colnames(counts_gex))
  # change tdTom naming (- will throw error because it'll try to separate into columns)
  rownames(counts_gex) <- gsub(pattern = "tdTomato-WPRE", 
                               replacement = "tdTomato.WPRE", 
                               x = rownames(counts_gex),
                               ignore.case=T)
  counts_gex <- counts_gex[,colnames(multiome)]
  
  multiome[["RNA_raw"]] <- CreateAssayObject(counts_gex)
  multiome <- subset(multiome, subset = nCount_RNA_raw > 200)
  
  # load cellbended gene expression matrix
  counts_cb <- Read10X_h5(paste0(path, i, "/outs/cellbended_gex_matrix_fpr001_seurat.h5"))
  #colnames(counts_cb) <- paste0(i, "_", colnames(counts_cb))
  # change tdTom naming (- will throw error because it'll try to separate into columns)
  rownames(counts_cb) <- gsub(pattern = "tdTomato-WPRE", 
                              replacement = "tdTomato.WPRE", 
                              x = rownames(counts_cb),
                              ignore.case=T)
  counts_cb <- counts_cb[,colnames(multiome)]
  
  multiome[["RNA"]] <- CreateAssayObject(counts_cb)
  unmerged_obj[[i]] <- multiome
  
}

# add atac qc metrics and percent mt
unmerged_obj <- lapply(X = unmerged_obj, FUN = function(x) {
  DefaultAssay(x) <- "peaks"
  x <- NucleosomeSignal(x)
  x <- TSSEnrichment(x)
  
  DefaultAssay(x) <- "RNA"
  x[["percent.mt"]] <- PercentageFeatureSet(x, pattern = "^mt-")
  
})


unmerged_obj <- lapply(X = unmerged_obj, FUN = function(x) {
  # process RNA and peaks
  DefaultAssay(x) <- "RNA"
  x <- x %>% 
    NormalizeData() %>% 
    FindVariableFeatures() %>% 
    ScaleData() %>% 
    RunPCA()
  
  DefaultAssay(x) <- "peaks"
  x <- x %>%
    FindTopFeatures(min.cutoff = 5) %>%
    RunTFIDF() %>%
    RunSVD()
  
})

# build a joint neighbor graph and UMAP visualization using both assays
unmerged_obj <- lapply(X = unmerged_obj, FUN = function(x) {
  sig_pc_pca <- length(which(x@reductions[["pca"]]@stdev > 2))
  sig_pc_lsi <- length(which(x@reductions[["lsi"]]@stdev > 2))
  x <- x %>%
    FindMultiModalNeighbors(
      reduction.list = list("pca", "lsi"), 
      dims.list = list(1:sig_pc_pca, 2:sig_pc_lsi),
      modality.weight.name = "RNA.weight") %>%
    RunUMAP(
      nn.name = "weighted.nn",
      reduction.name = "wnn.umap", 
      reduction.key = "wnnUMAP_") %>%
    FindClusters(
      graph.name = "wsnn", 
      algorithm = 3, 
      resolution = 0.3)
  
  x <- RunUMAP(x, reduction = 'pca', dims = 1:sig_pc_pca, assay = 'RNA', 
               reduction.name = 'rna.umap', reduction.key = 'rnaUMAP_')
  x <- RunUMAP(x, reduction = 'lsi', dims = 2:sig_pc_lsi, assay = 'peaks', 
               reduction.name = 'lsi.umap', reduction.key = 'lsiUMAP_')
  
  x <- FindClusters(x,
                    graph.name = "wsnn", 
                    algorithm = 3, 
                    resolution = 0.2)
  
})

# save unmerged objects
saveRDS(unmerged_obj, "unmerged_obj_multiome.rds")

