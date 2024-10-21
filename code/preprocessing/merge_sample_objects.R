source("../source.R")

object <- readRDS("unmerged_obj_multiome_clean.rds")

# merge objects
multiome <- merge(
  x = object[[1]],
  y = c(object[2:21])
)

rm(object)

multiome$experiment <- strsplit(multiome$sample, "_") %>% sapply("[[", 1)
multiome$condition <- strsplit(multiome$sample, "_") %>% sapply("[[", 2)
multiome$batch <- strsplit(multiome$sample, "_") %>% sapply("[[", 3)

#remove doublets
multiome <- subset(multiome, subset = DF_doublet == "Singlet")


# process RNA and peaks
DefaultAssay(multiome) <- "RNA"
multiome <- multiome %>% 
  NormalizeData() %>% 
  FindVariableFeatures() %>% 
  ScaleData() %>% 
  RunPCA()

DefaultAssay(multiome) <- "peaks"
multiome <- multiome%>%
  FindTopFeatures(min.cutoff = "q75") %>%
  RunTFIDF() %>%
  RunSVD()

# build a joint neighbor graph and UMAP visualization using both assays
multiome <- multiome %>%
  FindMultiModalNeighbors(
    reduction.list = list("pca", "lsi"), 
    dims.list = list(1:30, 2:30),
    modality.weight.name = "RNA.weight") %>%
  RunUMAP(
    nn.name = "weighted.nn",
    reduction.name = "wnn.umap", 
    reduction.key = "wnnUMAP_") %>%
  FindClusters(
    graph.name = "wsnn", 
    algorithm = 3, 
    resolution = 0.2)

multiome <- RunUMAP(multiome, reduction = 'pca', dims = 1:30, assay = 'RNA', 
                    reduction.name = 'rna.umap', reduction.key = 'rnaUMAP_')
multiome <- RunUMAP(multiome, reduction = 'lsi', dims = 2:30, assay = 'peaks', 
                    reduction.name = 'lsi.umap', reduction.key = 'lsiUMAP_')

# annotate clusters
Idents(multiome) <- "wsnn_res.0.1"
multiome <- RenameIdents(multiome, 
                         "0" = "Oligodendrocytes", 
                         "1" = "Microglia",
                         "2" = "Ependymal",
                         "3" = "Neurons_V",
                         "4" = "Macrophages",
                         "5" = "Astrocytes", "6" = "Astrocytes",
                         "7" = "Neurons_D", "8" = "Neurons_D", "13" = "Neurons_D", "15" = "Neurons_D","11" = "Neurons_D",
                         "9" = "OPCs",
                         "10" = "Perivascular", #Pericytes
                         "14" = "Perivascular", #Lepto
                         "12" = "Endothelial")

Idents(multiome) <- factor(x = Idents(multiome), 
                           levels = c("Neurons_V", "Neurons_D",
                                      "Astrocytes", "Ependymal",
                                      "OPCs", "Oligodendrocytes",
                                      "Microglia", "Macrophages", 
                                      "Perivascular", "Endothelial"))
multiome$cluster_ids <- Idents(multiome)


#clustering modalities (RNA only, ATAC only, multimodal)
#fig 1
multiome <- multiome %>% 
  FindNeighbors(
    reduction = "pca",  
    dims = 1:30,
    graph.name = "RNA.weight") %>%
  RunUMAP(
    dims = 1:30,
    reduction.name = "rna.umap", 
    reduction.key = "rnaUMAP_") %>%
  FindClusters(
    graph.name = "RNA.weight", 
    algorithm = 3, 
    resolution = 1.1)

multiome <- multiome %>%
  FindNeighbors(
    reduction = "lsi",  
    dims = 2:30,
    graph.name = "ATAC.weight") %>%  
  RunUMAP(
    dims = 2:30,
    reduction.name = "lsi.umap", 
    reduction.key = "lsiUMAP_") %>%
  FindClusters(
    graph.name = "ATAC.weight", 
    algorithm = 3, 
    resolution = 1.1)

saveRDS(multiome, "merged_multiome.rds")

