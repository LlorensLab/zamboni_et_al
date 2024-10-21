source("../source.R")

object <- readRDS("unmerged_obj_multiome_clean.rds")

# merge objects
combined <- merge(
  x = object[[1]],
  y = c(object[2:21])
)

rm(object)

combined$experiment <- strsplit(combined$sample, "_") %>% sapply("[[", 1)
combined$condition <- strsplit(combined$sample, "_") %>% sapply("[[", 2)
combined$batch <- strsplit(combined$sample, "_") %>% sapply("[[", 3)

#remove doublets
combined <- subset(combined, subset = DF_doublet == "Singlet")

# process RNA and mpeaks
DefaultAssay(combined) <- "RNA"
combined <- combined %>% 
  NormalizeData() %>% 
  FindVariableFeatures() %>% 
  ScaleData() %>% 
  RunPCA()

DefaultAssay(combined) <- "cpeaks"
combined <- combined%>%
  FindTopFeatures(min.cutoff = "q75") %>%
  RunTFIDF() %>%
  RunSVD()

# build a joint neighbor graph and UMAP visualization using both assays
combined <- combined %>%
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

combined <- RunUMAP(combined, reduction = 'pca', dims = 1:30, assay = 'RNA', 
                    reduction.name = 'rna.umap', reduction.key = 'rnaUMAP_')
combined <- RunUMAP(combined, reduction = 'lsi', dims = 2:30, assay = 'peaks', 
                    reduction.name = 'lsi.umap', reduction.key = 'lsiUMAP_')
combined <- FindClusters(combined,
                         graph.name = "wsnn", 
                         algorithm = 3, 
                         resolution = 0.2)

