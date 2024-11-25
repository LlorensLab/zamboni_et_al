source("code/source.R")

# enhancer aav reporter assay was mapped with cellranger with ab capture pipeline
cell_ids <- read.table("enhancer_aav/Uninjured/outs/filtered_feature_bc_matrix/barcodes.tsv.gz")
mat_raw <- Read10X_h5("enhancer_aav/Uninjured/outs/raw_feature_bc_matrix.h5")
uninj_obj <- CreateSeuratObject(mat_raw$`Gene Expression`, project = "Uninjured")
uninj_obj[["aavs"]] <- CreateAssayObject(mat_raw$`Antibody Capture`)
uninj_obj <- uninj_obj %>%
  subset(cells = cell_ids$V1) %>%
  NormalizeData()

cell_ids <- read.table("enhancer_aav/Injured/outs/filtered_feature_bc_matrix/barcodes.tsv.gz")
mat_raw <- Read10X_h5("enhancer_aav/Injured/outs/raw_feature_bc_matrix.h5")
inj_obj <- CreateSeuratObject(mat_raw$`Gene Expression`, project = "Injured")
inj_obj[["aavs"]] <- CreateAssayObject(mat_raw$`Antibody Capture`)
inj_obj <- inj_obj %>%
  subset(cells = cell_ids$V1) %>%
  NormalizeData()

# merge enhancer IREN12 (with 2 BCs) 
inj_obj@assays$aavs$data["IREN12",] <- inj_obj@assays$aavs$data["IREN12a",] + inj_obj@assays$aavs$data["IREN12b",]
uninj_obj@assays$aavs$data["IREN12",] <- uninj_obj@assays$aavs$data["IREN12a",] + uninj_obj@assays$aavs$data["IREN12b",]

# merge injured and uninjured dataset
obj <- merge(inj_obj, uninj_obj, merge.data = T) 
obj <- subset(obj, subset = nCount_RNA > 500 & nCount_RNA < 50000)

obj <- obj %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData() %>% 
  RunPCA()
ElbowPlot(obj)

obj <- obj %>%
  RunUMAP(dims = 1:15) %>%
  FindNeighbors(dims = 1:15) %>%
  FindClusters(resolution = 0.5)

FeaturePlot(obj, c("tdTomato-WPRE", "Aldh1l1"), max.cutoff = "q99", order = T)

# only retain astrocytes
astros <- subset(obj, idents = c("0", "1", "2", "3", "4", "7"))
astros <- astros %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData() %>% 
  RunPCA()
ElbowPlot(astros)

astros <- astros %>%
  RunUMAP(dims = 1:8) %>%
  FindNeighbors(dims = 1:8) %>%
  FindClusters(resolution = 0.1)

FeaturePlot(astros, "aavs_IREN5", max.cutoff = "q99", order = T)

saveRDS(astros, "enhancer_aav/astros_aav_obj_240625.rds")


