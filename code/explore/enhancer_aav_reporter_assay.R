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

irens <- paste0("aavs_IREN", c(5, 8, 11, 12, 24))
FeaturePlot(astros, irens, max.cutoff = "q99", order = T)

# count detected BCs per cell
Idents(astros) <- "orig.ident"
injury <- subset(astros, idents = "Injured")
meta_aavs <- injury@assays$aavs@counts[c(irens, "A1", "A5"),] %>% t() %>% as.data.frame()
meta_aavs_bin <- meta_aavs %>%
  mutate(across(everything(), ~ ifelse(. != 0, 1, 0)))
colnames(meta_aavs_bin) <- c(irens, "A1", "A5")
meta_aavs_bin <- meta_aavs_bin %>%
  mutate(RowSum = rowSums(across(everything()))) %>%
  mutate(RowSum_M = rowSums(across(1:5)))

table_df <- table(meta_aavs_bin$RowSum) %>% as.data.frame()
ggplot(table_df, aes(x=Var1, y=Freq)) + 
  geom_bar(stat = "identity")


# plot model interpretation (observed vs predicted accessibility)
# load profiles observed and predicted by each model (astrocytes, reactive astrocytes and microglia)
observed_m5 <- list(
  Astrocytes = "chrombpnet_interpret/models/Astrocytes_fl0/auxiliary/data_unstranded.bw",
  Reactive_astrocytes = "chrombpnet_interpret/models/Injury_Astrocytes_fl0/auxiliary/data_unstranded.bw",
  Microglia = "chrombpnet_interpret/models/Microglia_fl0/auxiliary/data_unstranded.bw")

predicted_m5 <- list(
  Astrocytes = "enhancer_virus_astrocytes_bigwig_chrombpnet.bw",
  Reactive_astrocytes = "enhancer_virus_bigwig_astrocytes_injury_chrombpnet.bw",
  Microglia = "enhancer_virus_microglia_bigwig_chrombpnet.bw")

coordinates <- read.table("enhancer_aav/Enhancer_coordinates_BCs.csv", sep = ";")
enhancers = list(coordinates$V1)
names(enhancers) <- irens

p_obs <- list()
p_pred <- list()
for(i in names(enhancers)){
  p_obs[[i]] <- BigwigTrack(
    enhancers[[i]], observed_m5,
    extend.upstream = 1000, extend.downstream = 1000,
    type = "coverage", y_label = i, bigwig.scale = "common") + 
    scale_fill_manual(values = c(palette[3], palette[4], palette[7]))
  
  p_pred[[i]] <- BigwigTrack(
    enhancers[[i]], predicted_m5,
    extend.upstream = 1000, extend.downstream = 1000,
    type = "line", y_label = i, bigwig.scale = "common") + 
    scale_color_manual(values = rep("black", 10))
  
}

# plot observed and predicted profiles
ggpubr::ggarrange(plotlist = p_obs, common.legend = T, ncol = 5)
ggpubr::ggarrange(plotlist = p_pred, common.legend = T, ncol = 5)



## interpret enhancer viruses
coordinates <- read.table("enhancer_aav/Enhancer_coordinates_BCs.csv", sep = ";")
chr <- sapply(strsplit(coordinates$V3, "-"), "[", 1) 
start <- sapply(strsplit(coordinates$V3, "-"), "[", 2) 
end <- sapply(strsplit(coordinates$V3, "-"), "[", 3) 
summit <- round(coordinates$V4/2, 0)
bed <- cbind(chr, start, end, 
             rep("-", 23), rep("-", 23), rep("-", 23),
             rep("-", 23), rep("-", 23), rep("-", 23),
             summit)
bed_file <- paste0("chrombpnet_interpret/enhancer_virus.bed")
write.table(bed, bed_file,row.names = F,col.names = F, sep="\t", quote=FALSE)



saveRDS(astros, "enhancer_aav/astros_aav_obj_240625.rds")


