source("source.R")

multiome <- readRDS("objects/multiome_clean_240925.rds")
Idents(multiome) <- "cluster_ids"

# prepare fragment file for bias model training 
#(take batch 1, encompassing unsorted data from all timepoints)
SplitFragments(
  multiome,
  assay = "peaks",
  group.by = "batch",
  idents = "1",
  outdir = "ML/celltype_fragments/",
  append = TRUE,
  buffer_length = 256L,
  verbose = TRUE, 
  file.suffix = "_input_bias"
)

# split fragment by celltype
SplitFragments(
  multiome,
  assay = "peaks",
  group.by = "cluster_ids",
  outdir = "ML/celltype_fragments/",
  append = TRUE,
  buffer_length = 256L,
  verbose = TRUE
)

# split fragment by celltype and injury
multiome$injury <- ifelse(multiome$condition == "Uninj", "U", "I")
table(multiome$injury)
multiome$cluster_ids_injury <- paste0(multiome$cluster_ids, "_", multiome$injury)
Idents(multiome) <- "cluster_ids_injury"
Idents(multiome) %>% unique()

SplitFragments(
  multiome,
  assay = "peaks",
  outdir = "ML/celltype_fragments_injury/",
  append = TRUE,
  buffer_length = 256L,
  verbose = TRUE
)


# prepare bed files for chrombpnet interpretation (e.g., cell type- and injury-specific regions)
celltype_injury_dars <- readRDS("markers/markers_Uvsothers_peaks_all.rds")
for(i in names(celltype_injury_dars)){
  celltype_injury_dars[[i]]$peak <- sapply(strsplit(rownames(celltype_injury_dars[[i]]), "[.]"), "[[", 2)
  celltype_injury_dars[[i]] <- celltype_injury_dars[[i]] %>%
    dplyr::filter(p_val < 0.05 & avg_log2FC > 0.5)
}

for(i in names(celltype_injury_dars)){
  chr <- sapply(strsplit(celltype_injury_dars[[i]]$peak, "-"), "[", 1) 
  start <- sapply(strsplit(celltype_injury_dars[[i]]$peak, "-"), "[", 2) 
  end <- sapply(strsplit(celltype_injury_dars[[i]]$peak, "-"), "[", 3) 
  bed <- cbind(chr, start, end)
  bed_file <- paste0(i, "_I_dars.bed")
  write.table(bed, bed_file,row.names = F,col.names = F, sep="\t", quote=FALSE)
}
