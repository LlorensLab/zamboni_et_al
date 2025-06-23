# code for heatmaps in Figure S11a
multiome <- readRDS("../sci_mus_multiome_only/multiome_obj_clean_240326.rds")
Idents(multiome) <- "cluster_ids_injury"
DefaultAssay(multiome) <- "peaks"

# load eregulon computed with standard SCENI+ pipeline
eregulon <- read.csv("~/Downloads/scenic/eRegulons_direct.tsv", sep = "\t")
eregulon$Region <- gsub(":", "-", eregulon$Region)

# load precomputed injury dars
marker_peaks_injury <- readRDS("markers/markers_injury_peaks.rds")
celltypes_injury <- names(marker_peaks_injury)

# find overlaps between precomputed injury dars and peaks used for model interpretation
gr_dars <- list()
for(i in celltypes_injury){
  marker_peaks_injury[[i]]$peak <- rownames(marker_peaks_injury[[i]])
  sig_markers <- marker_peaks_injury[[i]] %>%
    dplyr::filter(p_val < 0.05 & pct.2 < 0.1 & avg_log2FC > 1) 
  
  gr_dars[[i]] <- GRanges(seqnames = sapply(strsplit(sig_markers$peak, "-"), "[[", 1), 
                          ranges = paste0(
                            sapply(strsplit(sig_markers$peak, "-"), "[[", 2), "-",
                            sapply(strsplit(sig_markers$peak, "-"), "[[", 3)))
}

binary_matrices <- list()
for(i in celltypes_injury){
  genes_expressed <- FoldChange(multiome, ident.1 = i, assay = "RNA") %>%
    filter(pct.1 > 0.01)
  eregulon_celltype <- eregulon %>%
    filter(TF %in% rownames(genes_expressed) & 
             Gene %in% rownames(genes_expressed) & 
             regulation > 0)

  # find overlaps between regions in eregulons and precomputed injury dars
  overlaps <- subsetByOverlaps(GRanges(eregulon_celltype$Region), gr_dars[[i]])
  eregulons_peaks <- paste0(overlaps@seqnames, ":", overlaps@ranges)
  celltype_eregulons <- eregulon_celltype[eregulon_celltype$Region %in% eregulons_peaks,]
  
  eregulons_long <- celltype_eregulons %>%
    select(TF, Region) %>%
    distinct(Region, TF, .keep_all = TRUE) %>%  # Remove exact duplicates
    group_by(Region, TF) %>%
    summarise(Present = 1L, .groups = "drop")  # Ensure single value per combination
  
  # Create binary matrix
  binary_matrix <- eregulons_long %>%
    pivot_wider(names_from = TF, values_from = Present, values_fill = 0) %>%
    column_to_rownames("Region")
  
  # Order TFs (columns) by frequency
  tfbs_order <- colSums(binary_matrix, na.rm = TRUE) %>%
    sort(decreasing = TRUE) %>%
    names()
  
  # Reorder columns (TFs) by frequency
  binary_matrix <- binary_matrix[, tfbs_order]
  binary_matrix <- binary_matrix[rowSums(binary_matrix[])>0, ]
  
  # Function to order rows based on TF presence for each TF in order
  for (tf in tfbs_order) {
    binary_matrix <- binary_matrix[order(-binary_matrix[, tf]), ]  # Sort rows by presence of each TF
  }
  
  binary_matrices[[i]] <- binary_matrix
  
}

cols <- c("#8b167c", "#ff595e", "#ffcc00", "#f89f00", "#2C558F")
names(cols) <- c("Astrocytes", "Ependymal", "OPCs", "Oligodendrocytes", "Microglia")
hm <- list()
for(i in names(binary_matrices)){
  celltype <- strsplit(i, "_")[[1]][1]
  hm[[i]] <- pheatmap(binary_matrices[[i]][,1:15], cluster_rows = F, cluster_cols = F, color = c("white", cols[celltype]))
}

ggpubr::ggarrange(plotlist = hm[1:10], ncol = 4, nrow = 4)
