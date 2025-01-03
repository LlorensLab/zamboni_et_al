# interpretation injury
source("code/source.R")

celltypes_injury <- c("Astrocytes_I", "Astrocytes_U",
                      "Ependymal_I", "Ependymal_U",
                      "Microglia_I", "Microglia_U",
                      "OPCs_I", "OPCs_U",
                      "Oligodendrocytes_I", "Oligodendrocytes_U")

motifs_injury_per_region <- list()
injury_beds <- list()
modisco_injury_report <- list()
for(i in celltypes_injury){
  file_path = paste0("chrombpnet_interpret_injury/", i, "_injury_motif_counts_per_region.csv") # see jupyter to make file
  motifs_injury_per_region[[i]] <- read.csv(file_path, header = T, sep = ",")
  motifs_injury_per_region[[i]]$celltype <- i
  file_path = paste0("chrombpnet_interpret_injury/", i, "_injury_contrib.interpreted_regions.bed") # chrombpnet output
  injury_beds[[i]] <- read.delim(file_path, sep = "\t", header = F)
  
  html_file <- paste0("chrombpnet_interpret_injury/", i, "_injury_modisco_report/motifs.html") # chrombpnet output
  html_content <- read_html(html_file)
  modisco_injury_report[[i]] <- html_content %>% html_nodes("table") %>% .[[1]] %>% html_table()
  modisco_injury_report[[i]]$celltype <- i
  modisco_injury_report[[i]]$match0 <- sub("[.]", "_", modisco_injury_report[[i]]$match0)
  modisco_injury_report[[i]]$match1 <- sub("[.]", "_", modisco_injury_report[[i]]$match1)
  modisco_injury_report[[i]]$match2 <- sub("[.]", "_", modisco_injury_report[[i]]$match2)
  modisco_injury_report[[i]]$tf0 <- sapply(strsplit(modisco_injury_report[[i]]$match0, "_"), "[[", 1) %>% stringr::str_to_title()
  modisco_injury_report[[i]]$tf1 <- sapply(strsplit(modisco_injury_report[[i]]$match1, "_"), "[[", 1) %>% stringr::str_to_title()
  modisco_injury_report[[i]]$tf2 <- sapply(strsplit(modisco_injury_report[[i]]$match2, "_"), "[[", 1) %>% stringr::str_to_title()
  modisco_injury_report[[i]] <- replace_human_genes_with_mouse(modisco_injury_report[[i]], human_to_mouse_dict)
  modisco_injury_report[[i]]$total_peaks <- nrow(injury_beds[[i]])
  modisco_injury_report[[i]]$peaks_withpattern <- nrow(motifs_injury_per_region[[i]])
  
}

DefaultAssay(multiome) <- "RNA"
Idents(multiome) <- "cluster_ids_injury"
for(i in names(modisco_injury_report)){
  for(j in c("tf0", "tf1", "tf2")){
    tf <- intersect(as.vector(modisco_injury_report[[i]][[j]]), rownames(multiome))
    fc_tf <- FoldChange(multiome, features = tf, ident.1 = i)
    fc_tf$tf_expressed <- ifelse(fc_tf$pct.1 > 0.01, "YES", "NO")
    gene_expression_map <- setNames(fc_tf$tf_expressed, rownames(fc_tf))
    tf_expressed <- paste0(j, "_expressed")
    modisco_injury_report[[i]][[tf_expressed]] <- gene_expression_map[modisco_injury_report[[i]][[j]]] %>% as.vector()
    modisco_injury_report[[i]][[tf_expressed]][is.na(modisco_injury_report[[i]][[tf_expressed]])] <- "Unknown"
    
  }
}

modisco_injury_report_df_all <- do.call("rbind", modisco_injury_report)
modisco_injury_report_df_all$tf_expressed <- paste0(modisco_injury_report_df_all$tf0_expressed, "_",
                                                    modisco_injury_report_df_all$tf1_expressed, "_",
                                                    modisco_injury_report_df_all$tf2_expressed)

modisco_injury_report_df <- modisco_injury_report_df_all %>%
  dplyr::filter(!tf_expressed == "NO_NO_NO")


### find injury-enriched patterns
gr_list <- list()
for(i in names(injury_beds)){
  injury_beds[[i]]$peak <- paste0(injury_beds[[i]]$V1, "-", injury_beds[[i]]$V2, "-", injury_beds[[i]]$V3)
  
  peak_set <- injury_beds[[i]][,1:3]
  colnames(peak_set) <- c("chr", "start", "end")
  gr_list[[i]] <- makeGRangesFromDataFrame(peak_set)
}

# get patterns for each peak
modisco_injury_report_pos <- modisco_injury_report_df_all[grep("pos_",
                                                               modisco_injury_report_df_all$pattern),]
modisco_injury_report_pos$Pattern <- sapply(strsplit(modisco_injury_report_pos$pattern, "[.]"), "[[", 2)
modisco_injury_report_pos$motif <- paste0(modisco_injury_report_pos$tf0, "_",
                                          modisco_injury_report_pos$tf1, "_",
                                          modisco_injury_report_pos$tf2)
patterns_df <- list()
for(i in celltype){
  csv_path <- paste0("chrombpnet_interpret_injury/", i, "_example_idx.csv") #see jupyter notebook to create file
  csv <- read.csv(csv_path, header = T)
  modisco_injury_celltype <- modisco_injury_report_pos %>% dplyr::filter(celltype == i)
  csv_motif <- merge(csv, modisco_injury_celltype, by = "Pattern")
  
  # Summarize the data by concatenating patterns for each ID
  patterns_df[[i]] <- csv_motif %>%
    group_by(Example_idx) %>%
    summarise(patterns = paste(Pattern, collapse = ", "),
              motifs =  paste(motif, collapse = ", ")) %>%
    ungroup()
  
  patterns_df[[i]]$Example_idx <-  patterns_df[[i]]$Example_idx +1
}

# look at motifs enriched in the DARs for each celltype and injury condition
# get dars 
marker_peaks_injury <- list()
DefaultAssay(multiome) <- "peaks"
Idents(multiome) <- "cluster_ids_injury"
for(i in celltypes_injury){
  marker_peaks_injury[[i]] <- FindMarkers(multiome, ident.1 = i, 
                                          only.pos = T, 
                                          logfc.threshold = 0.5)
  
}
saveRDS(marker_peaks_injury, "markers/markers_injury_peaks.rds")
marker_peaks_injury <- readRDS("markers/markers_injury_peaks.rds")

# find overlaps between computed dars and peaks used for model interpretation
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

tables_patterns <- list()
for(i in celltypes_injury){
  overlaps <- subsetByOverlaps(gr_list[[i]], gr_dars[[i]])
  overlaps_i <- injury_beds[[i]][injury_beds[[i]]$peak %in% 
                                   paste0(overlaps@seqnames, "-", overlaps@ranges),]
  
  tables_patterns[[i]] <- table(patterns_df[[i]]$motifs[
    patterns_df[[i]]$Example_idx %in% as.numeric(rownames(overlaps_i))]) %>%
    as.data.frame()
  
}

# plot enriched motifs for each cell type and injury condition
# Figure 5c, S8
for(i in c("Astrocytes", "Ependymal", "Microglia", "OPCs", "Oligodendrocytes")){
  cell_injury <- paste0(i, "_I")
  cell_uninjured <- paste0(i, "_U")
  
  top_i_patterns <- tables_patterns[[cell_injury]] %>%
    top_n(30, Freq)
  p1 <- ggplot(top_i_patterns, aes(x=reorder(Var1, Freq), y=Freq)) + 
    geom_bar(stat = "identity") +
    coord_flip() +
    theme_minimal() 
  
  top_u_patterns <- tables_patterns[[cell_uninjured]] %>%
    top_n(30, Freq)
  p2 <- ggplot(top_u_patterns, aes(x=reorder(Var1, Freq), y=Freq)) + 
    geom_bar(stat = "identity") +
    coord_flip() +
    theme_minimal() 
  
  print(p1 | p2)
  
}
