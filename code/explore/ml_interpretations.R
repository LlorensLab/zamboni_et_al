# interpretations of models output
source("code/source.R")

# extracts qc metrics from chrombpnet json output file
extract_metrics_to_df <- function(json_file_path) {
  json_data <- fromJSON(json_file_path)
  
  spearmanr <- json_data$counts_metrics$peaks$spearmanr
  pearsonr <- json_data$counts_metrics$peaks$pearsonr
  mse <- json_data$counts_metrics$peaks$mse
  median_jsd <- json_data$profile_metrics$peaks$median_jsd
  median_norm_jsd <- json_data$profile_metrics$peaks$median_norm_jsd
  
  # Create a data frame
  df <- data.frame(
    spearmanr = spearmanr,
    pearsonr = pearsonr,
    mse = mse,
    median_jsd = median_jsd,
    median_norm_jsd = median_norm_jsd
  )
  
  return(df)
}

# qc metrics for each validation fold
# Figure 2c
metrics <- list()
path_to_cross_valid <- "chrombpnet_interpret/models/cross_validation/"
for(i in names(celltype_dars)){
  for(j in c("fl0", "fl1", "fl2", "fl3", "fl4")){
    cross_validation <- paste0(i, "_", j)
    json_file_path <- paste0(path_to_cross_valid, j, "/", i, "_", j, "/evaluation/chrombpnet_metrics.json")
    metrics[[cross_validation]] <- extract_metrics_to_df(json_file_path)
    metrics[[cross_validation]]$celltype <- i
    metrics[[cross_validation]]$fold <- j
  }
}

performance <- do.call("rbind", metrics)

summary_performance <- performance %>%
  group_by(celltype) %>%
  summarise(
    mean_correlation = mean(pearsonr),
    std_correlation = sd(pearsonr)
  )

# Barplot
ggplot(summary_performance, aes(x=celltype, y=mean_correlation)) + 
  geom_bar(stat = "identity") + 
  ylim(c(0,1)) +
  geom_errorbar(aes(ymin=mean_correlation, ymax=mean_correlation+std_correlation), width=.2,
                position=position_dodge(.9)) +
  theme_minimal() 


# plot observed and predicted accessibility across cell type models
# Figure 2b
dir_interpret = "chrombpnet_interpret/"
pred_file = "_clusterpeaks_contrib.counts_scores.bw"
clusters <- levels(multiome) # celltypes

predicted <- setNames(lapply(clusters, function(cluster) paste0(dir_interpret, cluster, pred_file)), clusters)

regions <- c("chr2-32852099-32852599", #neuron v
             "chr3-80187273-80187773", #neuron d
             "chr4-40584918-40585418", #astro
             "chr5-90178120-90178620", #epend
             "chr13-90846088-90846588", #opc
             "chr9-109536671-109537171", #oligo
             "chr8-89392714-89393214", #microglia
             "chr1-6464829-6465379", #periph
             "chr1-11842130-11842843", #perivascular
             "chr1-44892462-44892804" #endothelial
)

names(regions) <- clusters

p_coverage <- list()
for(i in names(regions)){
  p_coverage[[i]] <- CoveragePlot(
    multiome,
    regions[i],
    extend.upstream = 400,
    extend.downstream = 400) + 
    scale_fill_manual(values = palette)
}

p_predicted <- list()
for(i in names(regions)){
  p_predicted[[i]] <- BigwigTrack(
    regions[i],
    predicted,
    extend.upstream = 400,
    extend.downstream = 400,
    type = "line",
    y_label = i,
    bigwig.scale = "common") + 
    scale_color_manual(values = rep("black", 10))
}

ggpubr::ggarrange(plotlist = p_predicted, common.legend = T, ncol = 10)
ggpubr::ggarrange(plotlist = p_coverage, common.legend = T, ncol = 10)




## DARs interpret
# check modisco reports
modisco_report <- list()
for(i in levels(multiome)){
  html_file <- paste0("chrombpnet_interpret/", i, "/_modisco_clusterpeaks_report/motifs.html")
  html_content <- read_html(html_file)
  modisco_report[[i]] <- html_content %>% html_nodes("table") %>% .[[1]] %>% html_table()
  modisco_report[[i]]$celltype <- i
  modisco_report[[i]]$match0 <- sub("[.]", "_", modisco_report[[i]]$match0)
  modisco_report[[i]]$match1 <- sub("[.]", "_", modisco_report[[i]]$match1)
  modisco_report[[i]]$match2 <- sub("[.]", "_", modisco_report[[i]]$match2)
  modisco_report[[i]]$tf0 <- sapply(strsplit(modisco_report[[i]]$match0, "_"), "[[", 1) %>% stringr::str_to_title()
  modisco_report[[i]]$tf1 <- sapply(strsplit(modisco_report[[i]]$match1, "_"), "[[", 1) %>% stringr::str_to_title()
  modisco_report[[i]]$tf2 <- sapply(strsplit(modisco_report[[i]]$match2, "_"), "[[", 1) %>% stringr::str_to_title()
  modisco_report[[i]] <- replace_human_genes_with_mouse(modisco_report[[i]], human_to_mouse_dict)
}

# check if TFs are expressed in the cell type of interest (> 1% of cells)
for(i in names(modisco_report)){
  for(j in c("tf0", "tf1", "tf2")){
    tf <- intersect(as.vector(modisco_report[[i]][[j]]), rownames(multiome))
    fc_tf <- FoldChange(multiome, features = tf, ident.1 = i)
    fc_tf$tf_expressed <- ifelse(fc_tf$pct.1 > 0.01, "YES", "NO")
    gene_expression_map <- setNames(fc_tf$tf_expressed, rownames(fc_tf))
    tf_expressed <- paste0(j, "_expressed")
    modisco_report[[i]][[tf_expressed]] <- gene_expression_map[modisco_report[[i]][[j]]] %>% as.vector()
    modisco_report[[i]][[tf_expressed]][is.na(modisco_report[[i]][[tf_expressed]])] <- "Unknown"
  }
}

modisco_report_df_all <- do.call("rbind", modisco_report)
modisco_report_df_all$tf_expressed <- paste0(modisco_report_df_all$tf0_expressed, "_",
                                             modisco_report_df_all$tf1_expressed, "_",
                                             modisco_report_df_all$tf2_expressed)
# Filter out patterns without TF expression
modisco_report_df_all <- modisco_report_df_all %>%
  dplyr::filter(!tf_expressed == "NO_NO_NO")

motifs_freq <- modisco_report_df %>%
  group_by(celltype, motifs_ordered) %>%
  summarise(total_frequency = sum(freq_withpattern)) %>%
  ungroup()

# save modisco_report_df_all and manually curate to order the TF matches consistently across celltypes
# see meta for curated files
modisco_report_adj <- read.csv("modisco_report.csv", sep = ";", row.names = 1) 

# plot motif frequencies
motifs_freq_singles <- modisco_report_adj %>%
  group_by(celltype, tfs) %>%
  summarise(freq_withpattern = sum(freq_withpattern)) %>%
  top_n(20, freq_withpattern) %>%
  ungroup()
motifs_freq_singles$round_freq <- round(motifs_freq_singles$freq_withpattern, 3)

motifs_freq_wide <- motifs_freq_singles[, c("celltype", "tfs", "round_freq")] %>%
  tidyr::spread(key = tfs, value = round_freq, fill = 0)
celltype <- motifs_freq_wide$celltype
motifs_freq_wide <- motifs_freq_wide[,c(2:ncol(motifs_freq_wide))]
rownames(motifs_freq_wide) <- celltype

fun_color_range <- colorRampPalette(c("#2D558E", "grey99","#8A267C"))  # Create color generating function
hm_colors <- fun_color_range(100) 

# Figure 2d
pheatmap::pheatmap(t(motifs_freq_wide), cluster_rows=T, cluster_cols=T,
                   scale="row", clustering_method="ward.D2",
                   fontsize=6, color = hm_colors, 
                   border_color = NA)

# plot tf expression
tfs <- unique(unlist(strsplit(modisco_report_adj$motifs, "/")))

# Figure S5e
avg_rna_celltypes <- AverageExpression(multiome, assays = "RNA", return.seurat = T, group.by = "cluster_ids")
agg_rna <- avg_rna_celltypes@assays$RNA@data[tfs,] %>% as.matrix()
pheatmap::pheatmap(agg_rna, 
                   cluster_rows=T, cluster_cols=TRUE,
                   show_colnames = T,
                   scale="row", clustering_method="ward.D2",
                   fontsize=6, color = viridis(80, option = "inferno"), 
                   breaks = seq(-4,4, 0.1)) 


# plot number of motifs per regions and proportion of regions with/without motifs
motifs_per_region <- list()
for(i in unique(multiome$cluster_ids)){
  file_path = paste0("chrombpnet_out/interpret/", i, "_motif_counts_per_region.csv") #see jupyter note to create file
  motifs_per_region[[i]] <- read.csv(file_path, header = T, sep = ",")
  motifs_per_region[[i]]$celltype <- i
}

zero_patterns <- NULL
for(i in names(celltype_dars)){
  diff <- peaks$peak_n[peaks$celltype == i] - nrow(motifs_per_region[[i]])
  zero_patterns <- c(zero_patterns, diff)
}

peaks$zero_pattern <- zero_patterns 
peaks$percent_with <- (peaks$peak_n - peaks$zero_pattern)/peaks$peak_n *100

peaks$celltype <- factor(peaks$celltype, levels = names(celltype_dars))

# Figure S5a
p1 <- ggplot(peaks, aes(x=celltype, y=peak_n)) + 
  geom_bar(stat = "identity") +
  ylim(0,165000) +
  theme_minimal() 
p2 <- ggplot(peaks, aes(x=celltype, y=zero_pattern)) + 
  geom_bar(stat = "identity") +
  ylim(0,165000) +
  theme_minimal() 
p1 / p2

# plot motifs per region - Figure S5b
freq_df <- do.call("rbind", motifs_per_region)
table_df <- table(freq_df$celltype, freq_df$Count) %>% as.data.frame()
table_df$patterns <- ifelse(table_df$Var2 %in% c(4, 5, 6, 7, 8, 9), "4+", table_df$Var2)

ggplot(table_df, aes(x=patterns, y=Freq, fill = Var1)) + 
  geom_bar(stat="identity", width=0.9, position = "dodge") +
  scale_fill_manual(values = palette)

n_motifs <- NULL
for(i in names(celltype_dars)){n_motifs <- c(n_motifs, nrow(modisco_report[[i]]))}

n_motifs_df <- data.frame(clusters = as.factor(names(celltype_dars)),
                          n_motifs = n_motifs
)

ggplot(n_motifs_df, aes(x=reorder(clusters, names(celltype_dars)), y=n_motifs)) + 
  geom_bar(stat = "identity") +
  theme_minimal() 


## chromvar comparison - Figure S5d
# add motif information
multiome <- AddMotifs(
  object = multiome,
  genome = BSgenome.Mmusculus.UCSC.mm10,
  pfm = pfm
)

DefaultAssay(multiome) <- "peaks"
Idents(multiome) <- "cluster_ids"
multiome <- RegionStats(multiome, genome = BSgenome.Mmusculus.UCSC.mm10)
enriched.motifs <- list()
for(i in levels(multiome)){
  # match the overall GC content in the peak set
  meta.feature <- GetAssayData(multiome, assay = "peaks", layer = "meta.features")
  open.peaks <- AccessiblePeaks(multiome, idents = i)
  
  peaks.matched <- MatchRegionStats(
    meta.feature = meta.feature[open.peaks, ],
    query.feature = meta.feature[unique(celltype_dars[[i]]$peak), ],
    n = 50000
  )
  
  # test enrichment
  enriched.motifs[[i]] <- FindMotifs(
    object = multiome,
    features = unique(celltype_dars[[i]]$peak),
    background=peaks.matched
  )
  
}

DefaultAssay(multiome) <- "RNA"
for(i in names(celltype_dars)){
  # test enrichment
  enriched.motifs[[i]]$tfs <- sapply(strsplit(enriched.motifs[[i]][["motif.name"]], "[(]"), "[[", 1) %>%
    stringr::str_to_title() %>%
    replace_human_genes_with_mouse(human_to_mouse_dict)
  
}

expressed <- list()
for(i in names(celltype_dars)){
  significant_motifs <- enriched.motifs[[i]] %>%
    dplyr::filter(p.adjust < 0.05)
  tfs <- intersect(significant_motifs$tfs, rownames(multiome)) %>% as.character()
  expressed[[i]] <- FoldChange(multiome, features = tfs, ident.1 = i)
  expressed[[i]]$expressed <- ifelse(expressed[[i]]$pct.1 > 0.01, "YES", "NO")
  expressed[[i]]$gene <- rownames(expressed[[i]])
  print(table(expressed[[i]]$expressed))
}

expressed_df <- do.call("rbind", expressed)
expressed_df$expressed %>% table()

expressed_df_chromvar <- expressed_df[expressed_df$expressed == "YES",]
expressed_df_modisco <- modisco_report_adj[modisco_report_adj$tf_expressed == "YES",]
intersect(unique(expressed_df_chromvar$gene), unique(expressed_df_modisco$tfs)) %>% length()


