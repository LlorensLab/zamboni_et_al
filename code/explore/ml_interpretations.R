# interpretations of models output

source("code/source.R")
library(jsonlite)

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


