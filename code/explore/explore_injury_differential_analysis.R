# code for exploring differential analysis in injury:
# clustering of genes and peaks based on activation dynamics
# shared glial gene expression program
# AP1 transcription factor expression dynamics


#clustering of DEGs and DARs
# injury-induced genes and peaks are clustered using TCseq into modules with the same dynamics
# Figure 3
source("../source.R")
multiome <- readRDS("merged_multiome.rds")

glia <- c("Microglia", "OPCs", "Oligodendrocytes", "Ependymal", "Astrocytes")

markers_Uvsothers_rna <- readRDS("markers_Uvsothers_rna_all.rds")
markers_Uvsothers_peaks <- readRDS("markers_Uvsothers_peaks_all.rds")

#RNA
p_tclust <- list()
tclust <- list()
avg_celltype <- list()
for(i in glia){
  celltype <- subset(multiome, idents = i)
  Idents(celltype) <- "condition"
  DefaultAssay(celltype) <- "RNA"
  
  markers <- markers_Uvsothers_rna[[i]]
  avg_celltype[[i]] <- AverageExpression(celltype, return.seurat = T, assays = "RNA")
  genes <- intersect(markers$gene, rownames(avg_celltype[[i]]@assays$RNA@counts))
  
  tclust[[i]] <- timeclust(x = as.matrix(avg_celltype[[i]]@assays$RNA@counts[genes,]), algo="cm", k=6, standardize =T) #i use cmeans for fuzzy clustering as algorithm as it gives you membership scores
  p_tclust[[i]] <- timeclustplot(tclust[[i]], value = paste0("z-score(mean(vst))"), cols = 3)
}

### go analysis
#####
go <- list()
go_celltype <- list()
for(i in c("Astrocytes", "Ependymal", "Microglia", "OPCs", "Oligodendrocytes")){
  tclust_celltype <- tclust[[i]]
  
  for(j in 1:6){
    cluster_genes <- names(tclust_celltype@cluster)[tclust_celltype@cluster %in% j]
    gostres <- gost(query = cluster_genes,
                    organism = "mmusculus")
    gostres$result$log_p <- -log10(gostres$result$p_value)
    gostres$result$cluster <- j
    
    go_celltype[[j]] <- gostres$result
  }
  go[[i]] <- do.call("rbind", go_celltype)
  go[[i]]$celltype <- i
}
go_df <- do.call(rbind, go)

# filter out general terms (based on the number of genes assigned to the GO hit)
top_go <- list()
for(i in names(go)){
  top_go[[i]] <- go[[i]] %>% 
    dplyr::filter(source == "GO:BP") %>%
    dplyr::filter(term_size < 2000) %>%
    group_by(cluster) %>%
    top_n(10, log_p)
  top_go[[i]][["parents"]] <- NULL
  write.csv(as.data.frame(top_go[[i]]), paste0("top_go_tclust_", i, ".csv"))
}



## TCseq on ATAC
p_tclust <- list()
tclust <- list()
avg_celltype <- list()
for(i in glia){
  celltype <- subset(multiome, idents = i)
  Idents(celltype) <- "condition"
  DefaultAssay(celltype) <- "peaks"
  
  markers <- markers_Uvsothers_peaks[[i]]
  avg_celltype[[i]] <- AverageExpression(celltype, return.seurat = T, assays = "peaks")
  peak <- intersect(markers$peak, rownames(avg_celltype[[i]]@assays$peaks@counts))
  
  tclust[[i]] <- timeclust(x = as.matrix(avg_celltype[[i]]@assays$peaks@counts[peak,]), algo="cm", k=6, standardize =T) #i use cmeans for fuzzy clustering as algorithm as it gives you membership scores
  p_tclust[[i]] <- timeclustplot(tclust[[i]], value = paste0("z-score(mean(vst))"), cols = 3)
}


# group DARs in modules
astro_modules <- list(
  uninjured = c(
    names(tclust$Astrocytes@cluster)[tclust$Astrocytes@cluster %in% 4],
    names(tclust$Astrocytes@cluster)[tclust$Astrocytes@cluster %in% 1]
  ),
  acute = c(
    names(tclust$Astrocytes@cluster)[tclust$Astrocytes@cluster %in% 3],
    names(tclust$Astrocytes@cluster)[tclust$Astrocytes@cluster %in% 2],
    names(tclust$Astrocytes@cluster)[tclust$Astrocytes@cluster %in% 5]
  ),
  chronic = c(
    names(tclust$Astrocytes@cluster)[tclust$Astrocytes@cluster %in% 6]
  )
)

epend_modules <- list(
  uninjured = c(
    names(tclust$Ependymal@cluster)[tclust$Ependymal@cluster %in% 4]
  ),
  acute = c(
    names(tclust$Ependymal@cluster)[tclust$Ependymal@cluster %in% 1],
    names(tclust$Ependymal@cluster)[tclust$Ependymal@cluster %in% 6],
    names(tclust$Ependymal@cluster)[tclust$Ependymal@cluster %in% 5]
  ),
  chronic = c(
    names(tclust$Ependymal@cluster)[tclust$Ependymal@cluster %in% 3],
    names(tclust$Ependymal@cluster)[tclust$Ependymal@cluster %in% 2]
  )
)
ol_modules <- list(
  uninjured = c(
    names(tclust$OPCs@cluster)[tclust$OPCs@cluster %in% 4],
    names(tclust$Oligodendrocytes@cluster)[tclust$Oligodendrocytes@cluster %in% 6]
    
  ),
  acute = c(
    names(tclust$OPCs@cluster)[tclust$OPCs@cluster %in% 6],
    names(tclust$OPCs@cluster)[tclust$OPCs@cluster %in% 1],
    names(tclust$OPCs@cluster)[tclust$OPCs@cluster %in% 5],
    names(tclust$Oligodendrocytes@cluster)[tclust$Oligodendrocytes@cluster %in% 1],
    names(tclust$Oligodendrocytes@cluster)[tclust$Oligodendrocytes@cluster %in% 4],
    names(tclust$Oligodendrocytes@cluster)[tclust$Oligodendrocytes@cluster %in% 5],
    names(tclust$Oligodendrocytes@cluster)[tclust$Oligodendrocytes@cluster %in% 3]
    
  ),
  chronic = c(
    names(tclust$OPCs@cluster)[tclust$OPCs@cluster %in% 2],
    names(tclust$OPCs@cluster)[tclust$OPCs@cluster %in% 3],
    names(tclust$Oligodendrocytes@cluster)[tclust$Oligodendrocytes@cluster %in% 2]
  )
)
mg_modules <- list(
  uninjured = c(
    names(tclust$Microglia@cluster)[tclust$Microglia@cluster %in% 2],
    names(tclust$Microglia@cluster)[tclust$Microglia@cluster %in% 4]
  ),
  acute = c(
    names(tclust$Microglia@cluster)[tclust$Microglia@cluster %in% 3],
    names(tclust$Microglia@cluster)[tclust$Microglia@cluster %in% 5]
  ),
  chronic = c(
    names(tclust$Microglia@cluster)[tclust$Microglia@cluster %in% 6],
    names(tclust$Microglia@cluster)[tclust$Microglia@cluster %in% 1]
  )
  
)
modules_peaks <- list(astro_modules, epend_modules, ol_modules, mg_modules)
names(modules_peaks) <- c("astro_modules", "epend_modules", "ol_modules", "mg_modules")

# get enriched motifs for each module and cell type
enriched_motifs <- list()
enriched_celltype <- list()
DefaultAssay(multiome) <- "peaks"

pfm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE)
)
multiome <- AddMotifs(
  object = multiome,
  genome = BSgenome.Mmusculus.UCSC.mm10,
  pfm = pfm
)

for(j in names(modules_peaks)){
  celltype_module <- modules_peaks[[j]]
  
  for(i in c("uninjured", "acute", "chronic")){
    enriched_celltype[[i]] <- FindMotifs(
      object = multiome,
      features = celltype_module[[i]])
    enriched_celltype[[i]]$module <- i
  }
  enriched_motifs[[j]] <- do.call("rbind", enriched_celltype)
}

## shared glial injury programs
# upset plot to extract injury-enriched genes shared among glia cells
######
# Figure 4, Figure S7
list_markers <- list()
for(i in c("Astrocytes", "Ependymal", "Microglia", "Oligodendrocytes", "OPCs")){
  markers <- markers_Uvsothers_rna_all[[i]] %>%
    filter(abs(avg_log2FC) > 0.5 & p_val_adj < 0.05) %>%
    select(gene) 
  markers <- markers$gene %>% unique() %>% as.character()
  
  list_markers[[i]] <- markers
  
}
upset(fromList(list_markers), order.by = "freq")

# List of intersections 
df_int <- lapply(df1$gene,function(x){
  intersection <- df1 %>% 
    dplyr::filter(gene==x) %>% 
    arrange(path) %>% 
    pull("path") %>% 
    paste0(collapse = "|")
  data.frame(gene = x,int = intersection)
}) %>% 
  bind_rows()

pan_glia_injury <- df_int[df_int$int %in% "Astrocytes|Ependymal|Microglia|OPCs|Oligodendrocytes",]

#compute average expression for glia cells
glia_clusters <- c("Astrocytes", "Ependymal", "OPCs", "Oligodendrocytes", "Microglia")
glia <- subset(multiome, idents = glia_clusters)
avg_rna_clusters <- AverageExpression(glia, assays = "RNA", group.by = "cluster_ids", return.seurat = T)

DoHeatmap(avg_rna_clusters, pan_glia_injury$gene, draw.lines = F) +
  scale_fill_viridis(option = "inferno") +
  NoLegend() 

#get peaks linked to genes shared among glia cells
linkedpeaks_df <- do.call("rbind", linkedpeaks)
linked_panglia <- linkedpeaks_df %>%
  dplyr::filter(gene %in% pan_glia_injury$gene)

#plot average accessibility
avg_peak_clusters <- AverageExpression(glia, assays = "peaks", group.by = "cluster_ids", return.seurat = T)
DoHeatmap(avg_peaks_clusters, linked_panglia$peak, draw.lines = F) + 
  scale_fill_viridis(option = "mako") +
  NoLegend()

# pairwise comparison
pairwise_glia_injury <- df_int[df_int$int %in% c("Astrocytes|Ependymal",
                                                 "Astrocytes|Microglia",
                                                 "Astrocytes|OPCs",
                                                 "Ependymal|Microglia",
                                                 "Ependymal|OPCs"),]

DoHeatmap(avg_rna_clusters, pairwise_glia_injury$gene, draw.lines = F) +
  scale_fill_viridis(option = "inferno") +
  NoLegend()

# explore functional terms on DAVID



## explore dynamic expression of AP1 transcription factors
# Figure 4, S8
genes <- "Fos|Jun|Atf|Jdp|Batf|Maf"
all_ap1 <- rownames(multiome@assays$RNA)[grep(genes, rownames(multiome@assays$RNA))]

DefaultAssay(multiome) <- "RNA"
pdf("all_AP1_genes.pdf", width = 10, height = 10)
for(i in all_ap1){
  p1 <- FeaturePlot(multiome, i, order = T, split.by = "injury")
  p2 <- VlnPlot(multiome, i, group.by = "cluster_ids_timepoint") + NoLegend()
  print(p1 / p2)
}
dev.off()

Idents(glia) <- "cluster_ids_timepoint"
timepoints <- c("U", "1dpi", "3dpi", "7dpi", "28dpi")
combinations <- expand.grid(glia_clusters, timepoints)
levels(glia) <- paste0(combinations$Var1, "_", combinations$Var2)

avg_rna <- AverageExpression(glia, return.seurat = T, assays = "RNA")
df_ap1 <- avg_rna@assays$RNA$data[all_ap1,] %>% t() %>% as.data.frame()
df_ap1$celltype <- sapply(strsplit(rownames(df_ap1), "-"), "[[", 1)
df_ap1$timepoint <- sapply(strsplit(rownames(df_ap1), "-"), "[[", 2) %>%
  factor(levels = c("U", "1dpi", "3dpi", "7dpi", "28dpi"))

# only use AP1 TFs that change GEX upon injury
ap1_injury <- c("Fosl2", "Fosb", "Fos", "Fosl1", 
                "Jun", "Junb", 
                "Atf3", "Atf5", "Atf4", "Atf1", "Atf6b", 
                "Jdp2", 
                "Batf", 
                "Mafk", "Mafg", "Maff")
glia <- AddModuleScore(glia, features = list(ap1_injury), name = "ap1_injury")
DotPlot(glia, features = "ap1_injury1", cols = c("lightgrey", "black"))

ap1_plots <- list()
for(i in ap1_injury){
  ap1_plots[[i]] <- ggplot(df_ap1, aes_string(x="timepoint", y=i, group="celltype")) +
    geom_line(aes(color=celltype))+
    geom_point(aes(color=celltype))+
    scale_color_manual(values=c(palette[3], palette[4], palette[7], palette[6], palette[5]))+
    theme_classic()
  
}
ggarrange(plotlist = ap1_plots, common.legend = T, nrow = 4, ncol = 4)

