#clustering of DEGs and DARs
# injury-induced genes and peaks are clustered using TCseq into modules with the same dynamics

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
