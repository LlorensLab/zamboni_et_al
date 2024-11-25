##code for:
#differential analysis at the cell type level
#identifying peaks at the TSS vs distal peaks
#linking peaks (performed for each cell type)
#differential analysis in injury

#differential analysis
source("../source.R")
multiome <- readRDS("merged_multiome.rds")

clusters <- Idents(multiome) %>% unique()

# enriched genes and regions per cell type 
# run with wsnn_res_1.1 for the cell subtypes in fig1
DefaultAssay(multiome) <- "RNA"
for(i in clusters){
  markers_cl <- FindMarkers(multiome, 
                            ident.1 = i, 
                            only.pos = TRUE,
                            min.pct = 0.05
  )
  
  markers_cl$gene <- rownames(markers_cl) # add column "gene"
  markers_cl$comparison <- i
  markers_rna[[i]] <- markers_cl
}
names(markers_rna) <- clusters

DefaultAssay(multiome) <- "peaks"
for(i in clusters){
  markers_cl <- FindMarkers(multiome, 
                            ident.1 = i, 
                            only.pos = TRUE,
                            min.pct = 0.05
  )
  
  markers_cl$peak <- rownames(markers_cl)
  markers_cl$comparison <- i
  markers_peaks[[i]] <- markers_cl
}
names(markers_peaks) <- clusters


for(i in clusters){
  markers_peaks[[i]] <- markers_peaks[[i]] %>%
    dplyr::filter(p_val_adj < 0.05)
  
  markers_rna[[i]] <- markers_rna[[i]] %>%
    dplyr::filter(p_val_adj < 0.05)
}

#tss vs distal peaks
DefaultAssay(multiome) <- "peaks"

tss.positions <- GetTSSPositions(ranges = annotations)
tss.positions <- Extend(
  x = tss.positions,
  upstream = 1000,
  downstream = 100,
  from.midpoint = TRUE
) # obtain region around TSS

closest_tss <- ClosestFeature(multiome, rownames(multiome), annotation = tss.positions)
tss_peaks <- closest_tss %>% filter(distance == 0) # peaks falling within promoter region (between -1000 and +100 from TSS)
distal_peaks <- closest_tss %>% filter(distance > 0) # peaks falling outside promoter region - assigned as distal

#save 
saveRDS(markers_rna, "markers_rna_multiome_celltypes.rds")
saveRDS(markers_peaks, "markers_peaks_multiome_celltypes.rds")

# link peaks
linkedpeaks <- list()
multiome <- RegionStats(multiome, BSgenome.Acahirinus.ENSEMBL.1.softmasked)
for(i in levels(multiome)){
  celltype <- subset(multiome, idents = i)
  celltype <- LinkPeaks(celltype, 
                        peak.assay = "peaks", 
                        expression.assay = "RNA", 
                        distance = 500000)
  linkedpeaks[[i]] <- celltype@assays[["peaks"]]@links
  saveRDS(linkedpeaks[[i]], paste0(i, "_linkpeaks.rds"))
}


# differential analysis injury
multiome$condition <- factor(multiome$condition,
                             levels = c("U", "1dpi", "3dpi", "7dpi", "28dpi"))

markers_Uvsothers_peaks <- list()
markers_Uvsothers_peaks_all <- list()
for(i in clusters){
  celltype <- subset(multiome, ident = i)
  Idents(celltype) <- "condition"
  DefaultAssay(celltype) <- "peaks"
  
  for(j in c("1dpi", "3dpi", "7dpi", "28dpi")){
    markers_Uvsothers_peaks[[j]] <- FindMarkers(celltype, 
                                                ident.1 = j,
                                                ident.2 = "U",
                                                only.pos = FALSE,
                                                min.pct = 0.05)
    markers_Uvsothers_peaks[[j]]$comparison <- paste0(i, "_", j)
    markers_Uvsothers_peaks[[j]]$peak <- rownames(markers_Uvsothers_peaks[[j]])
    
  }
  markers_Uvsothers_peaks_all[[i]] <- do.call("rbind", markers_Uvsothers_peaks)
  rm(celltype)
}

markers_Uvsothers_rna <- list()
markers_Uvsothers_rna_all <- list()
for(i in clusters){
  celltype <- subset(multiome, ident = i)
  Idents(celltype) <- "condition"
  DefaultAssay(celltype) <- "RNA"
  
  for(j in c("1dpi", "3dpi", "7dpi", "28dpi")){
    markers_Uvsothers_rna[[j]] <- FindMarkers(celltype, 
                                              ident.1 = j,
                                              ident.2 = "U",
                                              only.pos = FALSE,
                                              min.pct = 0.05, 
                                              logfc.threshold = 0.32)
    markers_Uvsothers_rna[[j]]$comparison <- paste0(i, "_", j)
    markers_Uvsothers_rna[[j]]$gene <- rownames(markers_Uvsothers_rna[[j]])
    
  }
  markers_Uvsothers_rna_all[[i]] <- do.call("rbind", markers_Uvsothers_rna)
  rm(celltype)
}

for(i in clusters){
  markers_Uvsothers_peaks_all[[i]] <- markers_Uvsothers_peaks_all[[i]] %>%
    dplyr::filter(p_val < 0.05 & abs(avg_log2FC) > 0.5)
  
  markers_Uvsothers_rna_all[[i]] <- markers_Uvsothers_rna_all[[i]] %>%
    dplyr::filter(p_val < 0.05 & abs(avg_log2FC) > 0.5)
}


#save
saveRDS(markers_Uvsothers_rna_all, "markers_Uvsothers_rna_all.rds")
saveRDS(markers_Uvsothers_peaks_all, "markers_Uvsothers_peaks_all.rds")


