# Preprocessing with ArchR for peak calling

source("../source.R")
path <- "data/"
dirs <- list.dirs(path, recursive = F, full.names = F)

addArchRThreads(threads = 5) 
addArchRGenome("mm10")

inputFiles <- paste0(path, dirs, "/outs/atac_fragments.tsv.gz")
names(inputFiles) <- dirs

# load metadata
barcodeList <- lapply(seq_along(dirs), function(x){
  csvFiles <- paste0(path, dirs[x], "/outs/per_barcode_metrics.csv")
  if(file.exists(csvFiles)){
    meta <- read.csv(csvFiles, row.names = 1) %>%
      dplyr::filter(is_cell == 1) 
    return(rownames(meta))
    
  } else {
    csvFiles <- paste0(path, dirs[x], "/outs/singlecell.csv")
    meta <- read.csv(csvFiles, row.names = 1) %>%
      dplyr::filter(is__cell_barcode == 1)
    return(rownames(meta))
  }
})

names(barcodeList) <- dirs


# create arrow files
ArrowFiles <- createArrowFiles(
  inputFiles = inputFiles,
  sampleNames = names(inputFiles),
  validBarcodes = barcodeList,
  geneAnnotation = getGeneAnnotation(),
  genomeAnnotation = getGenomeAnnotation(),
  QCDir = "QualityControl",
  logFile = createLogFile("createArrows"),
  minTSS = 1, #Dont set this too high because you can always increase later
  minFrags = 300, 
  addTileMat = T,
  addGeneScoreMat = T,
  force = T
)

ArrowFiles <- paste0(dirs, ".arrow")
archrproj <- ArchRProject(
  ArrowFiles = ArrowFiles, 
  outputDirectory = "proj_archr",
  copyArrows = TRUE #This is recommended so that if you modify the Arrow files you have an original copy for later usage.
)

# dimensionality reduction
archrproj <- addIterativeLSI(
  ArchRProj = archrproj,
  useMatrix = "TileMatrix", 
  name = "IterativeLSI", 
  iterations = 2, 
  clusterParams = list( #See Seurat::FindClusters
    resolution = c(0.2), 
    sampleCells = 10000, 
    n.start = 10
  ), 
  varFeatures = 25000, 
  dimsToUse = 1:30
)

archrproj <- addHarmony(
  ArchRProj = archrproj,
  reducedDims = "IterativeLSI",
  name = "Harmony",
  groupBy = "Sample"
)

archrproj <- addClusters(
  input = archrproj,
  reducedDims = "IterativeLSI",
  method = "Seurat",
  name = "Clusters",
  resolution = 1
)

archrproj <- addUMAP(
  ArchRProj = archrproj, 
  reducedDims = "IterativeLSI", 
  name = "UMAP", 
  nNeighbors = 30, 
  minDist = 0.5, 
  metric = "cosine"
)

archrproj <- addGroupCoverages(ArchRProj = archrproj, groupBy = "Clusters")

# peak calling with MACS2
pathToMacs2 <- findMacs2()
archrproj <- addReproduciblePeakSet(
  ArchRProj = archrproj, 
  groupBy = "Clusters", 
  pathToMacs2 = pathToMacs2
)

archr_peaks <- getPeakSet(archrproj)

#create feature matrices
mtx_cluster <- list()
for(i in seq_along(dirs)){
  frag <- CreateFragmentObject(inputFiles[i], 
                               cells = barcodeList[[i]])
  # quantify counts in each peak
  mtx_cluster[[i]] <- FeatureMatrix(
    fragments = frag,
    features = archr_peaks,
    process_n = 20000
  )
}


# save objects
saveArchRProject(ArchRProj = archrproj, outputDirectory = "archrproj_multiome", load = FALSE)
saveRDS(archr_peaks, "archr_peaks_multiome.rds")
saveRDS(mtx_cluster, "featurematrix_archr_clusterpeaks_multiome.rds")

