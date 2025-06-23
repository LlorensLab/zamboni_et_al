# cut and tag analysis
setwd("~/Dropbox/Postdoc/sequencing/cutandtag")
sample_sheet <- read.table("Sample_sheet_run_250206.csv", sep = ";",skip = 8, header = T)

## Summarize the duplication information from the picard summary outputs.
sampleList = sample_sheet$I7_Index_ID[1:8]
sampleInfo_all = sample_sheet$Sample_Name[1:8]
names(sampleInfo_all) = sampleList
dupResult = list()
for(sample in sampleList){
  dupRes = read.table(paste0("picard_summary/", sample, "_picard.rmDup.txt"), header = TRUE, fill = TRUE)
  
  sampleInfo = strsplit(sampleInfo_all[[sample]], "_")[[1]]
  dupResult[[sample]] = data.frame(Ab = sampleInfo[1], 
                         Condition = sampleInfo[2], 
                         Replicate = sampleInfo[3], 
                         MappedFragNum = dupRes$READ_PAIRS_EXAMINED[1] %>% 
                           as.character %>% 
                           as.numeric, 
                         DuplicationRate = gsub(",", ".", dupRes$PERCENT_DUPLICATION[1]) %>% 
                           as.character %>% 
                           as.numeric * 100, 
                         EstimatedLibrarySize = dupRes$ESTIMATED_LIBRARY_SIZE[1] %>% 
                           as.character %>% 
                           as.numeric) #%>% 
    #mutate(UniqueFragNum = MappedFragNum * (1-DuplicationRate/100))  %>% rbind(dupResult, .)
}
dupResult_df <- do.call(rbind, dupResult)

#The estimated library sizes is proportional to the abundance of the targeted epitope and to the quality of the antibody used, while the estimated library sizes of IgG samples are expected to be very low.

# Get fragment size
##=== R command ===## 
## Collect the fragment size information
fragLen = c()
for(sample in sampleList){
  sampleInfo = strsplit(sampleInfo_all[[sample]], "_")[[1]]
  fragLen = read.table(paste0("fragmentLen/", sample, "_fragmentLen.txt"), header = FALSE) %>% 
    mutate(fragLen = V1 %>% as.numeric, 
           fragCount = V2 %>% as.numeric, 
           Weight = as.numeric(V2)/sum(as.numeric(V2)), 
           Ab = sampleInfo[1], 
           Condition = sampleInfo[2], 
           Replicate = sampleInfo[3], 
           sampleInfo = sample) %>% rbind(fragLen, .) 
}
fragLen$sampleInfo = factor(fragLen$sampleInfo, levels = sampleList)
fragLen$Ab_condition = paste0(fragLen$Ab, "_", fragLen$Condition)

## Generate the fragment size density plot (violin plot)
p1 <- fragLen %>% ggplot(aes(x = sampleInfo, y = fragLen, weight = Weight, fill = Ab_condition)) +
  geom_violin(bw = 5) +
  scale_y_continuous(breaks = seq(0, 800, 50)) +
  scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.9, option = "magma", alpha = 0.8) +
  scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
  theme_bw(base_size = 20) +
  ggpubr::rotate_x_text(angle = 20) +
  ylab("Fragment Length") +
  xlab("")

p2 <- fragLen %>% ggplot(aes(x = fragLen, y = fragCount, color = sampleInfo, group = sampleInfo, linetype = Replicate)) +
  geom_line(size = 1) +
  scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9, option = "magma") +
  theme_bw(base_size = 20) +
  xlab("Fragment Length") +
  ylab("Count") +
  coord_cartesian(xlim = c(0, 500))

ggarrange(p1, p2, ncol = 2, nrow=1, common.legend = TRUE, legend="bottom")

# with duplicates
######
fragLen = c()
for(sample in sampleList){
  sampleInfo = strsplit(sampleInfo_all[[sample]], "_")[[1]]
  fragLen = read.table(paste0("fragmentLen/", sample, "_wDups_fragmentLen.txt"), header = FALSE) %>% 
    mutate(fragLen = V1 %>% as.numeric, 
           fragCount = V2 %>% as.numeric, 
           Weight = as.numeric(V2)/sum(as.numeric(V2)), 
           Ab = sampleInfo[1], 
           Condition = sampleInfo[2], 
           Replicate = sampleInfo[3], 
           sampleInfo = sample) %>% rbind(fragLen, .) 
}
fragLen$sampleInfo = factor(fragLen$sampleInfo, levels = sampleList)
fragLen$Ab_condition = paste0(fragLen$Ab, "_", fragLen$Condition)

## Generate the fragment size density plot (violin plot)
fragLen %>% ggplot(aes(x = sampleInfo, y = fragLen, weight = Weight, fill = Ab_condition)) +
  geom_violin(bw = 5) +
  scale_y_continuous(breaks = seq(0, 800, 50)) +
  scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.9, option = "magma", alpha = 0.8) +
  scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
  theme_bw(base_size = 20) +
  ggpubr::rotate_x_text(angle = 20) +
  ylab("Fragment Length") +
  xlab("")

fragLen %>% ggplot(aes(x = fragLen, y = fragCount, color = Ab_condition, group = sampleInfo, linetype = Replicate)) +
  geom_line(size = 1) +
  scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9, option = "magma") +
  theme_bw(base_size = 20) +
  xlab("Fragment Length") +
  ylab("Count") +
  coord_cartesian(xlim = c(0, 500))
######

peakN = c()
peakWidth = c()
peakType = c("stringent", "relaxed")
for(sample in c("sox9_1dpi_f", "sox9_1dpi_m", "sox9_u_f", "sox9_u_m")){
    for(type in peakType){
      sampleInfo = strsplit(sample, "_")[[1]]
      peakInfo = read.table(paste0("peakCalling/", sample, "_seacr_", type, ".peaks.", type, ".bed"), header = FALSE, fill = TRUE)  %>% mutate(width = abs(V3-V2))
      peakN = data.frame(peakN = nrow(peakInfo), peakType = type, Condition = sampleInfo[2], Replicate = sampleInfo[3]) %>% rbind(peakN, .)
      peakWidth = data.frame(width = peakInfo$width, peakType = type, Condition = sampleInfo[2], Replicate = sampleInfo[3])  %>% rbind(peakWidth, .)
  }
}

condL = c("1dpi", "u")
repL = c("f", "m")
peakType = c("stringent", "relaxed")
peakOverlap = c()
for(type in peakType){
  for(cond in condL){
    overlap.gr = GRanges()
    for(rep in repL){
      peakInfo = read.table(paste0("peakCalling/sox9_", cond, "_", rep, "_seacr_", type, ".peaks.", type, ".bed"), header = FALSE, fill = TRUE)  %>% mutate(width = abs(V3-V2))
      peakInfo.gr = GRanges(peakInfo$V1, IRanges(start = peakInfo$V2, end = peakInfo$V3), strand = "*")
      if(length(overlap.gr) >0){
        overlap.gr = overlap.gr[findOverlaps(overlap.gr, peakInfo.gr)@from]
      }else{
        overlap.gr = peakInfo.gr
        
      }
    }
    peakOverlap = data.frame(peakReprod = length(overlap.gr), Condition = cond, peakType = type) %>% rbind(peakOverlap, .)
  }
}

peakReprod = left_join(peakN, peakOverlap, by = c("Condition", "peakType")) %>% mutate(peakReprodRate = peakReprod/peakN * 100)
peakReprod %>% dplyr::select(Condition, Replicate, peakType, peakN, peakReprodNum = peakReprod, peakReprodRate)

inPeakData = c()
## overlap with bam file to get count
samples = c("sox9_u_m", "sox9_u_f", "sox9_1dpi_m", "sox9_1dpi_f")
samples_names = c("Ad2.10", "Ad2.11", "Ad2.13", "Ad2.14")
for(sample in 1:4){
    peakRes = read.table(paste0("peakCalling/", samples[sample], "_seacr_relaxed.peaks.relaxed.bed"), header = FALSE, fill = TRUE)
    peak.gr = GRanges(seqnames = peakRes$V1, IRanges(start = peakRes$V2, end = peakRes$V3), strand = "*")
    bamFile = paste0("bam/", samples_names[sample], ".sorted.rmDup.mapped.bam")
    fragment_counts <- getCounts(bamFile, peak.gr, paired = TRUE, by_rg = FALSE, format = "bam")
    inPeakN = counts(fragment_counts)[,1] %>% sum
    sampleInfo = strsplit(samples[sample], "_")[[1]]
    inPeakData = rbind(inPeakData, data.frame(inPeakN = inPeakN, Condition = sampleInfo[2], Replicate = sampleInfo[3]))
  }
inPeakData$MappedFrag = dupResult_df$EstimatedLibrarySize[5:8]
inPeakData$FriPs = inPeakData$inPeakN/inPeakData$MappedFrag * 100
inPeakData

p1 <- peakN %>% ggplot(aes(x = Condition, y = peakN, fill = Condition)) +
  geom_boxplot() +
  geom_jitter(aes(color = Replicate), position = position_jitter(0.15)) +
  facet_grid(~peakType) +
  scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.55, option = "magma", alpha = 0.8) +
  scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
  theme_bw(base_size = 18) +
  ylab("Number of Peaks") +
  xlab("")

p2 <- peakWidth %>% ggplot(aes(x = Condition, y = width, fill = Condition)) +
  geom_violin() +
  facet_grid(Replicate~peakType) +
  scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.55, option = "magma", alpha = 0.8) +
  scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
  scale_y_continuous(trans = "log", breaks = c(400, 3000, 22000)) +
  theme_bw(base_size = 18) +
  ylab("Width of Peaks") +
  xlab("")

p3 <- peakReprod %>% ggplot(aes(x = Condition, y = peakReprodRate, fill = Condition, label = round(peakReprodRate, 2))) +
  geom_bar(stat = "identity") +
  geom_text(vjust = 0.1) +
  facet_grid(Replicate~peakType) +
  scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.55, option = "magma", alpha = 0.8) +
  scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
  theme_bw(base_size = 18) +
  ylab("% of Peaks Reproduced") +
  xlab("")

p4 <- inPeakData %>% ggplot(aes(x = Condition, y = FriPs, fill = Condition, label = round(FriPs, 2))) +
  geom_boxplot() +
  geom_jitter(aes(color = Replicate), position = position_jitter(0.15)) +
  scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.55, option = "magma", alpha = 0.8) +
  scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
  theme_bw(base_size = 18) +
  ylab("% of Fragments in Peaks") +
  xlab("")

ggarrange(p1, p2, p4, ncol = 3, nrow=1, common.legend = TRUE, legend="bottom")

bigwig <- list(
  sox9_u_m = "bigwigs/Ad2.10.sorted.rmDup.fragments.bw", 
  sox9_u_f = "bigwigs/Ad2.11.sorted.rmDup.fragments.bw",
  sox9_1dpi_m = "bigwigs/Ad2.13.sorted.rmDup.fragments.bw",
  sox9_1dpi_f = "bigwigs/Ad2.14.sorted.rmDup.fragments.bw")

BigwigTrack(
  region = "chr2-102872975-102873741",
  bigwig = bigwig,
  extend.upstream = 2000,
  extend.downstream = 2000,
  type = "coverage",
  bigwig.scale = "common")

BigwigTrack(
  region = "chr2-102872806-102874142",
  bigwig = bigwig,
  extend.upstream = 2000,
  extend.downstream = 2000,
  type = "coverage",
  bigwig.scale = "common")

peaks_sox9_1dpi_f_gr_stringent = get_granges("peakCalling/sox9_1dpi_f_seacr_stringent.peaks.stringent.bed")
peaks_sox9_1dpi_m_gr_stringent = get_granges("peakCalling/sox9_1dpi_m_seacr_stringent.peaks.stringent.bed")
peaks_sox9_u_f_gr_stringent = get_granges("peakCalling/sox9_u_f_seacr_stringent.peaks.stringent.bed")
peaks_sox9_u_m_gr_stringent = get_granges("peakCalling/sox9_u_m_seacr_stringent.peaks.stringent.bed")

subsetByOverlaps(peaks_sox9_1dpi_f_gr_stringent, astro_i_gr, type = "any", ignore.strand = T)
subsetByOverlaps(peaks_sox9_1dpi_m_gr_stringent, astro_i_gr, type = "any", ignore.strand = T)
subsetByOverlaps(peaks_sox9_u_f_gr_stringent, astro_i_gr, type = "any", ignore.strand = T)
subsetByOverlaps(peaks_sox9_u_m_gr_stringent, astro_i_gr, type = "any", ignore.strand = T)
#
#44k overlap (29k stringent)

#combined_sox9_1dpi_peaks <- reduce(x = c(peaks_sox9_1dpi_f_gr, peaks_sox9_1dpi_m_gr))
#combined_sox9_1dpi_peaks <- combined_sox9_1dpi_peaks[combined_sox9_1dpi_peaks@ranges@width < 2000]
#overlaps <- subsetByOverlaps(combined_sox9_1dpi_peaks, astro_i_sox9_gr, type = "any", ignore.strand = T)

subsetByOverlaps(peaks_sox9_1dpi_f_gr, astro_i_sox9_gr, type = "any", ignore.strand = T)
subsetByOverlaps(peaks_sox9_1dpi_m_gr, astro_i_sox9_gr, type = "any", ignore.strand = T)
subsetByOverlaps(peaks_sox9_1dpi_f_gr_stringent, astro_i_sox9_gr, type = "any", ignore.strand = T)
subsetByOverlaps(peaks_sox9_1dpi_m_gr_stringent, astro_i_sox9_gr, type = "any", ignore.strand = T)
subsetByOverlaps(peaks_sox9_u_f_gr_stringent, astro_i_sox9_gr, type = "any", ignore.strand = T)
subsetByOverlaps(peaks_sox9_u_m_gr_stringent, astro_i_sox9_gr, type = "any", ignore.strand = T)

BigwigTrack(
  region = "chr9-101445605-101447147",
  bigwig = bigwig,
  extend.upstream = 1000,
  extend.downstream = 1000,
  type = "coverage",
  bigwig.scale = "common")

patterns_ap1 <- paste0("pattern_", c(3, 13, 22, 26, 40))
astro_i_interpreted_ap1 <- astro_i_interpreted %>%
  filter(Pattern %in% patterns_ap1)
astro_i_interpreted_sox9_ap1 <- intersect(astro_i_interpreted_sox9$Example_idx, astro_i_interpreted_ap1$Example_idx)
astro_i_sox9_ap1_peaks <- astro_i_bed[astro_i_interpreted_sox9_ap1,]
astro_i_sox9_ap1_gr <- GRanges(seqnames = astro_i_sox9_ap1_peaks$V1, 
                           ranges = paste0(astro_i_sox9_ap1_peaks$V2, "-", astro_i_sox9_ap1_peaks$V3)
)

subsetByOverlaps(peaks_sox9_1dpi_f_gr, astro_i_sox9_ap1_gr, type = "any", ignore.strand = T)
subsetByOverlaps(peaks_sox9_1dpi_m_gr, astro_i_sox9_ap1_gr, type = "any", ignore.strand = T)
subsetByOverlaps(peaks_sox9_1dpi_f_gr_stringent, astro_i_sox9_ap1_gr, type = "any", ignore.strand = T)
subsetByOverlaps(peaks_sox9_1dpi_m_gr_stringent, astro_i_sox9_ap1_gr, type = "any", ignore.strand = T)
subsetByOverlaps(peaks_sox9_u_f_gr_stringent, astro_i_sox9_ap1_gr, type = "any", ignore.strand = T)
subsetByOverlaps(peaks_sox9_u_m_gr_stringent, astro_i_sox9_ap1_gr, type = "any", ignore.strand = T)

# 2676 regions with Sox9, of which 1009 (M) and 1015 (F) found overlap in the CUT&Tag
# 197 regions with Sox9/Ap1, of which 89 (M) and 75 (F) found overlap in the CUT&Tag
# total CUT&Tag: 69200 (M, 1dpi), 904521 (F, 1 dpi), 38471 (M, U) and 23901 (F, U)
