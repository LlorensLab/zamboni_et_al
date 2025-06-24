#!/bin/sh

#  cut_and_tag.sh
#  
#
#  Created by Margherita Zamboni on 2025-02-10.
#
# https://yezhengstat.github.io/CUTTag_tutorial

bowtie2-build genome.fa mus_genome

# trimming
zcat FT100037069/L01/FT100037069_L01_read_1.fq.gz | awk '{OFS=""}{if(NR%2 == 0){print substr($0,1,50)} else print;}'  | pigz -p20 > FT100037069_L01_r1_50bp.fastq.gz
zcat FT100037069/L01/FT100037069_L01_read_2.fq.gz | awk '{OFS=""}{if(NR%2 == 0){print substr($0,1,50)} else print;}'  | pigz -p20 > FT100037069_L01_r2_50bp.fastq.gz
zcat FT100037069/L01/FT100037069_L01_read_2.fq.gz | awk '{OFS=""}{if(NR%2 == 0){print substr($0,71,8)} else print;}'  | pigz -p20 > FT100037069_L01_i1.fastq.gz

export DEML=/mnt/storage/marzam/uppstore2018207/bioinfo_tools/deML_DNB/src
$DEML/deML --index samplesheet.txt \
-f FT100037069_L01_r1_50bp.fastq.gz -r FT100037069_L01_r2_50bp.fastq.gz \
-if1 FT100037069_L01_i1.fastq.gz \
--summary demux_stats.txt \
--outfile demux_50bp

mkdir aligned_50
for i in Ad2.2 Ad2.3 Ad2.7 Ad2.9 Ad2.10 Ad2.11 Ad2.13 Ad2.14;
do bowtie2 --local --very-sensitive --no-mixed --no-discordant --phred33 -I 10 -X 700 -x mus_genome/mus_genome -1 demux_50bp_${i}_r1.fq.gz -2 demux_50bp_${i}_r2.fq.gz -S aligned_50/${i}.sam;
done

##== linux command ==##
picardCMD="java -jar picard.jar"
mkdir picard_summary
mkdir aligned_sorted

for i in Ad2.2 Ad2.3 Ad2.7 Ad2.9 Ad2.10 Ad2.11 Ad2.13 Ad2.14;
do
## Sort by coordinate
java -jar ../picard.jar SortSam -I ${i}.sam -O ../aligned_sorted/${i}.sorted.sam --SORT_ORDER coordinate
done

cd ../
## mark duplicates
for i in Ad2.2 Ad2.3 Ad2.7 Ad2.9 Ad2.10 Ad2.11 Ad2.13 Ad2.14;
do
samtools addreplacerg -r "@RG\tID:RG1\tSM:SampleName\tPL:Illumina\tLB:Library.fa" -o aligned_sorted/${i}.sorted_added.sam aligned_sorted/${i}.sorted.sam

java -jar picard.jar MarkDuplicates -I aligned_sorted/${i}.sorted_added.sam -O aligned_sorted/${i}.sorted.dupMarked.sam --READ_NAME_REGEX null --METRICS_FILE picard_summary/${i}_picard.dupMark.txt
done
java -jar ../picard.jar MarkDuplicates -I Ad2.2.sorted.sam -O Ad2.2.sorted.dupMarked.sam --READ_NAME_REGEX null --METRICS_FILE Ad2.2_picard.dupMark.txt

## remove duplicates
for i in Ad2.2 Ad2.3 Ad2.7 Ad2.9 Ad2.10 Ad2.11 Ad2.13 Ad2.14;
do
java -jar picard.jar MarkDuplicates -I aligned_sorted/${i}.sorted_added.sam -O aligned_sorted/${i}.sorted.rmDup.sam --READ_NAME_REGEX null --REMOVE_DUPLICATES true --METRICS_FILE picard_summary/${i}_picard.rmDup.txt
done

for i in Ad2.2 Ad2.3 Ad2.7 Ad2.9 Ad2.10 Ad2.11 Ad2.13 Ad2.14;
do
## Sort by coordinate
java -jar picard.jar SortSam -I aligned_sorted/${i}.sorted.rmDup.sam -O aligned_sorted/${i}.sorted.rmDup.sorted.sam --SORT_ORDER queryname
done

# Get fragment size
##== linux command ==##
mkdir -p fragmentLen

## Extract the 9th column from the alignment sam file which is the fragment length
for i in Ad2.2 Ad2.3 Ad2.7 Ad2.9 Ad2.10 Ad2.11 Ad2.13 Ad2.14;
do
samtools view -F 0x04 aligned_sorted/${i}.sorted.rmDup.sam | awk -F'\t' 'function abs(x){return ((x < 0.0) ? -x : x)} {print abs($9)}' | sort | uniq -c | awk -v OFS="\t" '{print $2, $1/2}' >fragmentLen/${i}_fragmentLen.txt
done

##== linux command ==##
## Filter and keep the mapped read pairs
mkdir bam
mkdir bed
for i in Ad2.2 Ad2.3 Ad2.7 Ad2.9 Ad2.10 Ad2.11 Ad2.13 Ad2.14;
do
samtools view -bS -F 0x04 aligned_sorted/${i}.sorted.rmDup.sorted.sam >bam/${i}.sorted.rmDup.mapped.bam
bedtools bamtobed -i bam/${i}.sorted.rmDup.mapped.bam -bedpe >bed/${i}.sorted.rmDup.mapped.bed
awk '$1==$4 && $6-$2 < 1000 {print $0}' bed/${i}.sorted.rmDup.mapped.bed >bed/${i}.sorted.rmDup.clean.bed
cut -f 1,2,6 bed/${i}.sorted.rmDup.clean.bed | sort -k1,1 -k2,2n -k3,3n  >bed/${i}.sorted.rmDup.fragments.bed
done

# peak calling
# create bedgraphs
mkdir bedgraph
for i in Ad2.2 Ad2.3 Ad2.7 Ad2.9 Ad2.10 Ad2.11 Ad2.13 Ad2.14;
do
bedtools genomecov -bg -i bed/${i}.sorted.rmDup.fragments.bed -g mm10.chrom.sizes > bedgraph/${i}.sorted.rmDup.fragments.bedgraph
done

seacr="SEACR/SEACR_1.3.sh"
mkdir peakCalling

bash $seacr bedgraph/Ad2.10.sorted.rmDup.fragments.bedgraph bedgraph/Ad2.2.sorted.rmDup.fragments.bedgraph norm relaxed peakCalling/sox9_u_m_seacr_relaxed.peaks
bash $seacr bedgraph/Ad2.11.sorted.rmDup.fragments.bedgraph bedgraph/Ad2.3.sorted.rmDup.fragments.bedgraph norm relaxed peakCalling/sox9_u_f_seacr_relaxed.peaks
bash $seacr bedgraph/Ad2.13.sorted.rmDup.fragments.bedgraph bedgraph/Ad2.7.sorted.rmDup.fragments.bedgraph norm relaxed peakCalling/sox9_1dpi_m_seacr_relaxed.peaks
bash $seacr bedgraph/Ad2.14.sorted.rmDup.fragments.bedgraph bedgraph/Ad2.9.sorted.rmDup.fragments.bedgraph norm relaxed peakCalling/sox9_1dpi_f_seacr_relaxed.peaks

bash $seacr bedgraph/Ad2.10.sorted.rmDup.fragments.bedgraph bedgraph/Ad2.2.sorted.rmDup.fragments.bedgraph norm stringent peakCalling/sox9_u_m_seacr_stringent.peaks
bash $seacr bedgraph/Ad2.11.sorted.rmDup.fragments.bedgraph bedgraph/Ad2.3.sorted.rmDup.fragments.bedgraph norm stringent peakCalling/sox9_u_f_seacr_stringent.peaks
bash $seacr bedgraph/Ad2.13.sorted.rmDup.fragments.bedgraph bedgraph/Ad2.7.sorted.rmDup.fragments.bedgraph norm stringent peakCalling/sox9_1dpi_m_seacr_stringent.peaks
bash $seacr bedgraph/Ad2.14.sorted.rmDup.fragments.bedgraph bedgraph/Ad2.9.sorted.rmDup.fragments.bedgraph norm stringent peakCalling/sox9_1dpi_f_seacr_stringent.peaks

for i sox9_1dpi_f sox9_1dpi_m sox9_u_f sox9_u_m;
do
awk 'BEGIN {OFS="\t"} {print $1, $2, $3, "Peak_"NR, ".", ".", ".", ".", ".", int(($3-$2)/2)}' ${i}_seacr_stringent.peaks.stringent.bed > ${i}_seacr_stringent_10col.bed
done

awk 'BEGIN {OFS="\t"} {print $1, $2, $3, "Peak_"NR, ".", ".", ".", ".", ".", int(($3-$2)/2)}' sox9_1dpi_f_seacr_stringent.peaks.stringent.bed > sox9_1dpi_f_seacr_stringent_10col.bed
awk 'BEGIN {OFS="\t"} {print $1, $2, $3, "Peak_"NR, ".", ".", ".", ".", ".", int(($3-$2)/2)}' sox9_1dpi_m_seacr_stringent.peaks.stringent.bed > sox9_1dpi_m_seacr_stringent_10col.bed
awk 'BEGIN {OFS="\t"} {print $1, $2, $3, "Peak_"NR, ".", ".", ".", ".", ".", int(($3-$2)/2)}' sox9_u_f_seacr_stringent.peaks.stringent.bed > sox9_u_f_seacr_stringent_10col.bed
awk 'BEGIN {OFS="\t"} {print $1, $2, $3, "Peak_"NR, ".", ".", ".", ".", ".", int(($3-$2)/2)}' sox9_u_m_seacr_stringent.peaks.stringent.bed > sox9_u_m_seacr_stringent_10col.bed

mkdir bigwigs
for i in Ad2.2 Ad2.3 Ad2.7 Ad2.9 Ad2.10 Ad2.11 Ad2.13 Ad2.14;
do
bedGraphToBigWig bedgraph/${i}.sorted.rmDup.fragments.bedgraph mm10.chrom.sizes bigwigs/${i}.sorted.rmDup.fragments.bw
done

mkdir bigwig_samtools
for i in Ad2.2 Ad2.3 Ad2.7 Ad2.9 Ad2.10 Ad2.11 Ad2.13 Ad2.14;
do
samtools sort -o bam/${i}.mapped.bam bam/${i}.sorted.rmDup.mapped.bam
samtools index bam/${i}.mapped.bam
bamCoverage -b bam/${i}.mapped.bam -o bigwig_samtools/${i}_raw.bw
done

# index files
cd bedgraph
for i in Ad2.2 Ad2.3 Ad2.7 Ad2.9 Ad2.10 Ad2.11 Ad2.13 Ad2.14;
do
awk -v name="$i" 'BEGIN {OFS="\t"} {print $1, $2, $3, name, $4}' "${i}.sorted.rmDup.fragments.bedgraph" > "${i}.sorted.rmDup.fragments.tsv"

bgzip ${i}.sorted.rmDup.fragments.tsv
tabix -p bed ${i}.sorted.rmDup.fragments.tsv.gz
done

gunzip *.sorted.rmDup.fragments.tsv.gz
cat *.sorted.rmDup.fragments.tsv > cutandtag_merged_fragments.tsv
sort -k 1,1 -k2,2n cutandtag_merged_fragments.tsv > cutandtag_merged_fragments_sorted.tsv
bgzip cutandtag_merged_fragments_sorted.tsv
tabix -p bed cutandtag_merged_fragments_sorted.tsv.gz

## call peaks with macs2
mkdir macs2_peaks
cd bedgraph
for i in Ad2.2 Ad2.3 Ad2.7 Ad2.9 Ad2.10 Ad2.11 Ad2.13 Ad2.14;
do
macs2 callpeak -f AUTO -t ${i}.sorted.rmDup.fragments.tsv.gz -g 1.87e9 -B -p 0.01 --nomodel --extsize 200 --outdir ../macs2_peaks -n ${i};
done

cd ../macs2_peaks
for i in Ad2.2 Ad2.3 Ad2.7 Ad2.9 Ad2.10 Ad2.11 Ad2.13 Ad2.14;
do
awk '$10 >= 50 && $10 <= 2500' ${i}_peaks.narrowPeak > ${i}_peaks.filtered.bed;
done
