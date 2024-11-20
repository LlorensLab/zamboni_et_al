## install chrombpnet and dependencies
conda create -n chrombpnet python=3.8
conda activate chrombpnet
conda install -y -c conda-forge -c bioconda samtools bedtools ucsc-bedgraphtobigwig pybigwig meme
pip install chrombpnet

## run bias model
# sort and index the bgzipped file
sort -k 1,1V -k 2,2n 1_input_bias.bed > sorted_fragments_bias.tsv
bgzip -@ 20 sorted_fragments_bias.tsv
tabix -f -p bed sorted_fragments_bias.tsv.gz

macs2 callpeak -f AUTO -t sorted_fragments_bias.tsv.gz -g 1.87e9 -B -p 0.01 --nomodel --extsize 200 --outdir macs2_from_frag -n bias_fragments

# remove blacklist
bedtools slop -i mm10-blacklist.v2.bed.gz -g mm10.chrom.sizes -b 1057 > temp.bed
bedtools intersect -v -a bias_fragments.narrowPeak -b temp.bed  > bias_peaks_no_blacklist.bed

# make fold 0 splits
mkdir splits
chrombpnet prep splits -c mm10.chrom.sizes -tcr chr1 chr3 chr6 -vcr chr8 chr19 -op splits/fold_0

# prepare non-peak regions
chrombpnet prep nonpeaks -g genome.fa -p bias_peaks_no_blacklist.bed -c mm10.chrom.sizes -fl splits/fold_0.json -br mm10-blacklist.v2.bed.gz -o output_nonpeaks

# train bias model
chrombpnet bias pipeline \
        -ifrag sorted_fragments_bias.tsv.gz \
        -d "ATAC" \
        -g genome.fa \
        -c mm10.chrom.sizes \
        -p bias_peaks_no_blacklist.bed \
        -n output_nonpeaks_negatives.bed \
        -fl splits/fold_0.json \
        -b 0.5 \
        -o bias/ \
        -fp bias_model
        

## train cell type models (same for cell type-, injury-specific models) 
# only keep regions mapping to "chr"
for i in Astrocytes Ependymal Microglia Oligodendrocytes Macrophages Neurons_V Perivascular Endothelial OPCs Neurons_D
do
  grep '^chr' ML/celltype_fragments/${i}.bed > celltype_fragments_multiome_only/${i}_chr.bed
  macs2 callpeak -f AUTO \
  -t celltype_fragments_multiome_only/${i}_chr.bed \
  -g 1.87e9 -B -p 0.01 --nomodel --extsize 200 \
  --outdir macs2_celltype_peaks -n ${i}_fragments

done

mkdir filtered_peaks
mkdir negatives
bedtools slop -i mm10-blacklist.v2.bed.gz -g mm10.chr.sizes -b 1057 > temp.bed

# create peak sets without blacklist regions and prepare non-peak regions
for i in Astrocytes Ependymal Microglia Oligodendrocytes Macrophages Neurons_V Perivascular Endothelial OPCs Neurons_D
do
  bedtools intersect -v -a macs2_celltype_peaks/${i}_fragments_peaks.narrowPeak -b temp.bed  > filtered_peaks/${i}_peaks_no_blacklist.bed
  grep '^chr' filtered_peaks/${i}_peaks_no_blacklist.bed > filtered_peaks/${i}_peaks_no_blacklist_chr.bed
  chrombpnet prep nonpeaks -g genome.fa -p filtered_peaks/${i}_peaks_no_blacklist_chr.bed -c mm10.chr.sizes -fl splits/fold_0.json -br mm10-blacklist.v2.bed.gz -o negatives/${i}_nonpeaks
done

for i in Astrocytes Ependymal Microglia Oligodendrocytes Macrophages Neurons_V Perivascular Endothelial OPCs Neurons_D
do
  chrombpnet pipeline \
        -ifrag ML/celltype_fragments/${i}_chr.bed \
        -d "ATAC" \
        -g genome.fa \
        -c mm10.chr.sizes \
        -p filtered_peaks/${i}_peaks_no_blacklist.bed \
        -n negatives/A${i}_nonpeaks_negatives.bed \
        -fl splits/fold_0.json \
        -b small_frag_bias.h5 \
        -o models/${i} > logs/${i}_logs.txt 2>&1
done


# predict other sets of peaks (e.g., DARs enriched in specific cell types or in injury)
# obtain count and profile predictions
# compute base pair-level contribution scores
# identify important motifs and link to TFs
for i in Astrocytes Ependymal Microglia Oligodendrocytes Macrophages Neurons_V Perivascular Endothelial OPCs Neurons_D
do
  
  grep '^chr' dars_intersect_peaks/${i}_intersect_dars.bed > dars_intersect_peaks/${i}_intersect_dars_chr.bed
  
  chrombpnet pred_bw -bm models/${i}_fl0/models/bias_model_scaled.h5 \
    -cm models/${i}_fl0/models/chrombpnet.h5 \
    -r dars_intersect_peaks/${i}_intersect_dars.bed \
    -c mm10.chr.sizes -g genome.fa \
    -op interpret/{i}/{i}_dars_bigwig
    
  chrombpnet contribs_bw -m models/${i}_fl0/models/chrombpnet.h5 \
    -r dars_intersect_peaks/${i}_intersect_dars.bed \
    -c mm10.chrom.sizes -g genome.fa \
    -op interpret/{i}/{i}_dars_contrib
 
  modisco motifs -i interpret/{i}/{i}_dars_contrib.counts_scores.h5 -n 1000000 -o interpret/{i}/{i}_modisco.h5 -v
  modisco report -i interpret/{i}/{i}_modisco.h5 -o interpret/{i}/modisco_report/ -s interpret/{i}/modisco_report/ -m motifs.meme.txt

done

