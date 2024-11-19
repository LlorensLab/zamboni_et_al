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


