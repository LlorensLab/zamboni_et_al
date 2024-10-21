#!/bin/sh

#  mgi_demultiplex.sh
#  
#
#  Created by Margherita Zamboni on 2023-04-17.
#  

#[Sample Name]_S1_L00[Lane Number]_[Read Type]_001.fastq.gz

# GEX
# R1 28; R2 90; i7 10; i5 10 (edited)

# GEX raw files in V350145982
# read 1 doesn't need trimming
export DIR=/proj/uppstore2018207/Multiome_MGI_230414/fastq

cd $DIR/gex
zcat $DIR/V350145982/*/*read_1.fq.gz | pigz -p20 > V350145982_S1_L001_R1_001.fastq.gz &
zcat $DIR/V350145982/*/*read_2.fq.gz | awk '{OFS=""}{if(NR%2 == 0){print substr($0,1,90)} else print;}'  | pigz -p20 > V350145982_S1_L001_R2_001.fastq.gz &
zcat $DIR/V350145982/*/*read_2.fq.gz | awk '{OFS=""}{if(NR%2 == 0){print substr($0,91,10)} else print;}'  | pigz -p20 > V350145982_S1_L001_I1_001.fastq.gz &
zcat $DIR/V350145982/*/*read_2.fq.gz | awk '{OFS=""}{if(NR%2 == 0){print substr($0,101,10)} else print;}'  | pigz -p20 > V350145982_S1_L001_I2_001.fastq.gz

# ATAC
export DIR=/proj/uppstore2018207/Multiome_MGI_230414/fastq
export SEQTK=/proj/snic2022-22-388/bioinfo_tools/seqtk

cd $DIR/atac
# Read structure to trim the 4 ATAC samples to:
# R1 50; i7 8; i5 24; R2 50
# i5 should be revcom with seqtk (barcode)
zcat $DIR/V350145980/*/*read_1.fq.gz | pigz -p20 > V350145980_S1_L001_R1_001.fastq.gz &
zcat $DIR/V350145980/*/*read_2.fq.gz | awk '{OFS=""}{if(NR%2 == 0){print substr($0,1,50)} else print;}'  | pigz -p20 > V350145980_S1_L001_R2_001.fastq.gz &
zcat $DIR/V350145980/*/*read_2.fq.gz | awk '{OFS=""}{if(NR%2 == 0){print substr($0,51,8)} else print;}'  | pigz -p20 > V350145980_S1_L001_I1_001.fastq.gz &
zcat $DIR/V350145980/*/*read_2.fq.gz | awk '{OFS=""}{if(NR%2 == 0){print substr($0,59,24)} else print;}'  | $SEQTK/seqtk seq -r | pigz -p20 > V350145980_S1_L001_I2_001.fastq.gz


# demultiplex
module load gcc

export DIR=/proj/uppstore2018207/Multiome_MGI_230414
export DEML=/proj/snic2022-22-388/bioinfo_tools/deML_DNB/src

cd $DIR/fastq/atac
$DEML/deML --index $DIR/samplesheet_atac.txt \
-f V350145980_S1_L001_R1_001.fastq.gz -r V350145980_S1_L001_R2_001.fastq.gz \
-if1 V350145980_S1_L001_I1_001.fastq.gz -if2 V350145980_S1_L001_I2_001.fastq.gz \
--summary demux_stats.txt \
--outfile demux_

cd $DIR/fastq/gex
$DEML/deML --index $DIR/samplesheet_gex.txt \
-f V350145982_S1_L001_R1_001.fastq.gz -r V350145982_S1_L001_R2_001.fastq.gz \
-if1 V350145982_S1_L001_I1_001.fastq.gz -if2 V350145982_S1_L001_I2_001.fastq.gz \
--summary demux_stats.txt \
--outfile demux_
