#!/bin/bash

###################################################################################################################################
#This script generates a blacklist of regions in Arabidopsis that should be excluded from the ChIP-seq analysis for various reasons.
#These regions may correspond to over-sampled or under-samples regions. A reason for both could be highly repetitive sequence or cryptic duplications.
#Another reason could be that the ChIP produces an artefact in the controls due to unspecific "chromatin_stickyness". Since a different lot
#of the same antibody could behave differnetly, we find that it is better to create a custom blacklist for each experiment rather than using some published lists.
#####################################################################################################################################

working_directory=$(pwd)
path_alignments=${working_directory}/data/alignments_ChIP-seq/unique
path_reference=${working_directory}/data/reference
path_results=${working_directory}/results/ChIP-seq

window=200

#generate a bedfile of 200bp windows across the Arabidopsis genome unless the file already exists

if [ ! -f "{path_reference}/TAIR10.${window}bp.windows.bed" ]; then
bedtools makewindows -g ${path_reference}/TAIR10.chrlength_nuclear.txt -w ${window} >${path_reference}/TAIR10.${window}bp.windows.bed
fi
#the following calculates coverage across the genome for 200pb windows, the script summarizes coverag from the controls and inputs.
#BAMscale cov -t 4 -o ${path_alignments} --bed ${path_reference}/TAIR10.${window}bp.windows.bed --bam ${path_alignments}/controls.uniq.bam --bam ${path_alignments}/inputs.uniq.bam

#Now the scripts uses R to identify outlyer windows
#in R
cd ${path_alignments}
Rscript ${working_directory}/scripts/blacklist.r


#The blacklist contains many small windows which are often adjacent and cover larger regions. These should be merged. In this setting, windows will be merged if the distance in not more than 200 bp.

bedtools merge -d 400 -i ${path_alignments}/blacklist.bed >${path_results}/final.blacklist.bed
rm ${path_alignments}/blacklist.bed
#sed 's/Chr/chr/g' ${path_results}/final.blacklist.bed >${path_results}/chr_final.blacklist.bed

cd ${working_directory}





