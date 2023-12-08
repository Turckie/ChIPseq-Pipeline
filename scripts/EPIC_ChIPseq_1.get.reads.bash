#!/bin/bash

working_directory=$(pwd)
path_reads=${working_directory}/data/reads_ChIP-seq
path_reference=${working_directory}/data/reference
path_annotation=${working_directory}/data/annotation


mkdir -p ${path_reads}
mkdir -p ${path_reference}
mkdir -p ${path_annotation}

echo "downloading the read files from the European Nucleotide Archive"

file="$1"
#which column has the header "fastq_ftp"?
link=$(awk -F ',' 'NR==1{for(i=1;i<=NF;i++){if($i=="fastq_ftp"){col=i;break}}}{print $col}' "${working_directory}/${file}")
#now download the list of files:
while read line; do
adress=$(echo $line | cut -d '\t' -f ${link})
wget -P ${path_reads} $adress
done <(tail -n +2 "${working_directory}/${file}")

echo "done downloading the read files from the European Nucleotide Archive"

echo "downloading the Arabidopsis thaliana reference genome"

wget -P ${path_reference} https://www.arabidopsis.org/download_files/Genes/TAIR10_genome_release/TAIR10_chromosome_files/TAIR10_chr_all.fas.gz
gzip -d ${path_reference}/TAIR10_chr_all.fas.gz
samtools faidx ${path_reference}/TAIR10_chr_all.fas
awk 'BEGIN {FS="\t"}; {print $1 FS "0" FS $2}' ${path_reference}/TAIR10_chr_all.fas.fai > ${path_reference}/TAIR10_chr_all.bed
awk 'BEGIN {FS="\t"}; {print $1 FS $2}' ${path_reference}/TAIR10_chr_all.fas.fai > ${path_reference}/TAIR10.chrlength_all.txt
#if we want to focus on the nuclear genome
awk 'BEGIN {FS="\t"}; ($1 ~ /Chr[1-5]/) {print $1 FS "0" FS $2}' ${path_reference}/TAIR10_chr_all.fas.fai > ${path_reference}/TAIR10.chr_nuclear.bed
awk 'BEGIN {FS="\t"}; ($1 ~ /Chr[1-5]/) {print $1 FS $2}' ${path_reference}/TAIR10_chr_all.fas.fai > ${path_reference}/TAIR10.chrlength_nuclear.txt

echo "done downloading the Arabidopsis thaliana reference genome"

echo "downloading the Arabidopsis genome annotation file Araport11"

wget -P ${path_annotation} https://www.arabidopsis.org/download_files/Genes/Araport11_genome_release/Araport11_GFF3_genes_transposons.current.gff.gz
gzip -d  ${path_annotation}/Araport11_GFF3_genes_transposons.current.gff.gz
wget -P ${path_annotation} https://www.arabidopsis.org/download_files/Genes/Araport11_genome_release/Araport11_GTF_genes_transposons.current.gtf.gz
gzip -d  ${path_annotation}/Araport11_GTF_genes_transposons.current.gtf.gz

echo "done downloading the Arabidopsis genome annotation file Araport11"
