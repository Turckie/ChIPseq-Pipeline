#!/bin/bash
######################################################################
 ## File: run this as ./sripts/Custom_annotation_Homer.sh
 ## Description: Create custom genome annotation for Homer in the /share/homer placed in the micromamba or miniconda environment
 ## Date: 2021-04-28
 ## Author: Franziska Turck
######################################################################
 #####################################################################
 #  Author Modify description
 # Available genomes can be installed with a perl script found in ${install_path}/micromamba/envs/TRB_Krause_et_al/share/homer
 # TAIR10 annotation is available but not the newer Araport11
 #  The normal basic Homer program adds too much sequence to the promoter-TSS and TTS regions for it to make any sense for the small Arabidopsis genome
 # Therfore, this script generates a custom annotation file for homer, for Araport11 and for tair10. Since the default tair10 annotation contains many non-coding RNAs 
 # that may be true or not, we also add another custom file for TAIR10.
 #  These custom genomes declare TTS and TTS as plu/minus 100pb, 1kb promoter as -1000 to -100 from TSS and 3kb promoter from -3000 to -1000kb form TSS
 #..TAIR10 and tair10 should be the same, but the tair10 version has non-coding RNAs annotated that may or may not be true
#################################################################################################################################
#type which micromamba to find out the prefix to your micromamba
working_directory=$(pwd)
path_annotation=${working_directory}/data/annotation
path_custom_araport11=${MAMBA_ROOT_PREFIX}/envs/TRB_Krause_et_al/share/homer/data/genomes/Araport11/annotations/custom
#######################################################################################################################################################
#download the Aralort11 gtf file and save it in data/annotation
mkdir -p ${path_annotation}
wget -P ${path_annotation} https://www.arabidopsis.org/download_files/Genes/Araport11_genome_release/Araport11_GTF_genes_transposons.current.gtf.gz
gzip -d  ${path_annotation}/Araport11_GTF_genes_transposons.current.gtf.gz


#######################################################################################################################################################
#first install the preformated tair10 instance 
cd ${MAMBA_ROOT_PREFIX}/envs/TRB_Krause_et_al/share/homer
perl ./configureHomer.pl -install tair10
perl ./configureHomer.pl -install arabidopsis-0
perl ./configureHomer.pl -install arabidopsis-p
##############################################################################################################################################

#Araport11
#load the genome from the data folder to the correct homer repository

loadGenome.pl -force -name Araport11 -org arabidopsis -fasta ${working_directory}/data/reference/TAIR10_chr_all.fas -gtf ${working_directory}/data/annotation/Araport11_GTF_genes_transposons.current.gtf


#This parses a custom file that considers promoter-TSS and TTS 100bp plus minus start and end
parseGTF.pl ${working_directory}/data/annotation/Araport11_GTF_genes_transposons.current.gtf  ann -annTSSstartOffset -100 -annTSSendOffset 100 -annTTSstartOffset -100 -annTTSendOffset 100 > ${path_custom_araport11}/Araport11.annotations.txt

#now create a file that contains 1kb promoters 
parseGTF.pl ${working_directory}/data/annotation/Araport11_GTF_genes_transposons.current.gtf  ann -annTSSstartOffset -1000 -annTSSendOffset -100 |\
awk 'BEGIN {FS=OFS="\t"} $1 ~ /promoter-TSS/ {print $1,$2,$3,$4,$5,"1kb";}' -|\
sed 's/promoter-TSS/1kb-promoter/g' - > ${path_custom_araport11}/Araport11.1kbPromoter.annotations.txt
#now create a file that contains potential distal 3kb promoters 
parseGTF.pl ${working_directory}/data/annotation/Araport11_GTF_genes_transposons.current.gtf  ann -annTSSstartOffset -3000 -annTSSendOffset -1000 |\
awk 'BEGIN {FS=OFS="\t"} $1 ~ /promoter-TSS/ {print $1,$2,$3,$4,$5,"3kb";}' -|\
sed 's/promoter-TSS/3kb-promoter/g' - > ${path_custom_araport11}/Araport11.kbPromoter.annotations.txt
#extract all other annotations from the custom text file 

awk 'BEGIN {FS=OFS="\t"} {if (/exon/) { print $1,$2,$3,$4, $5,$6;}} ' ${path_custom_araport11}/Araport11.annotations.txt >${path_custom_araport11}/Araport11.exon.annotations.txt
awk 'BEGIN {FS=OFS="\t"} {if (/intron/) { print $1,$2,$3,$4, $5,$6;}} ' ${path_custom_araport11}/Araport11.annotations.txt >${path_custom_araport11}/Araport11.intron.annotations.txt
awk 'BEGIN {FS=OFS="\t"} {if (/promoter-TSS/) { print $1,$2,$3,$4, $5,$6;}} ' ${path_custom_araport11}/Araport11.annotations.txt >${path_custom_araport11}/Araport11.promoter-TSS.annotations.txt
awk 'BEGIN {FS=OFS="\t"} {if (/TTS/) { print $1,$2,$3,$4, $5,$6;}} ' ${path_custom_araport11}/Araport11.annotations.txt | sed 's/TTS/TTS-down/g' >${path_custom_araport11}/Araport11.TTS.annotations.txt

#now paste them together in the order in which the annotation should be assigned
cat ${path_custom_araport11}/Araport11.promoter-TSS.annotations.txt \
    ${path_custom_araport11}/Araport11.exon.annotations.txt \
    ${path_custom_araport11}/Araport11.intron.annotations.txt \
    ${path_custom_araport11}/Araport11.TTS.annotations.txt \
    ${path_custom_araport11}/Araport11.1kbPromoter.annotations.txt \
    ${path_custom_araport11}/Araport11.3kbPromoter.annotations.txt > ${path_custom_araport11}/Araport11.txt

assignGenomeAnnotation ${path_custom_araport11}/Araport11.txt ${path_custom_araport11}/Araport11.txt -prioritize ${path_custom_araport11}/Araport11.txt > stats.txt


