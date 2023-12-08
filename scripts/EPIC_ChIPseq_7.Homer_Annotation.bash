#!/bin/bash
######################################################################
 ## File: run this as bsub -q normal -R "rusage[mem=10240]" -M 12288 ./EPIC_ChIPseq_7.Homer_Annotation.bash
 ## Description: A new Homer annotation file for tair10 and Araprot11 was created using the homer installment in miniconda
 ## Date: 2021-05-19
 ## Author: Franziska Turck
######################################################################
 #####################################################################
 #  Author Modify description
 #2021-05-19 V01.000 Adapt script for IDR filter
 
######################################################################
 #variables
working_directory=$(pwd)
path_custom_araport11=${MAMBA_ROOT_PREFIX}/envs/TRB_Krause_et_al/share/homer/data/genomes/Araport11/annotations/custom
path_alignments=${working_directory}/data/alignments_ChIP-seq/unique
path_reference=${working_directory}/data/reference
path_results=${working_directory}/results/ChIP-seq
path_epic=${working_directory}/data/ChIP-seq/EPIC-IDR

genome_version=(Araport11)
declare -A ann
ann[Araport11]=${path_custom_araport11}/Araport11.txt

samples=(TRB3 TRB1 TRB2 merged TRB1_unique TRB2_unique TRB3_unique TRB1_2_private TRB1_3_private TRB2_3_private TRB123_common)
##############################################################################################################

cd ${path_epic}
mkdir -p HOMER

for i in  ${samples[@]}
do
 for a in  ${genome_version[@]}
 do
 awk '{OFS="\t"};{print $1,$2,$3,$1":"$2"-"$3"_IDRscore="$5,"NA","+"}' ${i}.bed| annotatePeaks.pl - ${a} -ann ${ann[${a}]} -annStats HOMER/${i}.${a}.stats | cut -f2,3,4,8 >HOMER/Annotate.${a}.${i}.txt
###############################################################################################################
#making a list file with AGIs
 cut -f4 HOMER/Annotate.${a}.${i}.txt | grep -o '\(A[TG0-9]*\)' | cut -c-9 | sed "1 d" |sort | uniq >HOMER/${a}.${i}.AGI
###############################################################################################################
#making a bed-like file as final file for the annotation
 awk -F"\t" 'BEGIN {OFS="\t"}{if (/promoter-TSS/) print $1,$2,$3,substr( $4, index($4,"A"), 9 ),"promoter-TSS",$4;
 	else if (/TTS-down/) print $1,$2,$3,substr( $4, index($4,"A"), 9 ),"TTS-down",$4 ;
 	else if (/exon/) print $1,$2,$3,substr( $4, index($4,"A"), 9 ),"genebody",$4;
 	else if (/intron/) print $1,$2,$3,substr( $4, index($4,"A"), 9 ),"genebody",$4 ;
 	else if (/1kb-promoter/) print $1,$2,$3,substr( $4, index($4,"A"), 9 ),"1kb-promoter",$4 ;
 	else if (/3kb-promoter/) print $1,$2,$3,substr( $4, index($4,"A"), 9 ),"3kb-promoter",$4 ;
 	else print $1,$2,$3,substr( $4, index($4,"A"), 9 ),"Intergenic",$4 ;}' HOMER/Annotate.${a}.${i}.txt >HOMER/${a}.${i}.final
awk -F"\t" -v s=1 'BEGIN {FS=OFS="\t"} (NR!=1) {print $1,$2-s,$3,$5,$4}' HOMER/${a}.${i}.final > HOMER/${a}.${i}.bed
rm >HOMER/Annotate.${a}.${i}.txt
###############################################################################################################
#counting annotation elements
 wc -l HOMER/${a}.${i}.AGI >>HOMER/Homer.analysis.txt
 awk -F"\t" '$5=="genebody" { count++ } END { print "genebody" , count }' HOMER/${a}.${i}.final >>HOMER/Homer.analysis.txt
 awk -F"\t" '$5=="promoter-TSS" { count++ } END { print "promoter-TSS" , count }' HOMER/${a}.${i}.final >>HOMER/Homer.analysis.txt
 awk -F"\t" '$5=="TTS-down" { count++ } END { print "TTS" , count }' HOMER/${a}.${i}.final >>HOMER/Homer.analysis.txt
 awk -F"\t" '$5=="1kb-promoter" { count++ } END { print "1kb-promoter" , count }' HOMER/${a}.${i}.final >>HOMER/Homer.analysis.txt
 awk -F"\t" '$5=="3kb-promoter" { count++ } END { print "3kb-promoter" , count }' HOMER/${a}.${i}.final >>HOMER/Homer.analysis.txt
 awk -F"\t" '$5=="Intergenic" { count++ } END { print "Intergenic" , count }' HOMER/${a}.${i}.final >>HOMER/Homer.analysis.txt
 done
done

rm ${path_results}/gene.category.counts.txt


for i in ${samples[@]}
do
wc -l HOMER/Araport11.${i}.AGI >>${path_results}/gene.category.counts.txt
wc -l HOMER/Araport11.${i}.bed >>${path_results}/peak.category.counts.txt
done
cd ${working_directory}