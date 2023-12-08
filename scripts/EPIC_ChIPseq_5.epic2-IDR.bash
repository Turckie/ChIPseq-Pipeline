#!/bin/bash
######################################################################
 ## Author: Franziska Turck
######################################################################

#paths
working_directory=$(pwd)
path_alignments=${working_directory}/data/alignments_ChIP-seq/unique
path_reference=${working_directory}/data/reference
path_results=${working_directory}/results/ChIP-seq
path_epic=${working_directory}/data/ChIP-seq/EPIC-IDR

mkdir -p ${path_epic}

#variables
control=controls
IP1=ERR833055
IP2=ERR833056
IP3=ERR11556841
IP4=ERR11556842
IP5=ERR11556843
IP6=ERR11556844
bs=100
fs=180
samples=(TRB1 TRB2 TRB3)
#merge Col-0 IPs as controls for TRB1 and TRB2/3

cd ${path_alignments}

#if [ ! -f "controlsTRB1.uniq.bam" ]; then
#samtools cat ERR833053.uniq.bam ERR833054.uniq.bam | samtools view -bhs 42.5 -  | samtools sort - -o controlsTRB1.uniq.bam
#samtools index controlsTRB1.uniq.bam
#fi

#if [ ! -f "controlsTRB23.uniq.bam" ]; then
#samtools cat ERR11556839.uniq.bam ERR11556840.uniq.bam | samtools view -bhs 42.5 -  | samtools sort - -o controlsTRB23.uniq.bam
#samtools index controlsTRB23.uniq.bam
#fi

#Predicting peaks using merged Col-ChIP replicates as control
echo "predict TRB1, TRB2 and TRB3 ChIP-seq"

#TRB1-R1
epic2 --treatment ${IP1}.uniq.bam --control controlsTRB1.uniq.bam --genome Tair10 --bin-size $bs --fragment-size $fs --gaps-allowed 0 --effective-genome-fraction 1 --chromsizes ${path_reference}/TAIR10.chrlength_nuclear.txt --output $path_epic/TRB1-R1

#TRB1-R2
epic2 --treatment ${IP2}.uniq.bam --control controlsTRB1.uniq.bam --genome Tair10 --bin-size $bs --fragment-size $fs --gaps-allowed 0 --effective-genome-fraction 1 --chromsizes ${path_reference}/TAIR10.chrlength_nuclear.txt --output ${path_epic}/TRB1-R2


#TRB2-R1
epic2 --treatment ${IP3}.uniq.bam --control controlsTRB23.uniq.bam --genome Tair10 --bin-size $bs --fragment-size $fs --gaps-allowed 0 --effective-genome-fraction 1 --chromsizes ${path_reference}/TAIR10.chrlength_nuclear.txt --output ${path_epic}/TRB2-R1

#TRB2-R2
epic2 --treatment ${IP4}.uniq.bam --control controlsTRB23.uniq.bam --genome Tair10 --bin-size $bs --fragment-size $fs --gaps-allowed 0 --effective-genome-fraction 1 --chromsizes ${path_reference}/TAIR10.chrlength_nuclear.txt --output ${path_epic}/TRB2-R2

#TRB3-R1
epic2 --treatment ${IP5}.uniq.bam --control controlsTRB23.uniq.bam --genome Tair10 --bin-size $bs --fragment-size $fs --gaps-allowed 0 --effective-genome-fraction 1 --chromsizes ${path_reference}/TAIR10.chrlength_nuclear.txt --output ${path_epic}/TRB3-R1

#TRB3-R2
epic2 --treatment ${IP6}.uniq.bam --control controlsTRB23.uniq.bam --genome Tair10 --bin-size $bs --fragment-size $fs --gaps-allowed 0 --effective-genome-fraction 1 --chromsizes ${path_reference}/TAIR10.chrlength_nuclear.txt --output ${path_epic}/TRB3-R2

###############################################################################################################################################################################################################################################################
#we want to filter the reads against the blacklist that was previously generated
#we would also somehow overlap the replicates in a meaningful was - this will be done using the IDR method implemented in the ChIP-seq pipeline
#from Encode.

cd ${working_directory}
echo "running IDR for experimental replicates"

for i in ${samples[@]}
do
 awk 'BEGIN{FS=OFS="\t"} NR>=2 {print $1,$2,$3,$1":"$2".."$3_"p"_$4_"fc"_$10,-(log($4)/log(10)),$6,$7,$8,$9,$10}' ${path_epic}/${i}-R1 | bedtools subtract -A -a - -b ${path_results}/final.blacklist.bed | bedtools intersect -wa -a - -b ${path_reference}/TAIR10.chr_nuclear.bed >tmp-R1
 awk 'BEGIN{FS=OFS="\t"} NR>=2 {print $1,$2,$3,$1":"$2".."$3_"p"_$4_"fc"_$10,-(log($4)/log(10)),$6,$7,$8,$9,$10}' ${path_epic}/${i}-R2 | bedtools subtract -A -a - -b ${path_results}/final.blacklist.bed | bedtools intersect -wa -a - -b ${path_reference}/TAIR10.chr_nuclear.bed >tmp-R2
 idr --samples ${path_epic}/tmp-R1 ${path_epic}/tmp-R2 --input-file-type bed  --rank 5 --peak-merge-method sum --output-file ${path_epic}/${i}-IDR.pValueBL --output-file-type bed --plot
 awk 'BEGIN{FS=OFS="\t"} ($5 >= 540) {print $0}' ${path_epic}/${i}-IDR.pValueBL > ${path_epic}/${i}-IDR.SigpValueBL
 awk '{OFS="\t"};{print $1,$2,$3,$1":"$2"-"$3"_IDRscore="$5,"NA","+"}' ${path_epic}/${i}-IDR.SigpValueBL > ${path_epic}/${i}.bed
 wc -l ${path_epic}/${i}-IDR.SigpValueBL | awk -v awkvar="$i" 'BEGIN {OFS="\t"} {print awkvar"_Nt", $1}' - >> ${path_results}/IDR-result.table.txt
 rm tmp*
done

cd ${path_epic}
#Generate a list of peaks using bedtools merge 
cat TRB1-IDR.SigpValueBL TRB2-IDR.SigpValueBL TRB3-IDR.SigpValueBL | sort -k1,1 -k2,2g | bedtools merge -i - | bedtools intersect -wa -a - -b ${path_reference}/TAIR10.chr_nuclear.bed > merged.bed
#Annotate peaks in merged file a TRB category

num=0
for i in ${samples[@]}
do
 num=$(( ${num} + 1 ))
 bedtools intersect -wa -c -a merged-IDR.SigpValueBL -b ${i}-IDR.SigpValueBL | cut -f4  > ${num}-out
done

paste merged-IDR.SigpValueBL *-out > tmp.txt
rm *-out

awk 'BEGIN {OFS=FS="\t"} $4>=1 && $5==0 && $6==0 {print $1,$2,$3,$4"-"$5"-"$6}' tmp.txt > TRB1_unique.bed
awk 'BEGIN {FS=OFS="\t"} $4==0 && $5>=1 && $6==0 {print $1,$2,$3,$4"-"$5"-"$6}' tmp.txt > TRB2_unique.bed
awk 'BEGIN {OFS=FS="\t"} $4==0 && $5==0 && $6>=1 {print $1,$2,$3,$4"-"$5"-"$6}' tmp.txt > TRB3_unique.bed
awk 'BEGIN {OFS=FS="\t"} $4>=1 && $5>=1 && $6==0 {print $1,$2,$3,$4"-"$5"-"$6}' tmp.txt > TRB1_2_private.bed
awk 'BEGIN {OFS=FS="\t"} $4>=1 && $5==0 && $6>=1 {print $1,$2,$3,$4"-"$5"-"$6}' tmp.txt > TRB1_3_private.bed
awk 'BEGIN {OFS=FS="\t"} $4==0 && $5>=1 && $6>=1 {print $1,$2,$3,$4"-"$5"-"$6}' tmp.txt > TRB2_3_private.bed
awk 'BEGIN {OFS=FS="\t"} $4>=1 && $5>=1 && $6>=1 {print $1,$2,$3,$4"-"$5"-"$6}' tmp.txt > TRB123_common.bed

rm tmp.txt

cd ${working_directory}




####################################################################
echo "Program Finished"
exit 0
