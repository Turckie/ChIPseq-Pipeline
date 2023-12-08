#!/bin/bash

# USAGE: sh EPIC3_IDR_pseudoreplicates.sh <input BAM rep1> <chip BAM rep1> <input BAM rep2> <chip BAM rep2> <NAME for IDR output>
# bsub -q normal -R "rusage[mem=10240]" -M 12288 ./EPIC_3_IDR_pseudoreplicates.sh 3546_E_run498_500.uniq.bam 3546_A_run498_500.uniq.bam 3546_F_run498_500.uniq.bam 3546_B_run498_500.uniq.bam TRB2
# bsub -q normal -R "rusage[mem=10240]" -M 12288 ./EPIC_3_IDR_pseudoreplicates.sh 3546_E_run498_500.uniq.bam 3546_C_run498_500.uniq.bam 3546_F_run498_500.uniq.bam 3546_D_run498_500.uniq.bam TRB3
# bsub -q multicore20 -n 4 -R "rusage[mem=20480]" -M 24576 ./EPIC_3_IDR_pseudoreplicates.sh 1255_A_run177_AAGGGAAT_L006_all.uniq.bam 1255_B_run177_CCTTCAAT_L006_all.uniq.bam 1255_D_run177_TTCAGCAT_L006_all.uniq.bam 1255_E_run177_AAGACGAT_L006_all.uniq.bam TRB1

#variables and paths
working_directory=$(pwd)
path_alignments=${working_directory}/data/alignments_ChIP-seq/unique
path_reference=${working_directory}/data/reference
path_results=${working_directory}/results/ChIP-seq
path_epic=${working_directory}/data/ChIP-seq/EPIC-IDR

pseudo=${path_epic}/PSEUDO

bs=100
fs=180

libraries=(ERR833055 ERR833056 ERR11556841 ERR11556842 ERR11556843 ERR11556844)
TRB1=(ERR833055 ERR833056)
TRB2=(ERR11556841 ERR11556842)
TRB3=(ERR11556843 ERR11556844)

samples=(TRB1 TRB2 TRB3)






#to evaluate the quality of the dataset, a IDR analysis is performed with pseudoreplicates and with internal replicates of each data set. As a rule of thumb suggested by the encode project, the difference in the muber of peaks identified by pseuoreplication or true replication should not be more than 2-fold. 
# Make Directories
mkdir -p ${pseudo}

echo "randomly splitting each bam file in half"
if [ ! -f "{path_alignments}/subsample"*".uniq.bam" ]; then
for i in ${libraries[@]}
do
${working_directory}/scripts/randomsplitbam.py -d ${path_alignments} -i ${i}.uniq.bam -f 0.5
samtools index ${path_alignments}/${i}.uniq.bam
done
fi



#merge subsamples from different replicates to Pseudoreplicates
#TRB1
echo "merging pseudoreplicates for TRB1"
for i in subsample1_0.5 subsample2_0.5
do
samtools cat ${path_alignments}/${i}_ERR833055.uniq.bam ${path_alignments}/${i}_ERR833056.uniq.bam | \
samtools sort - -o ${path_alignments}/${i}_TRB1.uniq.bam
samtools index ${path_alignments}/${i}_TRB1.uniq.bam
done

#TRB2
#echo "merging pseudoreplicates for TRB2"
for i in subsample1_0.5 subsample2_0.5
do
samtools cat ${path_alignments}/${i}_ERR11556841.uniq.bam ${path_alignments}/${i}_ERR11556842.uniq.bam | \
samtools sort - -o ${path_alignments}/${i}_TRB2.uniq.bam
samtools index ${path_alignments}/${i}_TRB2.uniq.bam
done

#TRB3
#echo "merging pseudoreplicates for TRB3"
for i in subsample1_0.5 subsample2_0.5
do
samtools cat ${path_alignments}/${i}_ERR11556843.uniq.bam ${path_alignments}/${i}_ERR11556844.uniq.bam | \
samtools sort - -o ${path_alignments}/${i}_TRB3.uniq.bam
samtools index ${path_alignments}/${i}_TRB3.uniq.bam
done


#Peak calling on pseudoreplicates
echo "Calling peaks for pseudoreplicates "

for r in subsample1_0.5 subsample2_0.5
 do
#macs2 callpeak -t ${tmpDir}/${TRB}.${r}.uniq.bam -c path_bam_unique/controls.uniq.bam -f BAM -g hs -n $macsDir/${NAME1}_pr -B -p 1e-3  2> $macsDir/${NAME1}_pr_macs2.log
epic2 --treatment ${path_alignments}/${r}_TRB1.uniq.bam --control ${path_alignments}/controlsTRB1.uniq.bam --genome Tair10 --bin-size ${bs} --fragment-size ${fs} --gaps-allowed 0 --effective-genome-fraction 1 --chromsizes ${path_reference}/TAIR10.chrlength_nuclear.txt --output ${pseudo}/${r}_TRB1
epic2 --treatment ${path_alignments}/${r}_TRB2.uniq.bam --control ${path_alignments}/controlsTRB23.uniq.bam --genome Tair10 --bin-size ${bs} --fragment-size ${fs} --gaps-allowed 0 --effective-genome-fraction 1 --chromsizes ${path_reference}/TAIR10.chrlength_nuclear.txt --output ${pseudo}/${r}_TRB2
epic2 --treatment ${path_alignments}/${r}_TRB3.uniq.bam --control ${path_alignments}/controlsTRB23.uniq.bam --genome Tair10 --bin-size ${bs} --fragment-size ${fs} --gaps-allowed 0 --effective-genome-fraction 1 --chromsizes ${path_reference}/TAIR10.chrlength_nuclear.txt --output ${pseudo}/${r}_TRB3
done

#Peak calling on subsamples
echo "Calling peaks for subsamples"

for i in ${TRB1[@]}
do
 for r in subsample1_0.5 subsample2_0.5
 do
#macs2 callpeak -t ${path_alignments}/${TRB}.${r}.uniq.bam -c path_bam_unique/controls.uniq.bam -f BAM -g hs -n $macsDir/${NAME1}_pr -B -p 1e-3  2> $macsDir/${NAME1}_pr_macs2.log
epic2 --treatment ${path_alignments}/${r}_${i}.uniq.bam --control ${path_alignments}/controlsTRB1.uniq.bam --genome Tair10 --bin-size ${bs} --fragment-size ${fs} --gaps-allowed 0 --effective-genome-fraction 1 --chromsizes ${path_reference}/TAIR10.chrlength_nuclear.txt --output ${pseudo}/${r}_${i}
 done
done

for i in ${TRB2[@]}
do
 for r in subsample1_0.5 subsample2_0.5
 do
#macs2 callpeak -t ${tmpDir}/${TRB}.${r}.uniq.bam -c path_bam_unique/controls.uniq.bam -f BAM -g hs -n $macsDir/${NAME1}_pr -B -p 1e-3  2> $macsDir/${NAME1}_pr_macs2.log
 epic2 --treatment ${path_alignments}/${r}_${i}.uniq.bam --control ${path_alignments}/controlsTRB23.uniq.bam --genome Tair10 --bin-size ${bs} --fragment-size ${fs} --gaps-allowed 0 --effective-genome-fraction 1 --chromsizes ${path_reference}/TAIR10.chrlength_nuclear.txt --output ${pseudo}/${r}_${i}
 done
done
 
for i in ${TRB3[@]}
do
 for r in subsample1_0.5 subsample2_0.5
 do
#macs2 callpeak -t ${tmpDir}/${TRB}.${r}.uniq.bam -c path_bam_unique/controls.uniq.bam -f BAM -g hs -n $macsDir/${NAME1}_pr -B -p 1e-3  2> $macsDir/${NAME1}_pr_macs2.log
 epic2 --treatment ${path_alignments}/${r}_${i}.uniq.bam --control ${path_alignments}/controlsTRB23.uniq.bam --genome Tair10 --bin-size ${bs} --fragment-size ${fs} --gaps-allowed 0 --effective-genome-fraction 1 --chromsizes ${path_reference}/TAIR10.chrlength_nuclear.txt --output ${pseudo}/${r}_${i}
 done
done

#IDR on pseudoreplicates
echo "Running IDR on pseudoreplicates..."
for TRB in  ${samples[@]}
do
 awk 'BEGIN{FS=OFS="\t"} NR>=2 {print $1,$2,$3,-(log($4)/log(10)), $5,$6,$7,$8,$9,$10}' ${pseudo}/subsample1_0.5_${TRB} | bedtools subtract -A -a - -b ${path_results}/final.blacklist.bed | bedtools intersect -wa -a - -b ${path_reference}/TAIR10.chr_nuclear.bed >${pseudo}/tmp-R1
 awk 'BEGIN{FS=OFS="\t"} NR>=2 {print $1,$2,$3,-(log($4)/log(10)), $5,$6,$7,$8,$9,$10}' ${pseudo}/subsample2_0.5_${TRB} | bedtools subtract -A -a - -b ${path_results}/final.blacklist.bed| bedtools intersect -wa -a - -b ${path_reference}/TAIR10.chr_nuclear.bed >${pseudo}/tmp-R2
 idr --samples ${pseudo}/tmp-R1 ${pseudo}/tmp-R2 --input-file-type bed  --rank 4 --peak-merge-method sum --output-file ${pseudo}/Pseudo_${TRB}-IDR.pValueBL --plot
 awk 'BEGIN{FS=OFS="\t"} ($5 >= 540) {print $0}' ${pseudo}/Pseudo_${TRB}-IDR.pValueBL >${pseudo}/Pseudo_${TRB}-IDR.SigpValueBL
 wc -l ${pseudo}/Pseudo_${TRB}-IDR.SigpValueBL | awk -v awkvar="$TRB" 'BEGIN {OFS="\t"} {print awkvar"_Np", $1}' - >>${path_results}/IDR-result.table.txt
 rm ${pseudo}/tmp-R1 ${pseudo}/tmp-R2
done

#IDR on subsamples
echo "Running IDR on subsamples..."
for TRB in ${libraries[@]}
do
 awk 'BEGIN{FS=OFS="\t"} NR>=2 {print $1,$2,$3,-(log($4)/log(10)), $5,$6 ,$7,$8,$9,$10}' ${pseudo}/subsample1_0.5_${TRB} | bedtools subtract -A -a - -b ${path_results}/final.blacklist.bed | bedtools intersect -wa -a - -b ${path_reference}/TAIR10.chr_nuclear.bed >${pseudo}/tmp-R1
 awk 'BEGIN{FS=OFS="\t"} NR>=2 {print $1,$2,$3,-(log($4)/log(10)), $5,$6 ,$7,$8,$9,$10}' ${pseudo}/subsample2_0.5_${TRB} | bedtools subtract -A -a - -b ${path_results}/final.blacklist.bed | bedtools intersect -wa -a - -b ${path_reference}/TAIR10.chr_nuclear.bed >${pseudo}/tmp-R2
 idr --samples ${pseudo}/tmp-R1 ${pseudo}/tmp-R2 --input-file-type bed  --rank 4 --peak-merge-method sum --output-file ${pseudo}/Pseudo_${TRB}-IDR.pValueBL --plot
 awk 'BEGIN{FS=OFS="\t"} ($5 >= 540) {print $0}' ${pseudo}/Pseudo_${TRB}-IDR.pValueBL >${pseudo}/Pseudo_${TRB}-IDR.SigpValueBL
 wc -l ${pseudo}/Pseudo_${TRB}-IDR.SigpValueBL | awk -v awkvar="$TRB" 'BEGIN {OFS="\t"} {print awkvar"_N", $1}' - >>${path_results}/IDR-result.table.txt
 rm ${pseudo}/tmp-R1 ${pseudo}/tmp-R2
done

