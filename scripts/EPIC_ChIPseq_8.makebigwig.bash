#!/bin/bash


working_directory=$(pwd)
path_results=${working_directory}/results/ChIP-seq
path_alignments=${working_directory}/data/alignments_ChIP-seq/unique
path_bigwig=${working_directory}/data/alignments_ChIP-seq/BIGWIG
#libraries=(ERR833055 ERR833056 ERR11556841 ERR11556842 ERR11556843 ERR11556844 controls controlsTRB1 controlsTRB23)
libraries=(ERR833056)

######################################################################################################################

mkdir -p ${path_bigwig}
                                                    
for i in ${libraries[@]}                                                  
do                                                              

bamCoverage -b ${path_alignments}/$i.uniq.bam -o ${path_bigwig}/${i}.bw -p 5 -bs 10 --blackListFileName ${path_results}/final.blacklist.bed --smoothLength 50 --extendReads 200 --ignoreDuplicates --normalizeUsing CPM
done
###############################################################################################################################
#namechange
file="Samples.ChIP-seq.txt"
while read line; do
old=$(echo $line | cut -d ' ' -f 9)
new=$(echo $line | cut -d ' ' -f 2)
mv ${path_bigwig}/${old}.bw ${path_bigwig}/${new}.bw
done < "${working_directory}/${file}"
                                                   
###############################################################################################################################
#subtract control from sample
for i in TRB1
do
 for r in R1 R2
 do
 bigwigCompare -b1 ${path_bigwig}/${i}_ChIP_${r}.bw -b2 ${path_bigwig}/controlsTRB1.bw --operation subtract -o ${path_bigwig}/${i}-control_ChIP_${r}.bw
 done
done

for i in TRB2 TRB3
do
 for r in R1 R2
 do
 bigwigCompare -b1 ${path_bigwig}/${i}_ChIP_${r}.bw -b2 ${path_bigwig}/controlsTRB23.bw --operation subtract -o ${path_bigwig}/${i}-control_ChIP_${r}.bw
 done
done