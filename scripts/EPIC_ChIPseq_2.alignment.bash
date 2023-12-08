#!/bin/bash

##############################################################################################################################
#This script aligns any compressed (ending with fastq.gz) fastq formated read file to a reference genome. The aligned reads are
#sorted and stored as .bam files. The script also creates a filtered version of the aligned reads so that ony those reads that map
#to an unique position to the genome are retained. Finally, the script counts the total reads, the mapped reads and the uniquely
#mapped reads and provides a text file with a summary.
##############################################################################################################################

working_directory=$(pwd)
path_reads=${working_directory}/data/reads_ChIP-seq
path_alignments=${working_directory}/data/alignments_ChIP-seq
path_reference=${working_directory}/data/reference
path_results=${working_directory}/results/ChIP-seq

#create the following structure of folders in the working directory  
#cd /biodata/dep_coupland/grp_turck/NGS_data_GC/GC_3546

mkdir -p ${path_alignments}/unique
mkdir -p ${path_results}


#the alignment is done with bwa

#If bwa is used for the first time on this genome, a genome index needs to be created. Uncomment and add the correct path to your genome fasta file.

#bwa index -p ${path_reference}/TAIR10bwaidx -a bwtsw ${path_reference}/TAIR10_chr_all.fas

#align with standard parameters, sort the aligned reads and save as .bam formated files
cd  ${path_reads}
for i in *fastq.gz
do
bwa mem -t 4 ${path_reference}/TAIR10bwaidx $i | samtools view -b - | samtools sort - -o ${path_alignments}/"${i%.fastq.gz}.bam" 
samtools index ${path_alignments}/"${i%.fastq.gz}.bam" 
done

#count the read and mapping statistics
cd ${path_alignments}

echo -e "file"'\t'"category"'\t'"number" > ${path_results}/Mapping.info.txt

for i in *.bam
do
samtools view -c ${i} > tmp
echo -e "${i}"'\t'"all_reads" >tmp2 
paste tmp2 tmp >>${path_results}/Mapping.info.txt
samtools view -c -F 260 ${i} > tmp
echo -e "${i}"'\t'"mapped_reads" >tmp2 
paste tmp2 tmp >>${path_results}/Mapping.info.txt
samtools view -c -q 10 ${i} >tmp
echo -e "${i}"'\t'"uniquely_mapped_reads" >tmp2
paste tmp2 tmp >>${path_results}/Mapping.info.txt
rm tmp tmp2
done




#filter the aligned reads to keep only uniquely aligned reads

for i in *.bam
do
samtools view -q 10 -b ${i} > ${path_alignments}/unique/${i%.bam}.uniq.bam
samtools index ${path_alignments}/unique/${i%.bam}.uniq.bam
done

cd ${working_directory}