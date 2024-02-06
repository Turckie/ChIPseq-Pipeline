#!/bin/bash

working_directory=$(pwd)
path_reads=${working_directory}/data/fastq
path_QC=${working_directory}/data/fastq_QC

mkdir -p ${path_reads}
mkdir -p ${path_QC}

echo "downloading the read files with wget"

file="$1"
#which column has the "link to the fastq files?
col_num=$(awk -F '\t' 'NR==1{for(i=1;i<=NF;i++)if($i=="link"){print i;exit}}' "${working_directory}/${file}")

#now download the files the fastq files:
cut -f$col_num ${working_directory}/${file} | tail -n +2 | while read line
do
    # Execute the line
    echo ${line}
    wget --no-check-certificate -P ${path_reads} ${line}
done

#which column has the QC files
col_num=$(awk -F '\t' 'NR==1{for(i=1;i<=NF;i++)if($i=="QC"){print i;exit}}' "${working_directory}/${file}")

#now download the files the fastq files:
cut -f$col_num ${working_directory}/${file} | tail -n +2 | while read line
do
    # Execute the line
    echo ${line}
    wget --no-check-certificate -P ${path_QC} ${line}
done
#unzip all GC files and delete the zipped files
for i in ${path_QC}/*zip
do
unzip ${i}
rm ${i}
done
#check if download is complete

awk 'NR==1 {for (i=1; i<=NF; i++) {f[$i] = i}}
{print $(f["file"]), $(f["md5sum"]) }' Samples.txt >md5sum.txt

for i in ${path_reads}/*fastq.gz
do
md5sum -c md5sum.txt $i >>${path_reads}/md5_check.txt
done
