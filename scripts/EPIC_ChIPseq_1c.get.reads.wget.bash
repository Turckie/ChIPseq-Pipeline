#!/bin/bash

working_directory=$(pwd)
path_reads=${working_directory}/data/fastq

mkdir -p ${path_reads}


echo "downloading the read files with wget"

file="$1"
#now download the files the fastq files:
do
    # Execute the line
    echo ${line}
    wget --no-check-certificate -P ${path_reads} ${line}
done

