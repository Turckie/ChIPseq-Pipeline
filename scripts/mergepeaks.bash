#!/bin/bash

working_directory=$(pwd)
path_reads=${working_directory}/data/reads_ChIP-seq
path_reference=${working_directory}/data/reference
path_annotation=${working_directory}/data/annotation
path_published=${working_directory}/data/published_data_for_analysis

#merge peak sets for CLF and SWN

cat ${path_published}/CLF.bed ${path_published}/SWN.bed |\
sort -k1,1 -k2,2g | bedtools merge -i - > ${path_published}/CLFSWN.bed