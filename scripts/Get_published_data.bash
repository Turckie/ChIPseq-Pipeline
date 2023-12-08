#!/bin/bash
working_directory=$(pwd)
path_reads=${working_directory}/data/reads_ChIP-seq
path_reference=${working_directory}/data/reference
path_annotation=${working_directory}/data/annotation
path_published=${working_directory}/data/published_data_for_analysis

mkdir -p ${path_published}



echo "Downloading Supplemental Data Set 1 from Zhou et al. Nature Genetics volume 50, pages638–644 (2018)"

wget -P ${path_published} https://static-content.springer.com/esm/art%3A10.1038%2Fs41588-018-0109-9/MediaObjects/41588_2018_109_MOESM3_ESM.xlsx


echo  "Downloading Supplemental Data Set 3 from Zhou et al. Nature Genetics volume 50, pages638–644 (2018)"

wget -P ${path_published} https://static-content.springer.com/esm/art%3A10.1038%2Fs41588-018-0109-9/MediaObjects/41588_2018_109_MOESM6_ESM.xlsx