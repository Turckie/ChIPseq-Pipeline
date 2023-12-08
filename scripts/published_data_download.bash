#!/bin/bash
#helperscrpt to download pbulished chromatin data


working_directory=$(pwd)
path_published_data=${working_directory}/data/published_data_for_analysis/

mkdir -p ${path_published_data}

cd ${path_published_data}
#CLF and SWN peaks
#https://onlinelibrary.wiley.com/action/downloadSupplement?doi=10.1002%2Fpld3.100&file=pld3100-sup-0003-DataS1.xlsx
#bed file extracted by hand

#H3K27me3 data from https://academic.oup.com/plcell/article/28/1/87/6098213#supplementary-data
#bed file extracted by hand


#H3K4me3 data from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE204681 (Jacobsen group paper on JMJ14 to TRB connection)
wget https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM6185nnn/GSM6185711/suppl/GSM6185711_H3K4me3-ChIPseq-Col0-Rep1_peaks.narrowPeak.gz
wget https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM6185nnn/GSM6185711/suppl/GSM6185708_H3K4me3-ChIPseq-Col0-Rep2_peaks.narrowPeak.gz

idr --samples GSM6185711_H3K4me3-ChIPseq-Col0-Rep1_peaks.narrowPeak GSM6185708_H3K4me3-ChIPseq-Col0-Rep2_peaks.narrowPeak --input-file-type narrowPeak  --peak-merge-method sum --output-file H3K4me3-IDR.pValueBL --output-file-type bed --plot
