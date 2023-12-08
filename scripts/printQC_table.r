#!/usr/bin/Rscript

library(tidyverse)
library(ggpubr)
library(forcats)
library("ggsci")
library("gridExtra")
library(reshape2)
library(gtools)
library(formattable)
library(data.table)


Mapping_sum_key <- read.delim("results/RNA-seq/Mapping_sum_key.txt")


seedlings<-ggplot(subset(Mapping_sum_key,type!="Input_reads"& tissue=="seedling"), aes(fct_inorder(sample)))+
geom_bar(aes(y = reads, fill = type), stat = "identity", position = "stack")+
labs(title = "seedlings",y = "Reads")+
scale_fill_npg()+
theme_classic()+
theme(axis.text.x = element_text(angle = 90),legend.position = "none",axis.title.x=element_blank())

pdf("results/RNA-seq/RNA-seq_Mapped_reads.pdf")
ggarrange(seedlings,
          labels = c("A"),
          common.legend = TRUE, legend="bottom",
          ncol = 1, nrow = 1)
dev.off()

#Make a complete table for a supplementary file:
table<-reshape2::dcast(Mapping_sum_key, sample+run+read_accession+tissue+genotype~type, value.var =c("reads"),sum )
table$fraction_mapped<-((table$Uniquely_mapped_reads+table$Many_mapped_reads+table$Multi_mapped_reads)/table$Input_reads)
table$fraction_mapped<-percent(table$fraction_mapped)
#this takes care of the problem with mixed letters and numbers for ordering, comes from gtools
table = data.table(table)
table<-table[gtools::mixedorder(as.character(sample))]
#this writes all the data to a plain txt file
write.table(table,"results/RNA-seq/RNA-seq_mapping_overview_table.txt", quote=FALSE, sep="\t", row.names=FALSE)


#this makes a pretty html version of the table, but only really works in interactive mode:

pretty<-formattable(table, list(
  fraction_mapped = color_tile("white", "orange"),
  area(col = c(Uniquely_mapped_reads, Multi_mapped_reads, Many_mapped_reads)) ~ normalize_bar("pink", 0.2)
))

