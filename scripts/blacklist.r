library(tidyverse)
#library(ggplot2)
#library(reshape2)
#library(dplyr)
#library(tidyr)


raw_coverages <- read_delim("Library_normalized_coverages.tsv", "\t", escape_double = FALSE, trim_ws = TRUE)
raw_coverages$sum<-(raw_coverages$controls.uniq.bam + raw_coverages$inputs.uniq.bam)
out <- boxplot.stats(raw_coverages$sum)$out
out_ind <- which(raw_coverages$sum %in% c(out))
inner_ind<-which(!raw_coverages$sum %in% c(out))
Outlier<-raw_coverages[out_ind, ]
Inner<-raw_coverages[inner_ind, ]
Outlier %>% separate(coordinate, c("chr", "pos"),":")->Outlier
Outlier %>% separate(pos, c("start", "end"),"-")->Outlier
write.table(Outlier[,c(1:3,6)],"blacklist.bed",row.names=FALSE,col.names=FALSE,sep="\t", eol="\n",quote = FALSE)


png("library_normalized_coverages.png")
boxplot(raw_coverages$sum, main="Coverage all windows", xlab="sum of library-normalized coverage")
dev.off()

png("library_normalized_coverages_after_blacklist.png")
boxplot(Inner$sum, main="Coverage all windows without blacklisted windows", xlab="sum of library-normalized coverage")
dev.off()

quit("no")