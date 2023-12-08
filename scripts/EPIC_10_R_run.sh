#!/bin/bash
######################################################################
 ## File:  ./EPIC_10_R_run.sh
 ## Description: scrpt to obtain fasta sequences for a list of bedfiles
 ## Date: 2021-04-29
 ## Author: Franziska Turck
######################################################################
 #####################################################################
#!/usr/bin/Rscript


library("RColorBrewer")
library("ggplot2")
library("ggsci")
library("ggpubr")
library("reshape2")
library("dplyr")
library("readr")

Araport11 <- read.table("data/ChIP-seq/RMASTER/Araport11.txt","\t",  header = TRUE, sep = "\t")

#TAIR10 <- read_delim("TAIR10.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
KKfull <- Araport11
KKfull$telobox[KKfull$telobox>=6]<- "more than 5"
KKfull$RMCCTAR[KKfull$RMCCTAR>=6]<- "more than 5"
KKfull$CRACCTA[KKfull$CRACCTA>=6]<- "more than 5"


KKfull$TRB[KKfull$TRB1>=1 & KKfull$TRB2>=1 & KKfull$TRB3>=1 ]<- "TRB123 common"
KKfull$TRB[KKfull$TRB1>=1 & KKfull$TRB2>=1 & KKfull$TRB3==0 ]<- "TRB1/2 private"
KKfull$TRB[KKfull$TRB1>=1 & KKfull$TRB2==0 & KKfull$TRB3==0 ]<- "TRB1 unique"
KKfull$TRB[KKfull$TRB1>=1 & KKfull$TRB2==0 & KKfull$TRB3>=1 ]<- "TRB1/3 private"
KKfull$TRB[KKfull$TRB1==0 & KKfull$TRB2==0 & KKfull$TRB3>=1 ]<- "TRB3 unique"
KKfull$TRB[KKfull$TRB1==0 & KKfull$TRB2>=1 & KKfull$TRB3>=1 ]<- "TRB2/3 private"
KKfull$TRB[KKfull$TRB1==0 & KKfull$TRB2>=1 & KKfull$TRB3==0 ]<- "TRB2 unique"
KKfull$TRB[KKfull$TRB1==0 & KKfull$TRB2==0 & KKfull$TRB3==0 ]<- "non target"
KKfull$length<-as.numeric(KKfull$stop)-as.numeric(KKfull$start)
kksmall<-subset(KKfull, type=="peak")

TRBmelt = melt(KKfull, id.vars = c("type", "H3K27me3", "Open", "SWN", "CLF", "H2AK121ub1","H3K27ac", "H3K4me3", "H3K36me3", "cluster", "annotation","length"), measure.vars = c("TRB1", "TRB2", "TRB3"))

##########################################

pdf(file = "results/TRB_fragment_length.pdf", width = 4, height = 4)
ggplot(KKfull,aes(length, color=type))+geom_histogram( )
dev.off()
##########################################

#overlap TRBs with chromatin features
a<-ggplot(subset(KKfull, length<=1500), aes(H3K27me3,fill=as.factor(H3K27me3))) + geom_bar(aes(type), position = position_fill(reverse = TRUE)) + scale_fill_npg() + theme_classic() + xlab("") + ylab("count proportion") + guides(fill = "none")
b<-ggplot(subset(KKfull, length<=1500), aes(Open,fill=as.factor(Open))) + geom_bar(aes(type), position = position_fill(reverse = TRUE)) + scale_fill_npg() + theme_classic() + xlab("") + ylab("count proportion") + guides(fill = "none")
c<-ggplot(subset(KKfull, length<=1500), aes(SWN,fill=as.factor(SWN))) + geom_bar(aes(type), position = position_fill(reverse = TRUE)) + scale_fill_npg() + theme_classic() + xlab("") + ylab("count proportion") + guides(fill = "none")
d<-ggplot(subset(KKfull, length<=1500), aes(CLF,fill=as.factor(CLF))) + geom_bar(aes(type), position = position_fill(reverse = TRUE)) + scale_fill_npg() + theme_classic() + xlab("") + ylab("count proportion") + guides(fill = "none")
e<-ggplot(subset(KKfull, length<=1500), aes(H2AK121ub1,fill=as.factor(H2AK121ub1))) + geom_bar(aes(type), position = position_fill(reverse = TRUE)) + scale_fill_npg() + theme_classic() + xlab("") + ylab("count proportion") + guides(fill = "none")
f<-ggplot(subset(KKfull, length<=1500), aes(H3K27ac,fill=as.factor(H3K27ac))) + geom_bar(aes(type), position = position_fill(reverse = TRUE)) + scale_fill_npg() + theme_classic() + xlab("") + ylab("count proportion") + guides(fill = "none")
g<-ggplot(subset(KKfull, length<=1500), aes(H3K4me3,fill=as.factor(H3K4me3))) + geom_bar(aes(type), position = position_fill(reverse = TRUE)) + scale_fill_npg() + theme_classic() + xlab("") + ylab("count proportion") + guides(fill = "none")
h<-ggplot(subset(KKfull, length<=1500), aes(H3K36me3,fill=as.factor(H3K36me3))) + geom_bar(aes(type), position = position_fill(reverse = TRUE)) + scale_fill_npg() + theme_classic() + xlab("") + ylab("count proportion") + guides(fill = "none")
i<-ggplot(subset(KKfull, length<=1500), aes(annotation,fill=annotation)) + geom_bar(aes(type), position = position_fill(reverse = TRUE)) + scale_fill_npg() + theme_classic() + xlab("") + ylab("count proportion") + guides(fill = "none")

aaa<-ggplot(subset(kksmall, length<=1500), aes(H3K27me3,fill=as.factor(H3K27me3))) + geom_bar(aes(TRB), position = position_fill(reverse = TRUE)) + scale_fill_npg() + theme_classic() + xlab("") + ylab("count proportion") + guides(fill = "none")
bbb<-ggplot(subset(kksmall, length<=1500), aes(Open,fill=as.factor(Open))) + geom_bar(aes(TRB), position = position_fill(reverse = TRUE)) + scale_fill_npg() + theme_classic() + xlab("") + ylab("count proportion") + guides(fill = "none")
ccc<-ggplot(subset(kksmall, length<=1500), aes(SWN,fill=as.factor(SWN))) + geom_bar(aes(TRB), position = position_fill(reverse = TRUE)) + scale_fill_npg() + theme_classic() + xlab("") + ylab("count proportion") + guides(fill = "none")
ddd<-ggplot(subset(kksmall, length<=1500), aes(CLF,fill=as.factor(CLF))) + geom_bar(aes(TRB), position = position_fill(reverse = TRUE)) + scale_fill_npg() + theme_classic() + xlab("") + ylab("count proportion") + guides(fill = "none")
eee<-ggplot(subset(kksmall, length<=1500), aes(H2AK121ub1,fill=as.factor(H2AK121ub1))) + geom_bar(aes(TRB), position = position_fill(reverse = TRUE)) + scale_fill_npg() + theme_classic() + xlab("") + ylab("count proportion") + guides(fill = "none")
fff<-ggplot(subset(kksmall, length<=1500), aes(H3K27ac,fill=as.factor(H3K27ac))) + geom_bar(aes(TRB), position = position_fill(reverse = TRUE)) + scale_fill_npg() + theme_classic() + xlab("") + ylab("count proportion") + guides(fill = "none")
ggg<-ggplot(subset(kksmall, length<=1500), aes(H3K4me3,fill=as.factor(H3K4me3))) + geom_bar(aes(TRB), position = position_fill(reverse = TRUE)) + scale_fill_npg() + theme_classic() + xlab("") + ylab("count proportion") + guides(fill = "none")
hhh<-ggplot(subset(kksmall, length<=1500), aes(H3K36me3,fill=as.factor(H3K36me3))) + geom_bar(aes(TRB), position = position_fill(reverse = TRUE)) + scale_fill_npg() + theme_classic() + xlab("") + ylab("count proportion") + guides(fill = "none")
iii<-ggplot(subset(kksmall, length<=1500), aes(annotation,fill=annotation)) + geom_bar(aes(TRB), position = position_fill(reverse = TRUE)) + scale_fill_npg() + theme_classic() + xlab("") + ylab("count proportion") + guides(fill = "none")

aa4<-ggplot(subset(kksmall, length<=1500), aes(cluster,fill=as.factor(H3K27me3))) + geom_bar(aes(cluster), position = position_fill(reverse = TRUE)) + scale_fill_npg() + theme_classic() + xlab("") + ylab("count proportion")
bb4<-ggplot(subset(kksmall, length<=1500), aes(Open,fill=as.factor(Open))) + geom_bar(aes(cluster), position = position_fill(reverse = TRUE)) + scale_fill_npg() + theme_classic() + xlab("") + ylab("count proportion")
cc4<-ggplot(subset(kksmall, length<=1500), aes(SWN,fill=as.factor(SWN))) + geom_bar(aes(cluster), position = position_fill(reverse = TRUE)) + scale_fill_npg() + theme_classic() + xlab("") + ylab("count proportion")
dd4<-ggplot(subset(kksmall, length<=1500), aes(CLF,fill=as.factor(CLF))) + geom_bar(aes(cluster), position = position_fill(reverse = TRUE)) + scale_fill_npg() + theme_classic() + xlab("") + ylab("count proportion")
ee4<-ggplot(subset(kksmall, length<=1500), aes(H2AK121ub1,fill=as.factor(H2AK121ub1))) + geom_bar(aes(cluster), position = position_fill(reverse = TRUE)) + scale_fill_npg() + theme_classic() + xlab("") + ylab("count proportion")
ff4<-ggplot(subset(kksmall, length<=1500), aes(H3K27ac,fill=as.factor(H3K27ac))) + geom_bar(aes(cluster), position = position_fill(reverse = TRUE)) + scale_fill_npg() + theme_classic() + xlab("") + ylab("count proportion")
gg4<-ggplot(subset(kksmall, length<=1500), aes(H3K4me3,fill=as.factor(H3K4me3))) + geom_bar(aes(cluster), position = position_fill(reverse = TRUE)) + scale_fill_npg() + theme_classic() + xlab("") + ylab("count proportion")
hh4<-ggplot(subset(kksmall, length<=1500), aes(H3K36me3,fill=as.factor(H3K36me3))) + geom_bar(aes(cluster), position = position_fill(reverse = TRUE)) + scale_fill_npg() + theme_classic() + xlab("") + ylab("count proportion")
ii4<-ggplot(subset(kksmall, length<=1500), aes(annotation,fill=annotation)) + geom_bar(aes(cluster), position = position_fill(reverse = TRUE)) + scale_fill_npg() + theme_classic() + xlab("") + ylab("count proportion")


aa<-ggplot(subset(TRBmelt, type=="peak" & value==1 & length<=1500), aes(variable, fill=as.factor(H3K27me3)))+ geom_bar(aes(variable), position = position_fill(reverse = TRUE)) + scale_fill_npg() + theme_classic() + xlab("") + ylab("count proportion") + guides(fill = "none")
bb<-ggplot(subset(TRBmelt, type=="peak" & value==1 & length<=1500), aes(variable, fill=as.factor(Open)))+ geom_bar(aes(variable), position = position_fill(reverse = TRUE)) + scale_fill_npg() + theme_classic() + xlab("") + ylab("count proportion") + guides(fill = "none")
cc<-ggplot(subset(TRBmelt, type=="peak" & value==1 & length<=1500), aes(variable, fill=as.factor(SWN)))+ geom_bar(aes(variable), position = position_fill(reverse = TRUE)) + scale_fill_npg() + theme_classic() + xlab("") + ylab("count proportion") + guides(fill = "none")
dd<-ggplot(subset(TRBmelt, type=="peak" & value==1 & length<=1500), aes(variable, fill=as.factor(CLF)))+ geom_bar(aes(variable), position = position_fill(reverse = TRUE)) + scale_fill_npg() + theme_classic() + xlab("") + ylab("count proportion") + guides(fill = "none")
ee<-ggplot(subset(TRBmelt, type=="peak" & value==1 & length<=1500), aes(variable, fill=as.factor(H2AK121ub1)))+ geom_bar(aes(variable), position = position_fill(reverse = TRUE)) + scale_fill_npg() + theme_classic() + xlab("") + ylab("count proportion") + guides(fill = "none")
ff<-ggplot(subset(TRBmelt, type=="peak" & value==1 & length<=1500), aes(variable, fill=as.factor(H3K27ac)))+ geom_bar(aes(variable), position = position_fill(reverse = TRUE)) + scale_fill_npg() + theme_classic() + xlab("") + ylab("count proportion") + guides(fill = "none")
gg<-ggplot(subset(TRBmelt, type=="peak" & value==1 & length<=1500), aes(variable, fill=as.factor(H3K4me3)))+ geom_bar(aes(variable), position = position_fill(reverse = TRUE)) + scale_fill_npg() + theme_classic() + xlab("") + ylab("count proportion") + guides(fill = "none")
hh<-ggplot(subset(TRBmelt, type=="peak" & value==1 & length<=1500), aes(variable, fill=as.factor(H3K36me3)))+ geom_bar(aes(variable), position = position_fill(reverse = TRUE)) + scale_fill_npg() + theme_classic() + xlab("") + ylab("count proportion") + guides(fill = "none")
ii<-ggplot(subset(TRBmelt, type=="peak" & value==1 & length<=1500), aes(variable, fill=annotation))+ geom_bar(aes(variable), position = position_fill(reverse = TRUE)) + scale_fill_npg() + theme_classic() + xlab("") + ylab("count proportion") + guides(fill = "none")




pdf(file = "results/ChIP-seq/chromatin_overlap.pdf", width = 12, height = 30)
ggarrange(a,aa,aaa,aa4,b,bb,bbb,bb4,c,cc,ccc,cc4,d,dd,ddd,dd4,e,ee,eee,ee4,f,ff,fff,ff4,g,gg,ggg,gg4,h,hh,hhh,hh4,i,ii,iii,ii4, ncol = 4, nrow = 9)
dev.off()









quit(save="no")


