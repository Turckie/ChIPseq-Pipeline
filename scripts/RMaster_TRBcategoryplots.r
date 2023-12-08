#!/usr/bin/Rscript


library("RColorBrewer")
library("ggplot2")
library("ggsci")
library("ggpubr")
library("reshape2")
library("dplyr")
library("readr")

KKTRB1 <- read.table("data/ChIP-seq/RMASTER/Araport11.TRB1_controls.txt","\t",  header = TRUE, sep = "\t")
KKTRB2 <- read.table("data/ChIP-seq/RMASTER/Araport11.TRB2_controls.txt","\t",  header = TRUE, sep = "\t")
KKTRB3 <- read.table("data/ChIP-seq/RMASTER/Araport11.TRB3_controls.txt","\t",  header = TRUE, sep = "\t")
KKTRB1_unique <- read.table("data/ChIP-seq/RMASTER/Araport11.TRB1_unique_controls.txt","\t",  header = TRUE, sep = "\t")
KKTRB2_3 <- read.table("data/ChIP-seq/RMASTER/Araport11.TRB2_3_private_controls.txt","\t",  header = TRUE, sep = "\t")

KKTRB1$category <- "TRB1"
KKTRB2$category <- "TRB2"
KKTRB3$category <- "TRB3"


KKTRB1$length<-as.numeric(KKTRB1$stop)-as.numeric(KKTRB1$start)
KKTRB2$length<-as.numeric(KKTRB2$stop)-as.numeric(KKTRB2$start)
KKTRB3$length<-as.numeric(KKTRB3$stop)-as.numeric(KKTRB3$start)
KKTRB1_unique$length<-as.numeric(KKTRB1_unique$stop)-as.numeric(KKTRB1_unique$start)
KKTRB2_3$length<-as.numeric(KKTRB2_3$stop)-as.numeric(KKTRB2_3$start)


TRBmelt = melt(KKfull, id.vars = c("type", "H3K27me3", "Open", "SWN", "CLF", "H2AK121ub1","H3K27ac", "H3K4me3", "H3K36me3", "cluster", "annotation","length"), measure.vars = c("TRB1", "TRB2", "TRB3"))

##########################################

pdf(file = "results/TRB_fragment_length.pdf", width = 4, height = 4)
ggplot(KKTRB1,aes(length, color=type))+geom_histogram( )
ggplot(KKTRB2,aes(length, color=type))+geom_histogram( )
ggplot(KKTRB3,aes(length, color=type))+geom_histogram( )
dev.off()
##########################################

#overlap TRBs with chromatin features
a<-ggplot(subset(KKTRB1, length<=1500), aes(H3K27me3,fill=as.factor(H3K27me3))) + geom_bar(aes(type), position = position_fill(reverse = TRUE)) + scale_fill_npg() + theme_classic() + xlab("") + ylab("count proportion") + guides(fill = "none")
b<-ggplot(subset(KKTRB1, length<=1500), aes(Open,fill=as.factor(Open))) + geom_bar(aes(type), position = position_fill(reverse = TRUE)) + scale_fill_npg() + theme_classic() + xlab("") + ylab("count proportion") + guides(fill = "none")
c<-ggplot(subset(KKTRB1, length<=1500), aes(SWN,fill=as.factor(SWN))) + geom_bar(aes(type), position = position_fill(reverse = TRUE)) + scale_fill_npg() + theme_classic() + xlab("") + ylab("count proportion") + guides(fill = "none")
d<-ggplot(subset(KKTRB1, length<=1500), aes(CLF,fill=as.factor(CLF))) + geom_bar(aes(type), position = position_fill(reverse = TRUE)) + scale_fill_npg() + theme_classic() + xlab("") + ylab("count proportion") + guides(fill = "none")
e<-ggplot(subset(KKTRB1, length<=1500), aes(H2AK121ub1,fill=as.factor(H2AK121ub1))) + geom_bar(aes(type), position = position_fill(reverse = TRUE)) + scale_fill_npg() + theme_classic() + xlab("") + ylab("count proportion") + guides(fill = "none")
f<-ggplot(subset(KKTRB1, length<=1500), aes(H3K27ac,fill=as.factor(H3K27ac))) + geom_bar(aes(type), position = position_fill(reverse = TRUE)) + scale_fill_npg() + theme_classic() + xlab("") + ylab("count proportion") + guides(fill = "none")
g<-ggplot(subset(KKTRB1, length<=1500), aes(H3K4me3,fill=as.factor(H3K4me3))) + geom_bar(aes(type), position = position_fill(reverse = TRUE)) + scale_fill_npg() + theme_classic() + xlab("") + ylab("count proportion") + guides(fill = "none")
h<-ggplot(subset(KKTRB1, length<=1500), aes(H3K36me3,fill=as.factor(H3K36me3))) + geom_bar(aes(type), position = position_fill(reverse = TRUE)) + scale_fill_npg() + theme_classic() + xlab("") + ylab("count proportion") + guides(fill = "none")
i<-ggplot(subset(KKTRB1, length<=1500), aes(annotation,fill=annotation)) + geom_bar(aes(type), position = position_fill(reverse = TRUE)) + scale_fill_npg() + theme_classic() + xlab("") + ylab("count proportion") + guides(fill = "none")

#overlap TRBs with chromatin features
aa<-ggplot(subset(KKTRB2, length<=1500), aes(H3K27me3,fill=as.factor(H3K27me3))) + geom_bar(aes(type), position = position_fill(reverse = TRUE)) + scale_fill_npg() + theme_classic() + xlab("") + ylab("count proportion") + guides(fill = "none")
bb<-ggplot(subset(KKTRB2, length<=1500), aes(Open,fill=as.factor(Open))) + geom_bar(aes(type), position = position_fill(reverse = TRUE)) + scale_fill_npg() + theme_classic() + xlab("") + ylab("count proportion") + guides(fill = "none")
cc<-ggplot(subset(KKTRB2, length<=1500), aes(SWN,fill=as.factor(SWN))) + geom_bar(aes(type), position = position_fill(reverse = TRUE)) + scale_fill_npg() + theme_classic() + xlab("") + ylab("count proportion") + guides(fill = "none")
dd<-ggplot(subset(KKTRB2, length<=1500), aes(CLF,fill=as.factor(CLF))) + geom_bar(aes(type), position = position_fill(reverse = TRUE)) + scale_fill_npg() + theme_classic() + xlab("") + ylab("count proportion") + guides(fill = "none")
ee<-ggplot(subset(KKTRB2, length<=1500), aes(H2AK121ub1,fill=as.factor(H2AK121ub1))) + geom_bar(aes(type), position = position_fill(reverse = TRUE)) + scale_fill_npg() + theme_classic() + xlab("") + ylab("count proportion") + guides(fill = "none")
ff<-ggplot(subset(KKTRB2, length<=1500), aes(H3K27ac,fill=as.factor(H3K27ac))) + geom_bar(aes(type), position = position_fill(reverse = TRUE)) + scale_fill_npg() + theme_classic() + xlab("") + ylab("count proportion") + guides(fill = "none")
gg<-ggplot(subset(KKTRB2, length<=1500), aes(H3K4me3,fill=as.factor(H3K4me3))) + geom_bar(aes(type), position = position_fill(reverse = TRUE)) + scale_fill_npg() + theme_classic() + xlab("") + ylab("count proportion") + guides(fill = "none")
hh<-ggplot(subset(KKTRB2, length<=1500), aes(H3K36me3,fill=as.factor(H3K36me3))) + geom_bar(aes(type), position = position_fill(reverse = TRUE)) + scale_fill_npg() + theme_classic() + xlab("") + ylab("count proportion") + guides(fill = "none")
ii<-ggplot(subset(KKTRB2, length<=1500), aes(annotation,fill=annotation)) + geom_bar(aes(type), position = position_fill(reverse = TRUE)) + scale_fill_npg() + theme_classic() + xlab("") + ylab("count proportion") + guides(fill = "none")



#overlap TRBs with chromatin features
aaa<-ggplot(subset(KKTRB3, length<=1500), aes(H3K27me3,fill=as.factor(H3K27me3))) + geom_bar(aes(type), position = position_fill(reverse = TRUE)) + scale_fill_npg() + theme_classic() + xlab("") + ylab("count proportion")
bbb<-ggplot(subset(KKTRB3, length<=1500), aes(Open,fill=as.factor(Open))) + geom_bar(aes(type), position = position_fill(reverse = TRUE)) + scale_fill_npg() + theme_classic() + xlab("") + ylab("count proportion")
ccc<-ggplot(subset(KKTRB3, length<=1500), aes(SWN,fill=as.factor(SWN))) + geom_bar(aes(type), position = position_fill(reverse = TRUE)) + scale_fill_npg() + theme_classic() + xlab("") + ylab("count proportion")
ddd<-ggplot(subset(KKTRB3, length<=1500), aes(CLF,fill=as.factor(CLF))) + geom_bar(aes(type), position = position_fill(reverse = TRUE)) + scale_fill_npg() + theme_classic() + xlab("") + ylab("count proportion")
eee<-ggplot(subset(KKTRB3, length<=1500), aes(H2AK121ub1,fill=as.factor(H2AK121ub1))) + geom_bar(aes(type), position = position_fill(reverse = TRUE)) + scale_fill_npg() + theme_classic() + xlab("") + ylab("count proportion")
fff<-ggplot(subset(KKTRB3, length<=1500), aes(H3K27ac,fill=as.factor(H3K27ac))) + geom_bar(aes(type), position = position_fill(reverse = TRUE)) + scale_fill_npg() + theme_classic() + xlab("") + ylab("count proportion")
ggg<-ggplot(subset(KKTRB3, length<=1500), aes(H3K4me3,fill=as.factor(H3K4me3))) + geom_bar(aes(type), position = position_fill(reverse = TRUE)) + scale_fill_npg() + theme_classic() + xlab("") + ylab("count proportion")
hhh<-ggplot(subset(KKTRB3, length<=1500), aes(H3K36me3,fill=as.factor(H3K36me3))) + geom_bar(aes(type), position = position_fill(reverse = TRUE)) + scale_fill_npg() + theme_classic() + xlab("") + ylab("count proportion")
iii<-ggplot(subset(KKTRB3, length<=1500), aes(annotation,fill=annotation)) + geom_bar(aes(type), position = position_fill(reverse = TRUE)) + scale_fill_npg() + theme_classic() + xlab("") + ylab("count proportion")


a1u<-ggplot(subset(KKTRB1_unique, length<=1500), aes(H3K27me3,fill=as.factor(H3K27me3))) + geom_bar(aes(type), position = position_fill(reverse = TRUE)) + scale_fill_npg() + theme_classic() + xlab("") + ylab("count proportion") + guides(fill = "none")
b1u<-ggplot(subset(KKTRB1_unique, length<=1500), aes(Open,fill=as.factor(Open))) + geom_bar(aes(type), position = position_fill(reverse = TRUE)) + scale_fill_npg() + theme_classic() + xlab("") + ylab("count proportion") + guides(fill = "none")
c1u<-ggplot(subset(KKTRB1_unique, length<=1500), aes(SWN,fill=as.factor(SWN))) + geom_bar(aes(type), position = position_fill(reverse = TRUE)) + scale_fill_npg() + theme_classic() + xlab("") + ylab("count proportion") + guides(fill = "none")
d1u<-ggplot(subset(KKTRB1_unique, length<=1500), aes(CLF,fill=as.factor(CLF))) + geom_bar(aes(type), position = position_fill(reverse = TRUE)) + scale_fill_npg() + theme_classic() + xlab("") + ylab("count proportion") + guides(fill = "none")
e1u<-ggplot(subset(KKTRB1_unique, length<=1500), aes(H2AK121ub1,fill=as.factor(H2AK121ub1))) + geom_bar(aes(type), position = position_fill(reverse = TRUE)) + scale_fill_npg() + theme_classic() + xlab("") + ylab("count proportion") + guides(fill = "none")
f1u<-ggplot(subset(KKTRB1_unique, length<=1500), aes(H3K27ac,fill=as.factor(H3K27ac))) + geom_bar(aes(type), position = position_fill(reverse = TRUE)) + scale_fill_npg() + theme_classic() + xlab("") + ylab("count proportion") + guides(fill = "none")
g1u<-ggplot(subset(KKTRB1_unique, length<=1500), aes(H3K4me3,fill=as.factor(H3K4me3))) + geom_bar(aes(type), position = position_fill(reverse = TRUE)) + scale_fill_npg() + theme_classic() + xlab("") + ylab("count proportion") + guides(fill = "none")
h1u<-ggplot(subset(KKTRB1_unique, length<=1500), aes(H3K36me3,fill=as.factor(H3K36me3))) + geom_bar(aes(type), position = position_fill(reverse = TRUE)) + scale_fill_npg() + theme_classic() + xlab("") + ylab("count proportion") + guides(fill = "none")
i1u<-ggplot(subset(KKTRB1_unique, length<=1500), aes(annotation,fill=annotation)) + geom_bar(aes(type), position = position_fill(reverse = TRUE)) + scale_fill_npg() + theme_classic() + xlab("") + ylab("count proportion") + guides(fill = "none")

a23<-ggplot(subset(KKTRB2_3, length<=1500), aes(H3K27me3,fill=as.factor(H3K27me3))) + geom_bar(aes(type), position = position_fill(reverse = TRUE)) + scale_fill_npg() + theme_classic() + xlab("") + ylab("count proportion")
b23<-ggplot(subset(KKTRB2_3, length<=1500), aes(Open,fill=as.factor(Open))) + geom_bar(aes(type), position = position_fill(reverse = TRUE)) + scale_fill_npg() + theme_classic() + xlab("") + ylab("count proportion")
c23<-ggplot(subset(KKTRB2_3, length<=1500), aes(SWN,fill=as.factor(SWN))) + geom_bar(aes(type), position = position_fill(reverse = TRUE)) + scale_fill_npg() + theme_classic() + xlab("") + ylab("count proportion")
d23<-ggplot(subset(KKTRB2_3, length<=1500), aes(CLF,fill=as.factor(CLF))) + geom_bar(aes(type), position = position_fill(reverse = TRUE)) + scale_fill_npg() + theme_classic() + xlab("") + ylab("count proportion")
e23<-ggplot(subset(KKTRB2_3, length<=1500), aes(H2AK121ub1,fill=as.factor(H2AK121ub1))) + geom_bar(aes(type), position = position_fill(reverse = TRUE)) + scale_fill_npg() + theme_classic() + xlab("") + ylab("count proportion")
f23<-ggplot(subset(KKTRB2_3, length<=1500), aes(H3K27ac,fill=as.factor(H3K27ac))) + geom_bar(aes(type), position = position_fill(reverse = TRUE)) + scale_fill_npg() + theme_classic() + xlab("") + ylab("count proportion")
g23<-ggplot(subset(KKTRB2_3, length<=1500), aes(H3K4me3,fill=as.factor(H3K4me3))) + geom_bar(aes(type), position = position_fill(reverse = TRUE)) + scale_fill_npg() + theme_classic() + xlab("") + ylab("count proportion")
h23<-ggplot(subset(KKTRB2_3, length<=1500), aes(H3K36me3,fill=as.factor(H3K36me3))) + geom_bar(aes(type), position = position_fill(reverse = TRUE)) + scale_fill_npg() + theme_classic() + xlab("") + ylab("count proportion")
i23<-ggplot(subset(KKTRB2_3, length<=1500), aes(annotation,fill=annotation)) + geom_bar(aes(type), position = position_fill(reverse = TRUE)) + scale_fill_npg() + theme_classic() + xlab("") + ylab("count proportion")


pdf(file = "results/ChIP-seq/chromatin_overlap.TRBs.pdf", width = 18, height = 30)
ggarrange(a,aa,aaa,b,bb,bbb,c,cc,ccc,d,dd,ddd,e,ee,eee,f,ff,fff,g,gg,ggg,h,hh,hhh,i,ii,iii, ncol = 3, nrow = 9)
dev.off()

pdf(file = "results/ChIP-seq/chromatin_overlap.TRBcats.pdf", width = 18, height = 30)
ggarrange(a1u,a23,b1u,b23,c1u,c23,d1u,d23,e1u,e23,f1u,f23,g1u,g23,h1u,h23,i1u,i23, ncol = 2, nrow = 9)
dev.off()
quit(save="no")


