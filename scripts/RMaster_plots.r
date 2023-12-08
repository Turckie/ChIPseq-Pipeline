#!/usr/bin/Rscript


library("RColorBrewer")
library("ggplot2")
library("ggsci")
library("ggpubr")
library("reshape2")
library("dplyr")
library("readr")

KKfull <- read.table("data/ChIP-seq/RMASTER/Araport11.txt","\t",  header = TRUE, sep = "\t")

#transform motif hit to factors
KKfull$CRACCTA <- as.numeric(KKfull$CGACCTA) + as.numeric(KKfull$CAACCTA)
KKfull$CRACCTA[as.numeric(KKfull$CRACCTA) >=6 ]<- "more than 5"
KKfull$AAACCCTA[as.numeric(KKfull$AAACCCTA) >=6 ]<- "more than 5"
KKfull$CGACCTA[as.numeric(KKfull$CGACCTA) >=6 ]<- "more than 5"
KKfull$CAACCTA[as.numeric(KKfull$CAACCTA) >=6 ]<- "more than 5"
KKfull$DWTTAGGKTKT[as.numeric(KKfull$DWTTAGGKTKT) >=6 ]<- "more than 5"
KKfull$MA1353.1[as.numeric(KKfull$MA1353.1) >=6 ]<- "more than 5"
KKfull$MA1352.1[as.numeric(KKfull$MA1352.1) >=6] <- "more than 5"
#KKfull$MA1354.1[as.numeric(KKfull$MA1354.1) >=6 ]<- "more than 5"
KKfull$MA1404.1[as.numeric(KKfull$MA1404.1) >=6 ]<- "more than 5"
KKfull$AGARRAAGAARRAGA[as.numeric(KKfull$AGARRAAGAARRAGA) >=6 ]<- "more than 5"


#MA1277.1 MA1402.1 MA1268.1 MA1416.1 MA1403.1 AAAAAAAAAAAAAAA 2-ARGCCCAWT

KKfull$length<-as.numeric(KKfull$stop)-as.numeric(KKfull$start)
KKfull$category <-factor(KKfull$category, levels=c( "TRB1", "TRB2", "TRB3","TRB1_unique","TRB1_2_private", "TRB123_common", "TRB2_unique","TRB2_3_private","TRB3_unique", "TRB1_3_private"))


##########################################

pdf(file = "results/TRB_fragment_length.pdf", width = 4, height = 12)
ggplot(KKfull,aes(length, color=type))+geom_histogram( )+ facet_grid(category ~ ., scales = "free_y")
dev.off()
##########################################

#overlap TRBs with chromatin features
a<-ggplot(subset(KKTRB, length<=1500), aes(category,fill=as.factor(H3K27me3))) + geom_bar(aes(type), position = position_fill(reverse = TRUE)) + scale_fill_npg() + facet_wrap(~ TRB, nrow=1) + theme_classic() + xlab("") + ylab("count proportion") + theme(legend.position="none")
b<-ggplot(subset(KKTRB, length<=1500), aes(category,fill=as.factor(CLFSWN))) + geom_bar(aes(type), position = position_fill(reverse = TRUE)) + scale_fill_npg() + facet_wrap(~ TRB, nrow=1) + theme_classic() + xlab("") + ylab("count proportion") + theme(legend.position="none")
c<-ggplot(subset(KKTRB, length<=1500), aes(category,fill=as.factor(H2AK121ub1))) + geom_bar(aes(type), position = position_fill(reverse = TRUE)) + scale_fill_npg() + facet_wrap(~ TRB, nrow=1) + theme_classic() + xlab("") + ylab("count proportion") + theme(legend.position="none")
d<-ggplot(subset(KKTRB, length<=1500), aes(category,fill=as.factor(Open))) + geom_bar(aes(type), position = position_fill(reverse = TRUE)) + scale_fill_npg() + facet_wrap(~ TRB, nrow=1) + theme_classic() + xlab("") + ylab("count proportion") + theme(legend.position="none")
e<-ggplot(subset(KKTRB, length<=1500), aes(category,fill=as.factor(H3K27ac))) + geom_bar(aes(type), position = position_fill(reverse = TRUE)) + scale_fill_npg() + facet_wrap(~ TRB, nrow=1) + theme_classic() + xlab("") + ylab("count proportion") + theme(legend.position="none")
f<-ggplot(subset(KKTRB, length<=1500), aes(category,fill=as.factor(H3K4me3))) + geom_bar(aes(type), position = position_fill(reverse = TRUE)) + scale_fill_npg() + facet_wrap(~ TRB, nrow=1) + theme_classic() + xlab("") + ylab("count proportion") + theme(legend.position="none")
g<-ggplot(subset(KKTRB, length<=1500), aes(category,fill=as.factor(H3K36me3))) + geom_bar(aes(type), position = position_fill(reverse = TRUE)) + scale_fill_npg() + facet_wrap(~ TRB, nrow=1) + theme_classic() + xlab("") + ylab("count proportion") + theme(legend.position="none")
h<-ggplot(subset(KKTRB, length<=1500), aes(category,fill=annotation)) + geom_bar(aes(type), position = position_fill(reverse = TRUE)) + scale_fill_npg() + facet_wrap(~ TRB, nrow=1) + theme_classic() + xlab("") + ylab("count proportion") + theme(legend.position="none")

aa<-ggplot(subset(KKfull, category %in% c("TRB1_unique", "TRB2_3_private") & length<=1500), aes(category,fill=as.factor(H3K27me3))) + geom_bar(aes(type), position = position_fill(reverse = TRUE)) + scale_fill_npg() + facet_wrap(~ category, nrow=1) + theme_classic() + xlab("") + ylab("count proportion") 
bb<-ggplot(subset(KKfull, category %in% c("TRB1_unique", "TRB2_3_private") & length<=1500), aes(category,fill=as.factor(CLFSWN))) + geom_bar(aes(type), position = position_fill(reverse = TRUE)) + scale_fill_npg() + facet_wrap(~ category, nrow=1) + theme_classic() + xlab("") + ylab("count proportion") 
cc<-ggplot(subset(KKfull, category %in% c("TRB1_unique", "TRB2_3_private") & length<=1500), aes(category,fill=as.factor(H2AK121ub1))) + geom_bar(aes(type), position = position_fill(reverse = TRUE)) + scale_fill_npg() + facet_wrap(~ category, nrow=1) + theme_classic() + xlab("") + ylab("count proportion") 
dd<-ggplot(subset(KKfull, category %in% c("TRB1_unique", "TRB2_3_private") & length<=1500), aes(category,fill=as.factor(Open))) + geom_bar(aes(type), position = position_fill(reverse = TRUE)) + scale_fill_npg() + facet_wrap(~ category, nrow=1) + theme_classic() + xlab("") + ylab("count proportion") 
ee<-ggplot(subset(KKfull, category %in% c("TRB1_unique", "TRB2_3_private") & length<=1500), aes(category,fill=as.factor(H3K27ac))) + geom_bar(aes(type), position = position_fill(reverse = TRUE)) + scale_fill_npg() + facet_wrap(~ category, nrow=1) + theme_classic() + xlab("") + ylab("count proportion") 
ff<-ggplot(subset(KKfull, category %in% c("TRB1_unique", "TRB2_3_private") & length<=1500), aes(category,fill=as.factor(H3K4me3))) + geom_bar(aes(type), position = position_fill(reverse = TRUE)) + scale_fill_npg() + facet_wrap(~ category, nrow=1) + theme_classic() + xlab("") + ylab("count proportion") 
gg<-ggplot(subset(KKfull, category %in% c("TRB1_unique", "TRB2_3_private") & length<=1500), aes(category,fill=as.factor(H3K36me3))) + geom_bar(aes(type), position = position_fill(reverse = TRUE)) + scale_fill_npg() + facet_wrap(~ category, nrow=1) + theme_classic() + xlab("") + ylab("count proportion") 
hh<-ggplot(subset(KKfull, category %in% c("TRB1_unique", "TRB2_3_private") & length<=1500), aes(category,fill=annotation)) + geom_bar(aes(type), position = position_fill(reverse = TRUE)) + scale_fill_npg() + facet_wrap(~ category, nrow=1) + theme_classic() + xlab("") + ylab("count proportion")

aaa<-ggplot(subset(KKfull,length<=1500), aes(category,fill=as.factor(H3K27me3))) + geom_bar(aes(type), position = position_fill(reverse = TRUE)) + scale_fill_npg() + facet_wrap(~ category, nrow=1) + theme_classic() + xlab("") + ylab("count proportion") 
bbb<-ggplot(subset(KKfull, length<=1500), aes(category,fill=as.factor(CLFSWN))) + geom_bar(aes(type), position = position_fill(reverse = TRUE)) + scale_fill_npg() + facet_wrap(~ category, nrow=1) + theme_classic() + xlab("") + ylab("count proportion") 
ccc<-ggplot(subset(KKfull, length<=1500), aes(category,fill=as.factor(H2AK121ub1))) + geom_bar(aes(type), position = position_fill(reverse = TRUE)) + scale_fill_npg() + facet_wrap(~ category, nrow=1) + theme_classic() + xlab("") + ylab("count proportion") 
ddd<-ggplot(subset(KKfull,   length<=1500), aes(category,fill=as.factor(Open))) + geom_bar(aes(type), position = position_fill(reverse = TRUE)) + scale_fill_npg() + facet_wrap(~ category, nrow=1) + theme_classic() + xlab("") + ylab("count proportion") 
eee<-ggplot(subset(KKfull,length<=1500), aes(category,fill=as.factor(H3K27ac))) + geom_bar(aes(type), position = position_fill(reverse = TRUE)) + scale_fill_npg() + facet_wrap(~ category, nrow=1) + theme_classic() + xlab("") + ylab("count proportion") 
fff<-ggplot(subset(KKfull, length<=1500), aes(category,fill=as.factor(H3K4me3))) + geom_bar(aes(type), position = position_fill(reverse = TRUE)) + scale_fill_npg() + facet_wrap(~ category, nrow=1) + theme_classic() + xlab("") + ylab("count proportion") 
ggg<-ggplot(subset(KKfull, length<=1500), aes(category,fill=as.factor(H3K36me3))) + geom_bar(aes(type), position = position_fill(reverse = TRUE)) + scale_fill_npg() + facet_wrap(~ category, nrow=1) + theme_classic() + xlab("") + ylab("count proportion") 
hhh<-ggplot(subset(KKfull, length<=1500), aes(category,fill=annotation)) + geom_bar(aes(type), position = position_fill(reverse = TRUE)) + scale_fill_npg() + facet_wrap(~ category, nrow=1) + theme_classic() + xlab("") + ylab("count proportion")

aaaa<-ggplot(subset(KKfull, category %in% c("TRB1_unique","TRB1_2_private", "TRB123_common", "TRB2_unique","TRB2_3_private","TRB3_unique", "TRB1_3_private") & length<=1500), aes(category,fill=as.factor(H3K27me3))) + geom_bar(aes(type), position = position_fill(reverse = TRUE)) + scale_fill_npg() + facet_wrap(~ category, nrow=1) + theme_classic() + xlab("") + ylab("count proportion") 
bbbb<-ggplot(subset(KKfull, category %in% c("TRB1_unique","TRB1_2_private", "TRB123_common", "TRB2_unique","TRB2_3_private","TRB3_unique", "TRB1_3_private") & length<=1500), aes(category,fill=as.factor(CLFSWN))) + geom_bar(aes(type), position = position_fill(reverse = TRUE)) + scale_fill_npg() + facet_wrap(~ category, nrow=1) + theme_classic() + xlab("") + ylab("count proportion") 
cccc<-ggplot(subset(KKfull, category %in% c("TRB1_unique","TRB1_2_private", "TRB123_common", "TRB2_unique","TRB2_3_private","TRB3_unique", "TRB1_3_private") & length<=1500), aes(category,fill=as.factor(H2AK121ub1))) + geom_bar(aes(type), position = position_fill(reverse = TRUE)) + scale_fill_npg() + facet_wrap(~ category, nrow=1) + theme_classic() + xlab("") + ylab("count proportion") 
dddd<-ggplot(subset(KKfull, category %in% c("TRB1_unique","TRB1_2_private", "TRB123_common", "TRB2_unique","TRB2_3_private","TRB3_unique", "TRB1_3_private") & length<=1500), aes(category,fill=as.factor(Open))) + geom_bar(aes(type), position = position_fill(reverse = TRUE)) + scale_fill_npg() + facet_wrap(~ category, nrow=1) + theme_classic() + xlab("") + ylab("count proportion") 
eeee<-ggplot(subset(KKfull, category %in% c("TRB1_unique","TRB1_2_private", "TRB123_common", "TRB2_unique","TRB2_3_private","TRB3_unique", "TRB1_3_private") & length<=1500), aes(category,fill=as.factor(H3K27ac))) + geom_bar(aes(type), position = position_fill(reverse = TRUE)) + scale_fill_npg() + facet_wrap(~ category, nrow=1) + theme_classic() + xlab("") + ylab("count proportion") 
ffff<-ggplot(subset(KKfull, category %in% c("TRB1_unique","TRB1_2_private", "TRB123_common", "TRB2_unique","TRB2_3_private","TRB3_unique", "TRB1_3_private") & length<=1500), aes(category,fill=as.factor(H3K4me3))) + geom_bar(aes(type), position = position_fill(reverse = TRUE)) + scale_fill_npg() + facet_wrap(~ category, nrow=1) + theme_classic() + xlab("") + ylab("count proportion") 
gggg<-ggplot(subset(KKfull, category %in% c("TRB1_unique","TRB1_2_private", "TRB123_common", "TRB2_unique","TRB2_3_private","TRB3_unique", "TRB1_3_private") & length<=1500), aes(category,fill=as.factor(H3K36me3))) + geom_bar(aes(type), position = position_fill(reverse = TRUE)) + scale_fill_npg() + facet_wrap(~ category, nrow=1) + theme_classic() + xlab("") + ylab("count proportion") 
hhhh<-ggplot(subset(KKfull, category %in% c("TRB1_unique","TRB1_2_private", "TRB123_common", "TRB2_unique","TRB2_3_private","TRB3_unique", "TRB1_3_private") & length<=1500), aes(category,fill=annotation)) + geom_bar(aes(type), position = position_fill(reverse = TRUE)) + scale_fill_npg() + facet_wrap(~ category, nrow=1) + theme_classic() + xlab("") + ylab("count proportion")



pdf(file = "results/ChIP-seq/chromatin_overlap_selection.pdf", width = 12, height = 20)
ggarrange(a,aa,b,bb,c,cc,d,dd,e,ee,f,ff,g,gg,h,hh, ncol = 2, nrow = 9)
dev.off()
pdf(file = "results/ChIP-seq/chromatin_overlap_TRBplain.pdf", width = 3, height = 20)
ggarrange(a,b,c,d,e,f,g,h, ncol = 1, nrow = 8)
dev.off()

pdf(file = "results/ChIP-seq/chromatin_overlapall.pdf", width = 12, height = 20)
ggarrange(aaa,bbb,ccc,ddd,eee,fff,ggg,hhh, ncol = 1, nrow = 8)
dev.off()

pdf(file = "results/ChIP-seq/chromatin_overlapallsubtypes.pdf", width = 12, height = 20)
ggarrange(aaaa,bbbb,cccc,dddd,eeee,ffff,gggg,hhhh, ncol = 1, nrow = 8)
dev.off()



#####################cismotifs

ma<-ggplot(subset(KKfull, category %in% c("TRB1", "TRB2", "TRB3") & length<=1500), aes(category,fill=as.factor(DWTTAGGKTKT))) + geom_bar(aes(type), position = position_fill(reverse = TRUE)) + scale_fill_npg() + facet_wrap(~ category, nrow=1) + theme_classic() + xlab("") + ylab("count proportion") #+ theme(legend.position="none")
mb<-ggplot(subset(KKfull, category %in% c("TRB1", "TRB2", "TRB3") &   length<=1500), aes(category,fill=as.factor(AAACCCTA))) + geom_bar(aes(type), position = position_fill(reverse = TRUE)) + scale_fill_npg() + facet_wrap(~ category, nrow=1) + theme_classic() + xlab("") + ylab("count proportion") #+ theme(legend.position="none")
mc<-ggplot(subset(KKfull, category %in% c("TRB1", "TRB2", "TRB3") &   length<=1500), aes(category,fill=as.factor(CRACCTA))) + geom_bar(aes(type), position = position_fill(reverse = TRUE)) + scale_fill_npg() + facet_wrap(~ category, nrow=1) + theme_classic() + xlab("") + ylab("count proportion") #+ theme(legend.position="none")

mma<-ggplot(subset(KKfull, category %in% c("TRB1_unique","TRB1_2_private", "TRB123_common", "TRB2_unique","TRB2_3_private","TRB3_unique", "TRB1_3_private") & length<=1500), aes(category,fill=as.factor(DWTTAGGKTKT))) + geom_bar(aes(type), position = position_fill(reverse = TRUE)) + scale_fill_npg() + facet_wrap(~ category, nrow=1) + theme_classic() + xlab("") + ylab("count proportion") #+ theme(legend.position="none")
mmb<-ggplot(subset(KKfull, category %in% c("TRB1_unique","TRB1_2_private", "TRB123_common", "TRB2_unique","TRB2_3_private","TRB3_unique", "TRB1_3_private") &  length<=1500), aes(category,fill=as.factor(AAACCCTA))) + geom_bar(aes(type), position = position_fill(reverse = TRUE)) + scale_fill_npg() + facet_wrap(~ category, nrow=1) + theme_classic() + xlab("") + ylab("count proportion") #+ theme(legend.position="none")
mmc<-ggplot(subset(KKfull, category %in% c("TRB1_unique","TRB1_2_private", "TRB123_common", "TRB2_unique","TRB2_3_private","TRB3_unique", "TRB1_3_private") &  length<=1500), aes(category,fill=as.factor(CRACCTA))) + geom_bar(aes(type), position = position_fill(reverse = TRUE)) + scale_fill_npg() + facet_wrap(~ category, nrow=1) + theme_classic() + xlab("") + ylab("count proportion") #+ theme(legend.position="none")


pdf(file = "results/ChIP-seq/motif_overlap_TRBplain.pdf", width = 6, height = 6)
ggarrange(ma,mb,mc, ncol = 1, nrow = 3)
dev.off()

pdf(file = "results/ChIP-seq/motif_overlap_TRBspecifc.pdf", width = 6, height = 6)
ggarrange(mma,mmb,mmc, ncol = 1, nrow = 3)
dev.off()


quit(save="no")

