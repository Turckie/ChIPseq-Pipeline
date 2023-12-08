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
KKfull$CRACCTA <-as.numeric(KKfull$CAACCTAA) + as.numeric(KKfull$CGACCTAA)
KKfull$CRACCTA[as.numeric(KKfull$CRACCTA) >=6 ]<- "more than 5"
KKfull$AAACCCTA[as.numeric(KKfull$AAACCCTA) >=6 ]<- "more than 5"
KKfull$CGACCTAA[as.numeric(KKfull$CGACCTAA) >=6 ]<- "more than 5"
KKfull$CAACCTAA[as.numeric(KKfull$CAACCTAA) >=6 ]<- "more than 5"
KKfull$DWTTAGGKTKT[as.numeric(KKfull$DWTTAGGKTKT) >=6 ]<- "more than 5"
KKfull$MA1353.1[as.numeric(KKfull$MA1353.1) >=6 ]<- "more than 5"
KKfull$MA1352.1[as.numeric(KKfull$MA1352.1) >=6] <- "more than 5"
#KKfull$MA1354.1[as.numeric(KKfull$MA1354.1) >=6 ]<- "more than 5"
KKfull$MA1404.1[as.numeric(KKfull$MA1404.1) >=6 ]<- "more than 5"
KKfull$AGARRAAGAARRAGA[as.numeric(KKfull$AGARRAAGAARRAGA) >=6 ]<- "more than 5"


#MA1277.1 MA1402.1 MA1268.1 MA1416.1 MA1403.1 AAAAAAAAAAAAAAA 2-ARGCCCAWT

KKfull$length<-as.numeric(KKfull$stop)-as.numeric(KKfull$start)
KKfull$category <-factor(KKfull$category, levels=c( "TRB1_unique","TRB1_2_private", "TRB123_common", "TRB2_unique","TRB2_3_private","TRB3_unique", "TRB1_3_private"))


#add a helper column for plain TRB peaks
KKTRB1<-subset(KKfull, category %in% c("TRB123_common", "TRB1_unique","TRB1_2_private","TRB1_3_private"))
KKTRB1$TRB <-"TRB1"
KKTRB2<-subset(KKfull, category %in% c("TRB123_common", "TRB2_unique","TRB1_2_private","TRB2_3_private"))
KKTRB2$TRB <-"TRB2"
KKTRB3<-subset(KKfull, category %in% c("TRB123_common", "TRB3_unique","TRB1_3_private","TRB2_3_private"))
KKTRB3$TRB <-"TRB3"

KKTRB<-rbind(KKTRB1,KKTRB2,KKTRB3)

telo_noCRACTA <-KKTRB %>% group_by (TRB,type) %>% summarize (count = sum (AAACCCTA != '0' & CRACCTA =='0'))
telo_noCRACTA$pie <- "telobox_noCRACCTA"
telo_CRACTA <- KKTRB %>% group_by (TRB,type) %>% summarize (count = sum (AAACCCTA != '0' & CRACCTA !='0'))
telo_CRACTA$pie <- "telobox_CRACCTA"
notelo_noCRACTA <-KKTRB %>% group_by (TRB,type) %>% summarize (count = sum (AAACCCTA =='0' & CRACCTA =='0'))
notelo_noCRACTA$pie <- "notelobox_noCRACCTA"
notelo_CRACTA <- KKTRB %>% group_by (TRB,type) %>% summarize (count = sum (AAACCCTA =='0' & CRACCTA !='0'))
notelo_CRACTA$pie <- "notelobox_CRACCTA"
pie<-rbind(telo_noCRACTA,telo_CRACTA,notelo_noCRACTA, notelo_CRACTA)

p1<- ggplot(subset(pie, type=="peak" & TRB=="TRB1"), aes(x = "TRB", y = count, fill = pie)) +
  geom_bar(stat = "identity", width = 3) +
  geom_text(aes(label = count), position = position_stack(vjust=0.5)) +
  labs(x = NULL, y = NULL) +
  coord_polar(theta = "y")

p2<- ggplot(subset(pie, type=="peak" & TRB=="TRB2"), aes(x = "TRB", y = count, fill = pie)) +
  geom_bar(stat = "identity", width = 3) +
  geom_text(aes(label = count), position = position_stack(vjust=0.5)) +
  labs(x = NULL, y = NULL) +
  coord_polar(theta = "y")

p3<- ggplot(subset(pie, type=="peak" & TRB=="TRB3"), aes(x = "TRB", y = count, fill = pie)) +
  geom_bar(stat = "identity", width = 3) +
  geom_text(aes(label = count), position = position_stack(vjust=0.5)) +
  labs(x = NULL, y = NULL) +
  coord_polar(theta = "y")
  

##########################################

pdf(file = "results/ChIP-seq/Piemotifs.pdf", width = 4, height = 12)
ggarrange(p1, p2,p3, ncol = 1, nrow = 3)
dev.off()


quit(save="no")

