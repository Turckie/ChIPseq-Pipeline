#!/usr/bin/Rscript

Araport11 <- read.table("Araport11.txt","\t",  header = TRUE, sep = "\t")

write.table(unique(subset(Araport11, TRB1==1 &type=="peak")[,5]), "Araport11.TRB1.AGI.txt", row.names=FALSE,col.names=FALSE,sep="\t", eol="\n",quote = FALSE)
write.table(unique(subset(Araport11, TRB2==1 &type=="peak")[,5]), "Araport11.TRB2.AGI.txt", row.names=FALSE,col.names=FALSE,sep="\t", eol="\n",quote = FALSE)
write.table(unique(subset(Araport11, TRB3==1 &type=="peak")[,5]), "Araport11.TRB3.AGI.txt", row.names=FALSE,col.names=FALSE,sep="\t", eol="\n",quote = FALSE)

write.table(unique(subset(Araport11, TRB1==1 &TRB2==0 &TRB3==0 &type=="peak")[,5]), "Araport11.TRB1_unique.AGI.txt", row.names=FALSE,col.names=FALSE,sep="\t", eol="\n",quote = FALSE)
write.table(unique(subset(Araport11, TRB1==0 &TRB2==1 &TRB3==0 &type=="peak")[,5]), "Araport11.TRB2_unique.AGI.txt", row.names=FALSE,col.names=FALSE,sep="\t", eol="\n",quote = FALSE)
write.table(unique(subset(Araport11, TRB1==0 &TRB2==0 &TRB3==1 &type=="peak")[,5]), "Araport11.TRB3_unique.AGI.txt", row.names=FALSE,col.names=FALSE,sep="\t", eol="\n",quote = FALSE)
write.table(unique(subset(Araport11, TRB1==1 &TRB2==1 &TRB3==1 &type=="peak")[,5]), "Araport11.common.AGI.txt", row.names=FALSE,col.names=FALSE,sep="\t", eol="\n",quote = FALSE)
write.table(unique(subset(Araport11, TRB1==1 &TRB2==1 &TRB3==0 &type=="peak")[,5]), "Araport11.TRB1_2_private.AGI.txt", row.names=FALSE,col.names=FALSE,sep="\t", eol="\n",quote = FALSE)
write.table(unique(subset(Araport11, TRB1==1 &TRB2==0 &TRB3==1 &type=="peak")[,5]), "Araport11.TRB1_3_private.AGI.txt", row.names=FALSE,col.names=FALSE,sep="\t", eol="\n",quote = FALSE)
write.table(unique(subset(Araport11, TRB1==0 &TRB2==1 &TRB3==1 &type=="peak")[,5]), "Araport11.TRB2_3_private.AGI.txt", row.names=FALSE,col.names=FALSE,sep="\t", eol="\n",quote = FALSE)

write.table(subset(Araport11, TRB1==1 &type=="peak")[,1:5], "Araport11.TRB1.bed", row.names=FALSE,col.names=FALSE,sep="\t", eol="\n",quote = FALSE)
write.table(subset(Araport11, TRB2==1 &type=="peak")[,1:5], "Araport11.TRB2.AGI.bed", row.names=FALSE,col.names=FALSE,sep="\t", eol="\n",quote = FALSE)
write.table(subset(Araport11, TRB3==1 &type=="peak")[,1:5], "Araport11.TRB3.AGI.bed", row.names=FALSE,col.names=FALSE,sep="\t", eol="\n",quote = FALSE)

write.table(subset(Araport11, TRB1==1 &TRB2==0 &TRB3==0 &type=="peak")[,1:5], "Araport11.TRB1_unique.AGI.bed", row.names=FALSE,col.names=FALSE,sep="\t", eol="\n",quote = FALSE)
write.table(subset(Araport11, TRB1==0 &TRB2==1 &TRB3==0 &type=="peak")[,1:5], "Araport11.TRB2_unique.AGI.bed", row.names=FALSE,col.names=FALSE,sep="\t", eol="\n",quote = FALSE)
write.table(subset(Araport11, TRB1==0 &TRB2==0 &TRB3==1 &type=="peak")[,1:5], "Araport11.TRB3_unique.AGI.bed", row.names=FALSE,col.names=FALSE,sep="\t", eol="\n",quote = FALSE)
write.table(subset(Araport11, TRB1==1 &TRB2==1 &TRB3==1 &type=="peak")[,1:5], "Araport11.common.AGI.bed", row.names=FALSE,col.names=FALSE,sep="\t", eol="\n",quote = FALSE)
write.table(subset(Araport11, TRB1==1 &TRB2==1 &TRB3==0 &type=="peak")[,1:5], "Araport11.TRB1_2_private.AGI.bed", row.names=FALSE,col.names=FALSE,sep="\t", eol="\n",quote = FALSE)
write.table(subset(Araport11, TRB1==1 &TRB2==0 &TRB3==1 &type=="peak")[,1:5], "Araport11.TRB1_3_private.AGI.bed", row.names=FALSE,col.names=FALSE,sep="\t", eol="\n",quote = FALSE)
write.table(subset(Araport11, TRB1==0 &TRB2==1 &TRB3==1 &type=="peak")[,1:5], "Araport11.TRB2_3_private.AGI.bed", row.names=FALSE,col.names=FALSE,sep="\t", eol="\n",quote = FALSE)




#write.table(subset(Araport11, k5=="cluster_1" &type=="peak")[,5], "Araport11.k5_k1.AGI.txt", row.names=FALSE,col.names=FALSE,sep="\t", eol="\n",quote = FALSE)
#write.table(subset(Araport11, k5=="cluster_2" &type=="peak")[,5], "Araport11.k5_k2.AGI.txt", row.names=FALSE,col.names=FALSE,sep="\t", eol="\n",quote = FALSE)
#write.table(subset(Araport11, k5=="cluster_3" &type=="peak")[,5], "Araport11.k5_k3.AGI.txt", row.names=FALSE,col.names=FALSE,sep="\t", eol="\n",quote = FALSE)
#write.table(subset(Araport11, k5=="cluster_4" &type=="peak")[,5], "Araport11.k5_k4.AGI.txt", row.names=FALSE,col.names=FALSE,sep="\t", eol="\n",quote = FALSE)
#write.table(subset(Araport11, k5=="cluster_5" &type=="peak")[,5], "Araport11.k5_k5.AGI.txt", row.names=FALSE,col.names=FALSE,sep="\t", eol="\n",quote = FALSE)
#write.table(subset(Araport11, k4=="cluster_1" &type=="peak")[,5], "Araport11.k4_k1.AGI.txt", row.names=FALSE,col.names=FALSE,sep="\t", eol="\n",quote = FALSE)
#write.table(subset(Araport11, k4=="cluster_2" &type=="peak")[,5], "Araport11.k4_k2.AGI.txt", row.names=FALSE,col.names=FALSE,sep="\t", eol="\n",quote = FALSE)
#write.table(subset(Araport11, k4=="cluster_3" &type=="peak")[,5], "Araport11.k4_k3.AGI.txt", row.names=FALSE,col.names=FALSE,sep="\t", eol="\n",quote = FALSE)
#write.table(subset(Araport11, k4=="cluster_4" &type=="peak")[,5], "Araport11.k4_k4.AGI.txt", row.names=FALSE,col.names=FALSE,sep="\t", eol="\n",quote = FALSE)

quit(save="no")