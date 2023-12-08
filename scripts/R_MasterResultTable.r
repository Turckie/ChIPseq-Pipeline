#!/usr/bin/Rscript

library("writexl")

KKfull <- read.table("data/ChIP-seq/RMASTER/Araport11.txt","\t",  header = TRUE, sep = "\t")

TRBcategories<-subset(KKfull, type=="peak" & category %in% c("TRB1_unique","TRB2_unique","TRB3_unique","TRB1_2_private", "TRB1_3_private","TRB2_3_private","TRB123_common"))

## ----Export data to Excel--------------------------------------------------------------------------------------------------------------------------


wb <- createWorkbook()

addWorksheet(wb, sheet = "TRB_target_categories", gridLines = TRUE)



writeData(wb, sheet = "TRB_target_categories", x = as.matrix(TRBcategories), rowNames = FALSE, withFilter = TRUE)


saveWorkbook(wb, "results/ChIP-seq/Supplemental_Data_2_ChIPseq.xlsx", overwrite = TRUE)



quit(save="no")