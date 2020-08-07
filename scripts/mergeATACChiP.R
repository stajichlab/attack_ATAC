library(plyr)
library(R.utils)

# goal: I used 03.1emptysite.sh to calculate coverage for ATAC and ChIP data
# Here I want to combine these datasets to the table I made with empty sites and also genomic features in mPingPromoter.bed
# just to mention, because of how things were analyzed, when it says A123 empty sites (mPing in A119 but not in A123), so data from A123 was used for A123 emptysites

setwd ("/rhome/ysun/bigdata/epigenome/attack_ATAC/results/Closest/")

Main = read.table("/rhome/ysun/bigdata/epigenome/attack_ATAC/results/Closest/mPingPromoter.bed", sep = '\t', header = FALSE, stringsAsFactors = FALSE)
Main = rename(Main, c("V1"="mPingChr" , "V2"="mPingSTART" , "V3"="mPingSTOP", "V4"="mPingDetails" , "V5"="featureChr" , "V6"="featureSTART" , "V7"="featureSTOP", "V8"="position", "V9"="strand" , "V10"="target_id"))


ATAC_A119 = read.table("/rhome/ysun/bigdata/epigenome/attack_ATAC/coverage/processed_empty/ATAC/ATAC_A119.A119_empty.100.regions.bed", sep = '\t', header = FALSE, stringsAsFactors = FALSE)            
ATAC_A119 = rename(ATAC_A119, c("V1"="mPingChr" , "V2"="mPingSTART" , "V3"="mPingSTOP", "V4"="mPingDetails" , "V5"="ATAC_A119.coverage"))
ATAC_A119 = ATAC_A119[,4:5]
ATAC_A123 = read.table("/rhome/ysun/bigdata/epigenome/attack_ATAC/coverage/processed_empty/ATAC/ATAC_A123.A123_empty.100.regions.bed", sep = '\t', header = FALSE, stringsAsFactors = FALSE)            
ATAC_A123 = rename(ATAC_A123, c("V1"="mPingChr" , "V2"="mPingSTART" , "V3"="mPingSTOP", "V4"="mPingDetails" , "V5"="ATAC_A123.coverage"))
ATAC_A123 = ATAC_A123[,4:5]
ATAC_EG4 = read.table("/rhome/ysun/bigdata/epigenome/attack_ATAC/coverage/processed_empty/ATAC/ATAC_EG4.EG4_empty.100.regions.bed", sep = '\t', header = FALSE, stringsAsFactors = FALSE)            
ATAC_EG4 = rename(ATAC_EG4, c("V1"="mPingChr" , "V2"="mPingSTART" , "V3"="mPingSTOP", "V4"="mPingDetails" , "V5"="ATAC_EG4.coverage"))
ATAC_EG4 = ATAC_EG4[,4:5]
ATAC_HEG4 = read.table("/rhome/ysun/bigdata/epigenome/attack_ATAC/coverage/processed_empty/ATAC/ATAC_HEG4.HEG4_empty.100.regions.bed", sep = '\t', header = FALSE, stringsAsFactors = FALSE)
ATAC_HEG4 = rename(ATAC_HEG4, c("V1"="mPingChr" , "V2"="mPingSTART" , "V3"="mPingSTOP", "V4"="mPingDetails" , "V5"="ATAC_HEG4.coverage"))
ATAC_HEG4 = ATAC_HEG4[,4:5]
ChIP27R1_A119 = read.table("/rhome/ysun/bigdata/epigenome/attack_ATAC/coverage/processed_empty/H3K27me3/H3K27me3_A119.rep1.A119_empty.100.regions.bed", sep = '\t', header = FALSE, stringsAsFactors = FALSE)  
ChIP27R1_A119 = rename(ChIP27R1_A119, c("V1"="mPingChr" , "V2"="mPingSTART" , "V3"="mPingSTOP", "V4"="mPingDetails" , "V5"="ChIP27R1_A119.coverage"))
ChIP27R1_A119 = ChIP27R1_A119[,4:5]
ChIP27R1_A123 = read.table("/rhome/ysun/bigdata/epigenome/attack_ATAC/coverage/processed_empty/H3K27me3/H3K27me3_A123.rep1.A123_empty.100.regions.bed", sep = '\t', header = FALSE, stringsAsFactors = FALSE) 
ChIP27R1_A123 = rename(ChIP27R1_A123, c("V1"="mPingChr" , "V2"="mPingSTART" , "V3"="mPingSTOP", "V4"="mPingDetails" , "V5"="ChIP27R1_A123.coverage"))
ChIP27R1_A123 = ChIP27R1_A123[,4:5]
ChIP27R1_EG4 = read.table("/rhome/ysun/bigdata/epigenome/attack_ATAC/coverage/processed_empty/H3K27me3/H3K27me3_EG4.rep1.EG4_empty.100.regions.bed", sep = '\t', header = FALSE, stringsAsFactors = FALSE)   
ChIP27R1_EG4 = rename(ChIP27R1_EG4, c("V1"="mPingChr" , "V2"="mPingSTART" , "V3"="mPingSTOP", "V4"="mPingDetails" , "V5"="ChIP27R1_EG4.coverage"))
ChIP27R1_EG4 = ChIP27R1_EG4[,4:5]
ChIP27R1_HEG4 = read.table("/rhome/ysun/bigdata/epigenome/attack_ATAC/coverage/processed_empty/H3K27me3/H3K27me3_HEG4.rep1.HEG4_empty.100.regions.bed", sep = '\t', header = FALSE, stringsAsFactors = FALSE) 
ChIP27R1_HEG4 = rename(ChIP27R1_HEG4, c("V1"="mPingChr" , "V2"="mPingSTART" , "V3"="mPingSTOP", "V4"="mPingDetails" , "V5"="ChIP27R1_HEG4.coverage"))
ChIP27R1_HEG4 = ChIP27R1_HEG4[,4:5]

ChIP27R2_A119 = read.table("/rhome/ysun/bigdata/epigenome/attack_ATAC/coverage/processed_empty/H3K27me3/H3K27me3_A119.rep2.A119_empty.100.regions.bed", sep = '\t', header = FALSE, stringsAsFactors = FALSE)   
ChIP27R2_A119 = rename(ChIP27R2_A119, c("V1"="mPingChr" , "V2"="mPingSTART" , "V3"="mPingSTOP", "V4"="mPingDetails" , "V5"="ChIP27R2_A119.coverage"))
ChIP27R2_A119 = ChIP27R2_A119[,4:5]
ChIP27R2_A123 = read.table("/rhome/ysun/bigdata/epigenome/attack_ATAC/coverage/processed_empty/H3K27me3/H3K27me3_A123.rep2.A123_empty.100.regions.bed", sep = '\t', header = FALSE, stringsAsFactors = FALSE)   
ChIP27R2_A123 = rename(ChIP27R2_A123, c("V1"="mPingChr" , "V2"="mPingSTART" , "V3"="mPingSTOP", "V4"="mPingDetails" , "V5"="ChIP27R2_A123.coverage"))
ChIP27R2_A123 = ChIP27R2_A123[,4:5]
ChIP27R2_EG4 = read.table("/rhome/ysun/bigdata/epigenome/attack_ATAC/coverage/processed_empty/H3K27me3/H3K27me3_EG4.rep2.EG4_empty.100.regions.bed", sep = '\t', header = FALSE, stringsAsFactors = FALSE)   
ChIP27R2_EG4 = rename(ChIP27R2_EG4, c("V1"="mPingChr" , "V2"="mPingSTART" , "V3"="mPingSTOP", "V4"="mPingDetails" , "V5"="ChIP27R2_EG4.coverage"))
ChIP27R2_EG4 = ChIP27R2_EG4[,4:5]
ChIP27R2_HEG4 = read.table("/rhome/ysun/bigdata/epigenome/attack_ATAC/coverage/processed_empty/H3K27me3/H3K27me3_HEG4.rep2.HEG4_empty.100.regions.bed", sep = '\t', header = FALSE, stringsAsFactors = FALSE) 
ChIP27R2_HEG4 = rename(ChIP27R2_HEG4, c("V1"="mPingChr" , "V2"="mPingSTART" , "V3"="mPingSTOP", "V4"="mPingDetails" , "V5"="ChIP27R2_HEG4.coverage"))
ChIP27R2_HEG4 = ChIP27R2_HEG4[,4:5]

ChIP36R1_A119 = read.table("/rhome/ysun/bigdata/epigenome/attack_ATAC/coverage/processed_empty/H3K36me3/H3K36me3_A119.rep1.A119_empty.100.regions.bed", sep = '\t', header = FALSE, stringsAsFactors = FALSE)  
ChIP36R1_A119 = rename(ChIP36R1_A119, c("V1"="mPingChr" , "V2"="mPingSTART" , "V3"="mPingSTOP", "V4"="mPingDetails" , "V5"="ChIP36R1_A119.coverage"))
ChIP36R1_A119 = ChIP36R1_A119[,4:5]
ChIP36R1_A123 = read.table("/rhome/ysun/bigdata/epigenome/attack_ATAC/coverage/processed_empty/H3K36me3/H3K36me3_A123.rep1.A123_empty.100.regions.bed", sep = '\t', header = FALSE, stringsAsFactors = FALSE) 
ChIP36R1_A123 = rename(ChIP36R1_A123, c("V1"="mPingChr" , "V2"="mPingSTART" , "V3"="mPingSTOP", "V4"="mPingDetails" , "V5"="ChIP36R1_A123.coverage"))
ChIP36R1_A123 = ChIP36R1_A123[,4:5]
ChIP36R1_EG4 = read.table("/rhome/ysun/bigdata/epigenome/attack_ATAC/coverage/processed_empty/H3K36me3/H3K36me3_EG4.rep1.EG4_empty.100.regions.bed", sep = '\t', header = FALSE, stringsAsFactors = FALSE)   
ChIP36R1_EG4 = rename(ChIP36R1_EG4, c("V1"="mPingChr" , "V2"="mPingSTART" , "V3"="mPingSTOP", "V4"="mPingDetails" , "V5"="ChIP36R1_EG4.coverage"))
ChIP36R1_EG4 = ChIP36R1_EG4[,4:5]
ChIP36R1_HEG4 = read.table("/rhome/ysun/bigdata/epigenome/attack_ATAC/coverage/processed_empty/H3K36me3/H3K36me3_HEG4.rep1.HEG4_empty.100.regions.bed", sep = '\t', header = FALSE, stringsAsFactors = FALSE) 
ChIP36R1_HEG4 = rename(ChIP36R1_HEG4, c("V1"="mPingChr" , "V2"="mPingSTART" , "V3"="mPingSTOP", "V4"="mPingDetails" , "V5"="ChIP36R1_HEG4.coverage"))
ChIP36R1_HEG4 = ChIP36R1_HEG4[,4:5]

ChIP36R2_A119 = read.table("/rhome/ysun/bigdata/epigenome/attack_ATAC/coverage/processed_empty/H3K36me3/H3K36me3_A119.rep2.A119_empty.100.regions.bed", sep = '\t', header = FALSE, stringsAsFactors = FALSE)   
ChIP36R2_A119 = rename(ChIP36R2_A119, c("V1"="mPingChr" , "V2"="mPingSTART" , "V3"="mPingSTOP", "V4"="mPingDetails" , "V5"="ChIP36R2_A119.coverage"))
ChIP36R2_A119 = ChIP36R2_A119[,4:5]
ChIP36R2_A123 = read.table("/rhome/ysun/bigdata/epigenome/attack_ATAC/coverage/processed_empty/H3K36me3/H3K36me3_A123.rep2.A123_empty.100.regions.bed", sep = '\t', header = FALSE, stringsAsFactors = FALSE)   
ChIP36R2_A123 = rename(ChIP36R2_A123, c("V1"="mPingChr" , "V2"="mPingSTART" , "V3"="mPingSTOP", "V4"="mPingDetails" , "V5"="ChIP36R2_A123.coverage"))
ChIP36R2_A123 = ChIP36R2_A123[,4:5]
ChIP36R2_EG4 = read.table("/rhome/ysun/bigdata/epigenome/attack_ATAC/coverage/processed_empty/H3K36me3/H3K36me3_EG4.rep2.EG4_empty.100.regions.bed", sep = '\t', header = FALSE, stringsAsFactors = FALSE)   
ChIP36R2_EG4 = rename(ChIP36R2_EG4, c("V1"="mPingChr" , "V2"="mPingSTART" , "V3"="mPingSTOP", "V4"="mPingDetails" , "V5"="ChIP36R2_EG4.coverage"))
ChIP36R2_EG4 = ChIP36R2_EG4[,4:5]
ChIP36R2_HEG4 = read.table("/rhome/ysun/bigdata/epigenome/attack_ATAC/coverage/processed_empty/H3K36me3/H3K36me3_HEG4.rep2.HEG4_empty.100.regions.bed", sep = '\t', header = FALSE, stringsAsFactors = FALSE) 
ChIP36R2_HEG4 = rename(ChIP36R2_HEG4, c("V1"="mPingChr" , "V2"="mPingSTART" , "V3"="mPingSTOP", "V4"="mPingDetails" , "V5"="ChIP36R2_HEG4.coverage"))
ChIP36R2_HEG4 = ChIP36R2_HEG4[,4:5]

ChIP56R1_A119 = read.table("/rhome/ysun/bigdata/epigenome/attack_ATAC/coverage/processed_empty/H3K56ac/H3K56ac_A119.rep1.A119_empty.100.regions.bed", sep = '\t', header = FALSE, stringsAsFactors = FALSE)  
ChIP56R1_A119 = rename(ChIP56R1_A119, c("V1"="mPingChr" , "V2"="mPingSTART" , "V3"="mPingSTOP", "V4"="mPingDetails" , "V5"="ChIP56R1_A119.coverage"))
ChIP56R1_A119 = ChIP56R1_A119[,4:5]
ChIP56R1_A123 = read.table("/rhome/ysun/bigdata/epigenome/attack_ATAC/coverage/processed_empty/H3K56ac/H3K56ac_A123.rep1.A123_empty.100.regions.bed", sep = '\t', header = FALSE, stringsAsFactors = FALSE) 
ChIP56R1_A123 = rename(ChIP56R1_A123, c("V1"="mPingChr" , "V2"="mPingSTART" , "V3"="mPingSTOP", "V4"="mPingDetails" , "V5"="ChIP56R1_A123.coverage"))
ChIP56R1_A123 = ChIP56R1_A123[,4:5]
ChIP56R1_EG4 = read.table("/rhome/ysun/bigdata/epigenome/attack_ATAC/coverage/processed_empty/H3K56ac/H3K56ac_EG4.rep1.EG4_empty.100.regions.bed", sep = '\t', header = FALSE, stringsAsFactors = FALSE)   
ChIP56R1_EG4 = rename(ChIP56R1_EG4, c("V1"="mPingChr" , "V2"="mPingSTART" , "V3"="mPingSTOP", "V4"="mPingDetails" , "V5"="ChIP56R1_EG4.coverage"))
ChIP56R1_EG4 = ChIP56R1_EG4[,4:5]
ChIP56R1_HEG4 = read.table("/rhome/ysun/bigdata/epigenome/attack_ATAC/coverage/processed_empty/H3K56ac/H3K56ac_HEG4.rep1.HEG4_empty.100.regions.bed", sep = '\t', header = FALSE, stringsAsFactors = FALSE) 
ChIP56R1_HEG4 = rename(ChIP56R1_HEG4, c("V1"="mPingChr" , "V2"="mPingSTART" , "V3"="mPingSTOP", "V4"="mPingDetails" , "V5"="ChIP56R1_HEG4.coverage"))
ChIP56R1_HEG4 = ChIP56R1_HEG4[,4:5]
ChIP56R2_A119 = read.table("/rhome/ysun/bigdata/epigenome/attack_ATAC/coverage/processed_empty/H3K56ac/H3K56ac_A119.rep2.A119_empty.100.regions.bed", sep = '\t', header = FALSE, stringsAsFactors = FALSE)   
ChIP56R2_A119 = rename(ChIP56R2_A119, c("V1"="mPingChr" , "V2"="mPingSTART" , "V3"="mPingSTOP", "V4"="mPingDetails" , "V5"="ChIP56R2_A119.coverage"))
ChIP56R2_A119 = ChIP56R2_A119[,4:5]
ChIP56R2_A123 = read.table("/rhome/ysun/bigdata/epigenome/attack_ATAC/coverage/processed_empty/H3K56ac/H3K56ac_A123.rep2.A123_empty.100.regions.bed", sep = '\t', header = FALSE, stringsAsFactors = FALSE)   
ChIP56R2_A123 = rename(ChIP56R2_A123, c("V1"="mPingChr" , "V2"="mPingSTART" , "V3"="mPingSTOP", "V4"="mPingDetails" , "V5"="ChIP56R2_A123.coverage"))
ChIP56R2_A123 = ChIP56R2_A123[,4:5]
ChIP56R2_EG4 = read.table("/rhome/ysun/bigdata/epigenome/attack_ATAC/coverage/processed_empty/H3K56ac/H3K56ac_EG4.rep2.EG4_empty.100.regions.bed", sep = '\t', header = FALSE, stringsAsFactors = FALSE)   
ChIP56R2_EG4 = rename(ChIP56R2_EG4, c("V1"="mPingChr" , "V2"="mPingSTART" , "V3"="mPingSTOP", "V4"="mPingDetails" , "V5"="ChIP56R2_EG4.coverage"))
ChIP56R2_EG4 = ChIP56R2_EG4[,4:5]
ChIP56R2_HEG4 = read.table("/rhome/ysun/bigdata/epigenome/attack_ATAC/coverage/processed_empty/H3K56ac/H3K56ac_HEG4.rep2.HEG4_empty.100.regions.bed", sep = '\t', header = FALSE, stringsAsFactors = FALSE) 
ChIP56R2_HEG4 = rename(ChIP56R2_HEG4, c("V1"="mPingChr" , "V2"="mPingSTART" , "V3"="mPingSTOP", "V4"="mPingDetails" , "V5"="ChIP56R2_HEG4.coverage"))
ChIP56R2_HEG4 = ChIP56R2_HEG4[,4:5]
############################################################################################################

joinmain = left_join(Main, ATAC_A119, by = "mPingDetails")
joinA123 = left_join(joinmain, ATAC_A123, by = "mPingDetails")
joinEG4 = left_join(joinA123, ATAC_EG4, by = "mPingDetails")
joinHEG4 = left_join(joinEG4, ATAC_HEG4, by = "mPingDetails")
joinChIP27R1_A119 = left_join(joinHEG4, ChIP27R1_A119, by = "mPingDetails")
joinChIP27R1_A123 = left_join(joinChIP27R1_A119, ChIP27R1_A123, by = "mPingDetails")
joinChIP27R1_EG4 = left_join(joinChIP27R1_A123, ChIP27R1_EG4, by = "mPingDetails")
joinChIP27R1_HEG4 = left_join(joinChIP27R1_EG4, ChIP27R1_HEG4, by = "mPingDetails")
joinChIP27R2_A119 = left_join(joinChIP27R1_HEG4, ChIP27R2_A119, by = "mPingDetails")
joinChIP27R2_A123 = left_join(joinChIP27R2_A119, ChIP27R2_A123, by = "mPingDetails")
joinChIP27R2_EG4 = left_join(joinChIP27R2_A123, ChIP27R2_EG4, by = "mPingDetails")
joinChIP27R2_HEG4 = left_join(joinChIP27R2_EG4, ChIP27R2_HEG4, by = "mPingDetails")
joinChIP36R1_A119 = left_join(joinChIP27R2_HEG4, ChIP36R1_A119, by = "mPingDetails")
joinChIP36R1_A123 = left_join(joinChIP36R1_A119, ChIP36R1_A123, by = "mPingDetails")
joinChIP36R1_EG4 = left_join(joinChIP36R1_A123, ChIP36R1_EG4, by = "mPingDetails")
joinChIP36R1_HEG4 = left_join(joinChIP36R1_EG4, ChIP36R1_HEG4, by = "mPingDetails")
joinChIP36R2_A119 = left_join(joinChIP36R1_HEG4, ChIP36R2_A119, by = "mPingDetails")
joinChIP36R2_A123 = left_join(joinChIP36R2_A119, ChIP36R2_A123, by = "mPingDetails")
joinChIP36R2_EG4 = left_join(joinChIP36R2_A123, ChIP36R2_EG4, by = "mPingDetails")
joinChIP36R2_HEG4 = left_join(joinChIP36R2_EG4, ChIP36R2_HEG4, by = "mPingDetails")
joinChIP56R1_A119 = left_join(joinChIP36R2_HEG4, ChIP56R1_A119, by = "mPingDetails")
joinChIP56R1_A123 = left_join(joinChIP56R1_A119, ChIP56R1_A123, by = "mPingDetails")
joinChIP56R1_EG4 = left_join(joinChIP56R1_A123, ChIP56R1_EG4, by = "mPingDetails")
joinChIP56R1_HEG4 = left_join(joinChIP56R1_EG4, ChIP56R1_HEG4, by = "mPingDetails")
joinChIP56R2_A119 = left_join(joinChIP56R1_HEG4, ChIP56R2_A119, by = "mPingDetails")
joinChIP56R2_A123 = left_join(joinChIP56R2_A119, ChIP56R2_A123, by = "mPingDetails")
joinChIP56R2_EG4 = left_join(joinChIP56R2_A123, ChIP56R2_EG4, by = "mPingDetails")

joinFINAL = left_join(joinChIP56R2_EG4, ChIP56R2_HEG4, by = "mPingDetails")

write.table(joinFINAL,"emptysiteATACCHIP.txt",sep="\t",row.names=FALSE)
