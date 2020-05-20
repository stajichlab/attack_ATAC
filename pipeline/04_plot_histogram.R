#!/usr/bin/env Rscript

setwd("/rhome/ysun/bigdata/epigenome/attack_ATAC/coverage/processed_empty")

ATAC.A119_empty.1000_window = read.table("ATAC.A119_empty.1000_window.regions.bed.gz",header=F)    
ATAC.A119_empty.500_window = read.table("ATAC.A119_empty.500_window.regions.bed.gz",header=F)    
ATAC.A123_empty.1000_window = read.table("ATAC.A123_empty.1000_window.regions.bed.gz",header=F)    
ATAC.A123_empty.500_window = read.table("ATAC.A123_empty.500_window.regions.bed.gz",header=F)    
ATAC.EG4_empty.1000_window = read.table("ATAC.EG4_empty.1000_window.regions.bed.gz",header=F)    
ATAC.EG4_empty.500_window = read.table("ATAC.EG4_empty.500_window.regions.bed.gz",header=F)    
ATAC.HEG4_empty.1000_window = read.table("ATAC.HEG4_empty.1000_window.regions.bed.gz",header=F)    
ATAC.HEG4_empty.500_window = read.table("ATAC.HEG4_empty.500_window.regions.bed.gz",header=F)    
H3K27me3.A119_empty.1000_window = read.table("H3K27me3.A119_empty.1000_window.regions.bed.gz",header=F)    
H3K27me3.A119_empty.500_window = read.table("H3K27me3.A119_empty.500_window.regions.bed.gz",header=F)    
H3K27me3.A123_empty.1000_window = read.table("H3K27me3.A123_empty.1000_window.regions.bed.gz",header=F)    
H3K27me3.A123_empty.500_window = read.table("H3K27me3.A123_empty.500_window.regions.bed.gz",header=F)    
H3K27me3.EG4_empty.1000_window = read.table("H3K27me3.EG4_empty.1000_window.regions.bed.gz",header=F)    
H3K27me3.EG4_empty.500_window = read.table("H3K27me3.EG4_empty.500_window.regions.bed.gz",header=F)    
H3K27me3.HEG4_empty.1000_window = read.table("H3K27me3.HEG4_empty.1000_window.regions.bed.gz",header=F)    
H3K27me3.HEG4_empty.500_window = read.table("H3K27me3.HEG4_empty.500_window.regions.bed.gz",header=F)    
H3K36me3.A119_empty.1000_window = read.table("H3K36me3.A119_empty.1000_window.regions.bed.gz",header=F)    
H3K36me3.A119_empty.500_window = read.table("H3K36me3.A119_empty.500_window.regions.bed.gz",header=F)    
H3K36me3.A123_empty.1000_window = read.table("H3K36me3.A123_empty.1000_window.regions.bed.gz",header=F)    
H3K36me3.A123_empty.500_window = read.table("H3K36me3.A123_empty.500_window.regions.bed.gz",header=F)    
H3K36me3.EG4_empty.1000_window = read.table("H3K36me3.EG4_empty.1000_window.regions.bed.gz",header=F)    
H3K36me3.EG4_empty.500_window = read.table("H3K36me3.EG4_empty.500_window.regions.bed.gz",header=F)    
H3K36me3.HEG4_empty.1000_window = read.table("H3K36me3.HEG4_empty.1000_window.regions.bed.gz",header=F)    
H3K36me3.HEG4_empty.500_window = read.table("H3K36me3.HEG4_empty.500_window.regions.bed.gz",header=F)    
H3K56ac.A119_empty.1000_window = read.table("H3K56ac.A119_empty.1000_window.regions.bed.gz",header=F)    
H3K56ac.A119_empty.500_window = read.table("H3K56ac.A119_empty.500_window.regions.bed.gz",header=F)    
H3K56ac.A123_empty.1000_window = read.table("H3K56ac.A123_empty.1000_window.regions.bed.gz",header=F)    
H3K56ac.A123_empty.500_window = read.table("H3K56ac.A123_empty.500_window.regions.bed.gz",header=F)    
H3K56ac.EG4_empty.1000_window = read.table("H3K56ac.EG4_empty.1000_window.regions.bed.gz",header=F)    
H3K56ac.EG4_empty.500_window = read.table("H3K56ac.EG4_empty.500_window.regions.bed.gz",header=F)    
H3K56ac.HEG4_empty.1000_window = read.table("H3K56ac.HEG4_empty.1000_window.regions.bed.gz",header=F)    
H3K56ac.HEG4_empty.500_window = read.table("H3K56ac.HEG4_empty.500_window.regions.bed.gz",header=F)    


colnames(ATAC.A119_empty.1000_window) <- c("CHR","START","END","NAME","COVERAGE")
colnames(ATAC.A119_empty.500_window) <- c("CHR","START","END","NAME","COVERAGE")
colnames(ATAC.A123_empty.1000_window) <- c("CHR","START","END","NAME","COVERAGE")
colnames(ATAC.A123_empty.500_window) <- c("CHR","START","END","NAME","COVERAGE")
colnames(ATAC.EG4_empty.1000_window) <- c("CHR","START","END","NAME","COVERAGE")
colnames(ATAC.EG4_empty.500_window) <- c("CHR","START","END","NAME","COVERAGE")
colnames(ATAC.HEG4_empty.1000_window) <- c("CHR","START","END","NAME","COVERAGE")
colnames(ATAC.HEG4_empty.500_window) <- c("CHR","START","END","NAME","COVERAGE")
colnames(H3K27me3.A119_empty.1000_window) <- c("CHR","START","END","NAME","COVERAGE")
colnames(H3K27me3.A119_empty.500_window) <- c("CHR","START","END","NAME","COVERAGE")
colnames(H3K27me3.A123_empty.1000_window) <- c("CHR","START","END","NAME","COVERAGE")
colnames(H3K27me3.A123_empty.500_window) <- c("CHR","START","END","NAME","COVERAGE")
colnames(H3K27me3.EG4_empty.1000_window) <- c("CHR","START","END","NAME","COVERAGE")
colnames(H3K27me3.EG4_empty.500_window) <- c("CHR","START","END","NAME","COVERAGE")
colnames(H3K27me3.HEG4_empty.1000_window) <- c("CHR","START","END","NAME","COVERAGE")
colnames(H3K27me3.HEG4_empty.500_window) <- c("CHR","START","END","NAME","COVERAGE")
colnames(H3K36me3.A119_empty.1000_window) <- c("CHR","START","END","NAME","COVERAGE")
colnames(H3K36me3.A119_empty.500_window) <- c("CHR","START","END","NAME","COVERAGE")
colnames(H3K36me3.A123_empty.1000_window) <- c("CHR","START","END","NAME","COVERAGE")
colnames(H3K36me3.A123_empty.500_window) <- c("CHR","START","END","NAME","COVERAGE")
colnames(H3K36me3.EG4_empty.1000_window) <- c("CHR","START","END","NAME","COVERAGE")
colnames(H3K36me3.EG4_empty.500_window) <- c("CHR","START","END","NAME","COVERAGE")
colnames(H3K36me3.HEG4_empty.1000_window) <- c("CHR","START","END","NAME","COVERAGE")
colnames(H3K36me3.HEG4_empty.500_window) <- c("CHR","START","END","NAME","COVERAGE")
colnames(H3K56ac.A119_empty.1000_window) <- c("CHR","START","END","NAME","COVERAGE")
colnames(H3K56ac.A119_empty.500_window) <- c("CHR","START","END","NAME","COVERAGE")
colnames(H3K56ac.A123_empty.1000_window) <- c("CHR","START","END","NAME","COVERAGE")
colnames(H3K56ac.A123_empty.500_window) <- c("CHR","START","END","NAME","COVERAGE")
colnames(H3K56ac.EG4_empty.1000_window) <- c("CHR","START","END","NAME","COVERAGE")
colnames(H3K56ac.EG4_empty.500_window) <- c("CHR","START","END","NAME","COVERAGE")
colnames(H3K56ac.HEG4_empty.1000_window) <- c("CHR","START","END","NAME","COVERAGE")
colnames(H3K56ac.HEG4_empty.500_window) <- c("CHR","START","END","NAME","COVERAGE")
########################################################################################################################
pdf("ATAC.A119_empty.1000_window.pdf")
hist(ATAC.A119_empty.1000_window$COVERAGE, 
     main="Histogram of ATAC data for empty sites in A119 with bins 1kb wide", 
     xlab="Coverage", 
     ylab="Frequency, bin = 100", 
     border="blue", 
     col="green",
     xlim=c(0,15),
     ylim=c(0,60),
     las=1, 
     breaks=100)
dev.off()
########################################################################################################################
pdf("ATAC.A123_empty.1000_window.pdf")
hist(ATAC.A123_empty.1000_window$COVERAGE, 
     main="Histogram of ATAC data for empty sites in A123 with bins 1kb wide", 
     xlab="Coverage", 
     ylab="Frequency, bin = 100", 
     border="blue", 
     col="orange",
     xlim=c(0,15),
     ylim=c(0,60),
     las=1, 
     breaks=100)
dev.off()
########################################################################################################################
pdf("ATAC.EG4_empty.1000_window.pdf")
hist(ATAC.EG4_empty.1000_window$COVERAGE, 
     main="Histogram of ATAC data for empty sites in EG4 with bins 1kb wide", 
     xlab="Coverage", 
     ylab="Frequency, bin = 100", 
     border="blue", 
     col="yellow",
     xlim=c(0,15),
     ylim=c(0,60),
     las=1, 
     breaks=100)
dev.off()
########################################################################################################################
pdf("ATAC.HEG4_empty.1000_window.pdf")
hist(ATAC.HEG4_empty.1000_window$COVERAGE, 
     main="Histogram of ATAC data for empty sites in HEG4 with bins 1kb wide", 
     xlab="Coverage", 
     ylab="Frequency, bin = 100", 
     border="blue", 
     col="red",
     xlim=c(0,15),
     ylim=c(0,60),
     las=1, 
     breaks=100)
dev.off()
########################################################################################################################
########################################################################################################################
pdf("ATAC.A119_empty.500_window.pdf")
hist(ATAC.A119_empty.500_window$COVERAGE, 
     main="Histogram of ATAC data for empty sites in A119 with bins .5kb wide", 
     xlab="Coverage", 
     ylab="Frequency, bin = 100", 
     border="blue", 
     col="green",
     xlim=c(0,15),
     ylim=c(0,60),
     las=1, 
     breaks=100)
dev.off()
########################################################################################################################
pdf("ATAC.A123_empty.500_window.pdf")
hist(ATAC.A123_empty.500_window$COVERAGE, 
     main="Histogram of ATAC data for empty sites in A123 with bins .5kb wide", 
     xlab="Coverage", 
     ylab="Frequency, bin = 100", 
     border="blue", 
     col="orange",
     xlim=c(0,15),
     ylim=c(0,60),
     las=1, 
     breaks=100)
dev.off()
########################################################################################################################
pdf("ATAC.EG4_empty.500_window.pdf")
hist(ATAC.EG4_empty.500_window$COVERAGE, 
     main="Histogram of ATAC data for empty sites in EG4 with bins .5kb wide", 
     xlab="Coverage", 
     ylab="Frequency, bin = 100", 
     border="blue", 
     col="yellow",
     xlim=c(0,15),
     ylim=c(0,60),
     las=1, 
     breaks=100)
dev.off()
########################################################################################################################
pdf("ATAC.HEG4_empty.500_window.pdf")
hist(ATAC.HEG4_empty.500_window$COVERAGE, 
     main="Histogram of ATAC data for empty sites in HEG4 with bins .5kb wide", 
     xlab="Coverage", 
     ylab="Frequency, bin = 100", 
     border="blue", 
     col="red",
     xlim=c(0,15),
     ylim=c(0,60),
     las=1, 
     breaks=100)
dev.off()
########################################################################################################################
########################################################################################################################
pdf("H3K27me3.A119_empty.1000_window.pdf")
hist(H3K27me3.A119_empty.1000_window$COVERAGE, 
     main="Histogram of H3K27me3 data for empty sites in A119 with bins 1kb wide", 
     xlab="Coverage", 
     ylab="Frequency, bin = 100", 
     border="blue", 
     col="green",
     xlim=c(0,80),
     ylim=c(0,15),
     las=1, 
     breaks=100)
dev.off()
########################################################################################################################
pdf("H3K27me3.A123_empty.1000_window.pdf")
hist(H3K27me3.A123_empty.1000_window$COVERAGE, 
     main="Histogram of H3K27me3 data for empty sites in A123 with bins 1kb wide", 
     xlab="Coverage", 
     ylab="Frequency, bin = 100", 
     border="blue", 
     col="orange",
     xlim=c(0,80),
     ylim=c(0,15),
     las=1, 
     breaks=100)
dev.off()
########################################################################################################################
pdf("H3K27me3.EG4_empty.1000_window.pdf")
hist(H3K27me3.EG4_empty.1000_window$COVERAGE, 
     main="Histogram of H3K27me3 data for empty sites in EG4 with bins 1kb wide", 
     xlab="Coverage", 
     ylab="Frequency, bin = 100", 
     border="blue", 
     col="yellow",
     xlim=c(0,80),
     ylim=c(0,15),
     las=1, 
     breaks=100)
dev.off()
########################################################################################################################
pdf("H3K27me3.HEG4_empty.1000_window.pdf")
hist(H3K27me3.HEG4_empty.1000_window$COVERAGE, 
     main="Histogram of H3K27me3 data for empty sites in HEG4 with bins 1kb wide", 
     xlab="Coverage", 
     ylab="Frequency, bin = 100", 
     border="blue", 
     col="red",
     xlim=c(0,80),
     ylim=c(0,15),
     las=1, 
     breaks=100)
dev.off()
########################################################################################################################
########################################################################################################################
pdf("H3K27me3.A119_empty.500_window.pdf")
hist(H3K27me3.A119_empty.500_window$COVERAGE, 
     main="Histogram of H3K27me3 data for empty sites in A119 with bins .5kb wide", 
     xlab="Coverage", 
     ylab="Frequency, bin = 100", 
     border="blue", 
     col="green",
     xlim=c(0,80),
     ylim=c(0,15),
     las=1, 
     breaks=100)
dev.off()
########################################################################################################################
pdf("H3K27me3.A123_empty.500_window.pdf")
hist(H3K27me3.A123_empty.500_window$COVERAGE, 
     main="Histogram of H3K27me3 data for empty sites in A123 with bins .5kb wide", 
     xlab="Coverage", 
     ylab="Frequency, bin = 100", 
     border="blue", 
     col="orange",
     xlim=c(0,80),
     ylim=c(0,15),
     las=1, 
     breaks=100)
dev.off()
########################################################################################################################
pdf("H3K27me3.EG4_empty.500_window.pdf")
hist(H3K27me3.EG4_empty.500_window$COVERAGE, 
     main="Histogram of H3K27me3 data for empty sites in EG4 with bins .5kb wide", 
     xlab="Coverage", 
     ylab="Frequency, bin = 100", 
     border="blue", 
     col="yellow",
     xlim=c(0,80),
     ylim=c(0,15),
     las=1, 
     breaks=100)
dev.off()
########################################################################################################################
pdf("H3K27me3.HEG4_empty.500_window.pdf")
hist(H3K27me3.HEG4_empty.500_window$COVERAGE, 
     main="Histogram of H3K27me3 data for empty sites in HEG4 with bins .5kb wide", 
     xlab="Coverage", 
     ylab="Frequency, bin = 100", 
     border="blue", 
     col="red",
     xlim=c(0,80),
     ylim=c(0,15),
     las=1, 
     breaks=100)
dev.off()
########################################################################################################################

########################################################################################################################
pdf("H3K36me3.A119_empty.1000_window.pdf")
hist(H3K36me3.A119_empty.1000_window$COVERAGE, 
     main="Histogram of H3K36me3 data for empty sites in A119 with bins 1kb wide", 
     xlab="Coverage", 
     ylab="Frequency, bin = 100", 
     border="blue", 
     col="green",
     xlim=c(0,80),
     ylim=c(0,15),
     las=1, 
     breaks=100)
dev.off()
########################################################################################################################
pdf("H3K36me3.A123_empty.1000_window.pdf")
hist(H3K36me3.A123_empty.1000_window$COVERAGE, 
     main="Histogram of H3K36me3 data for empty sites in A123 with bins 1kb wide", 
     xlab="Coverage", 
     ylab="Frequency, bin = 100", 
     border="blue", 
     col="orange",
     xlim=c(0,80),
     ylim=c(0,15),
     las=1, 
     breaks=100)
dev.off()
########################################################################################################################
pdf("H3K36me3.EG4_empty.1000_window.pdf")
hist(H3K36me3.EG4_empty.1000_window$COVERAGE, 
     main="Histogram of H3K36me3 data for empty sites in EG4 with bins 1kb wide", 
     xlab="Coverage", 
     ylab="Frequency, bin = 100", 
     border="blue", 
     col="yellow",
     xlim=c(0,80),
     ylim=c(0,15),
     las=1, 
     breaks=100)
dev.off()
########################################################################################################################
pdf("H3K36me3.HEG4_empty.1000_window.pdf")
hist(H3K36me3.HEG4_empty.1000_window$COVERAGE, 
     main="Histogram of H3K36me3 data for empty sites in HEG4 with bins 1kb wide", 
     xlab="Coverage", 
     ylab="Frequency, bin = 100", 
     border="blue", 
     col="red",
     xlim=c(0,80),
     ylim=c(0,15),
     las=1, 
     breaks=100)
dev.off()
########################################################################################################################
########################################################################################################################
pdf("H3K36me3.A119_empty.500_window.pdf")
hist(H3K36me3.A119_empty.500_window$COVERAGE, 
     main="Histogram of H3K36me3 data for empty sites in A119 with bins .5kb wide", 
     xlab="Coverage", 
     ylab="Frequency, bin = 100", 
     border="blue", 
     col="green",
     xlim=c(0,80),
     ylim=c(0,15),
     las=1, 
     breaks=100)
dev.off()
########################################################################################################################
pdf("H3K36me3.A123_empty.500_window.pdf")
hist(H3K36me3.A123_empty.500_window$COVERAGE, 
     main="Histogram of H3K36me3 data for empty sites in A123 with bins .5kb wide", 
     xlab="Coverage", 
     ylab="Frequency, bin = 100", 
     border="blue", 
     col="orange",
     xlim=c(0,80),
     ylim=c(0,15),
     las=1, 
     breaks=100)
dev.off()
########################################################################################################################
pdf("H3K36me3.EG4_empty.500_window.pdf")
hist(H3K36me3.EG4_empty.500_window$COVERAGE, 
     main="Histogram of H3K36me3 data for empty sites in EG4 with bins .5kb wide", 
     xlab="Coverage", 
     ylab="Frequency, bin = 100", 
     border="blue", 
     col="yellow",
     xlim=c(0,80),
     ylim=c(0,15),
     las=1, 
     breaks=100)
dev.off()
########################################################################################################################
pdf("H3K36me3.HEG4_empty.500_window.pdf")
hist(H3K36me3.HEG4_empty.500_window$COVERAGE, 
     main="Histogram of H3K36me3 data for empty sites in HEG4 with bins .5kb wide", 
     xlab="Coverage", 
     ylab="Frequency, bin = 100", 
     border="blue", 
     col="red",
     xlim=c(0,80),
     ylim=c(0,15),
     las=1, 
     breaks=100)
dev.off()
########################################################################################################################
########################################################################################################################
pdf("H3K56ac.A119_empty.1000_window.pdf")
hist(H3K56ac.A119_empty.1000_window$COVERAGE, 
     main="Histogram of H3K56ac data for empty sites in A119 with bins 1kb wide", 
     xlab="Coverage", 
     ylab="Frequency, bin = 100", 
     border="blue", 
     col="green",
     xlim=c(0,80),
     ylim=c(0,15),
     las=1, 
     breaks=100)
dev.off()
########################################################################################################################
pdf("H3K56ac.A123_empty.1000_window.pdf")
hist(H3K56ac.A123_empty.1000_window$COVERAGE, 
     main="Histogram of H3K56ac data for empty sites in A123 with bins 1kb wide", 
     xlab="Coverage", 
     ylab="Frequency, bin = 100", 
     border="blue", 
     col="orange",
     xlim=c(0,80),
     ylim=c(0,15),
     las=1, 
     breaks=100)
dev.off()
########################################################################################################################
pdf("H3K56ac.EG4_empty.1000_window.pdf")
hist(H3K56ac.EG4_empty.1000_window$COVERAGE, 
     main="Histogram of H3K56ac data for empty sites in EG4 with bins 1kb wide", 
     xlab="Coverage", 
     ylab="Frequency, bin = 100", 
     border="blue", 
     col="yellow",
     xlim=c(0,80),
     ylim=c(0,15),
     las=1, 
     breaks=100)
dev.off()
########################################################################################################################
pdf("H3K56ac.HEG4_empty.1000_window.pdf")
hist(H3K56ac.HEG4_empty.1000_window$COVERAGE, 
     main="Histogram of H3K56ac data for empty sites in HEG4 with bins 1kb wide", 
     xlab="Coverage", 
     ylab="Frequency, bin = 100", 
     border="blue", 
     col="red",
     xlim=c(0,80),
     ylim=c(0,15),
     las=1, 
     breaks=100)
dev.off()
########################################################################################################################
########################################################################################################################
pdf("H3K56ac.A119_empty.500_window.pdf")
hist(H3K56ac.A119_empty.500_window$COVERAGE, 
     main="Histogram of H3K56ac data for empty sites in A119 with bins .5kb wide", 
     xlab="Coverage", 
     ylab="Frequency, bin = 100", 
     border="blue", 
     col="green",
     xlim=c(0,80),
     ylim=c(0,15),
     las=1, 
     breaks=100)
dev.off()
########################################################################################################################
pdf("H3K56ac.A123_empty.500_window.pdf")
hist(H3K56ac.A123_empty.500_window$COVERAGE, 
     main="Histogram of H3K56ac data for empty sites in A123 with bins .5kb wide", 
     xlab="Coverage", 
     ylab="Frequency, bin = 100", 
     border="blue", 
     col="orange",
     xlim=c(0,80),
     ylim=c(0,15),
     las=1, 
     breaks=100)
dev.off()
########################################################################################################################
pdf("H3K56ac.EG4_empty.500_window.pdf")
hist(H3K56ac.EG4_empty.500_window$COVERAGE, 
     main="Histogram of H3K56ac data for empty sites in EG4 with bins .5kb wide", 
     xlab="Coverage", 
     ylab="Frequency, bin = 100", 
     border="blue", 
     col="yellow",
     xlim=c(0,80),
     ylim=c(0,15),
     las=1, 
     breaks=100)
dev.off()
########################################################################################################################
pdf("H3K56ac.HEG4_empty.500_window.pdf")
hist(H3K56ac.HEG4_empty.500_window$COVERAGE, 
     main="Histogram of H3K56ac data for empty sites in HEG4 with bins .5kb wide", 
     xlab="Coverage", 
     ylab="Frequency, bin = 100", 
     border="blue", 
     col="red",
     xlim=c(0,80),
     ylim=c(0,15),
     las=1, 
     breaks=100)
dev.off()
########################################################################################################################








