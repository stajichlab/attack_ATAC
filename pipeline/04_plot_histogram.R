#!/usr/bin/env Rscript

setwd("/coverage/processed_empty/ATAC")

# create list of all .csv files in folder 
file_list <- list.files(path=bedgz, pattern="bed.gz$") 

# forloop to loop through dataset
for(file.name in file_list) {
	#read in file
  temp.data <- read.table(file.name)
  #save pdf under this title
  pdf(paste0(gsub("\\.bed\\.gz", "", file.name), ".pdf"))
  # make histogram with column 5
  hist(temp.data$V5,
       main=paste("Histogram from", file.name),
       xlab="Coverage",
       ylab="Frequency, bin = 100",
       border="blue",
       col="green",
       xlim=c(0,15),
       ylim=c(0,60),
       las=1,
       breaks=100)
  dev.off()
}

# remove tempfiles
rm(file.name, temp.data)

####################################################################################

setwd("/coverage/processed_empty/H3K27me3")

# create list of all .csv files in folder 
file_list <- list.files(path=bedgz, pattern="bed.gz$") 

# forloop to loop through dataset
for(file.name in file_list) {
	#read in file
  temp.data <- read.table(file.name)
  #save pdf under this title
  pdf(paste0(gsub("\\.bed\\.gz", "", file.name), ".pdf"))
  # make histogram with column 5
  hist(temp.data$V5,
       main=paste("Histogram from", file.name),
       xlab="Coverage",
       ylab="Frequency, bin = 100",
       border="blue",
       col="red",
       xlim=c(0,80),
       ylim=c(0,15),
       las=1,
       breaks=100)
  dev.off()
}

# remove tempfiles
rm(file.name, temp.data)

####################################################################################
setwd("/coverage/processed_empty/H3K36me3")

# create list of all .csv files in folder 
file_list <- list.files(path=bedgz, pattern="bed.gz$") 

# forloop to loop through dataset
for(file.name in file_list) {
	#read in file
  temp.data <- read.table(file.name)
  #save pdf under this title
  pdf(paste0(gsub("\\.bed\\.gz", "", file.name), ".pdf"))
  # make histogram with column 5
  hist(temp.data$V5,
       main=paste("Histogram from", file.name),
       xlab="Coverage",
       ylab="Frequency, bin = 100",
       border="blue",
       col="yellow",
       xlim=c(0,80),
       ylim=c(0,15),
       las=1,
       breaks=100)
  dev.off()
}

# remove tempfiles
rm(file.name, temp.data)
####################################################################################
setwd("/coverage/processed_empty/H3K56ac")

# create list of all .csv files in folder 
file_list <- list.files(path=bedgz, pattern="bed.gz$") 

# forloop to loop through dataset
for(file.name in file_list) {
	#read in file
  temp.data <- read.table(file.name)
  #save pdf under this title
  pdf(paste0(gsub("\\.bed\\.gz", "", file.name), ".pdf"))
  # make histogram with column 5
  hist(temp.data$V5,
       main=paste("Histogram from", file.name),
       xlab="Coverage",
       ylab="Frequency, bin = 100",
       border="blue",
       col="orange",
       xlim=c(0,80),
       ylim=c(0,15),
       las=1,
       breaks=100)
  dev.off()
}

# remove tempfiles
rm(file.name, temp.data)
####################################################################################
####################################################################################
####################################################################################
setwd("/coverage/processed_windows/ATAC")

# create list of all .csv files in folder 
file_list <- list.files(path=bedgz, pattern="bed.gz$") 

# forloop to loop through dataset
for(file.name in file_list) {
	#read in file
  temp.data <- read.table(file.name)
  #save pdf under this title
  pdf(paste0(gsub("\\.bed\\.gz", "", file.name), ".pdf"))
  # make histogram with column 5
  hist(temp.data$V5,
       main=paste("Histogram from", file.name),
       xlab="Coverage",
       ylab="Frequency, bin = 100",
       border="blue",
       col="green",
       xlim=c(0,15),
       ylim=c(0,60),
       las=1,
       breaks=100)
  dev.off()
}

# remove tempfiles
rm(file.name, temp.data)

####################################################################################

setwd("/coverage/processed_windows/H3K27me3")

# create list of all .csv files in folder 
file_list <- list.files(path=bedgz, pattern="bed.gz$") 

# forloop to loop through dataset
for(file.name in file_list) {
	#read in file
  temp.data <- read.table(file.name)
  #save pdf under this title
  pdf(paste0(gsub("\\.bed\\.gz", "", file.name), ".pdf"))
  # make histogram with column 5
  hist(temp.data$V5,
       main=paste("Histogram from", file.name),
       xlab="Coverage",
       ylab="Frequency, bin = 100",
       border="blue",
       col="red",
       xlim=c(0,80),
       ylim=c(0,15),
       las=1,
       breaks=100)
  dev.off()
}

# remove tempfiles
rm(file.name, temp.data)

####################################################################################
setwd("/coverage/processed_windows/H3K36me3")

# create list of all .csv files in folder 
file_list <- list.files(path=bedgz, pattern="bed.gz$") 

# forloop to loop through dataset
for(file.name in file_list) {
	#read in file
  temp.data <- read.table(file.name)
  #save pdf under this title
  pdf(paste0(gsub("\\.bed\\.gz", "", file.name), ".pdf"))
  # make histogram with column 5
  hist(temp.data$V5,
       main=paste("Histogram from", file.name),
       xlab="Coverage",
       ylab="Frequency, bin = 100",
       border="blue",
       col="yellow",
       xlim=c(0,80),
       ylim=c(0,15),
       las=1,
       breaks=100)
  dev.off()
}

# remove tempfiles
rm(file.name, temp.data)
####################################################################################
setwd("/coverage/processed_windows/H3K56ac")

# create list of all .csv files in folder 
file_list <- list.files(path=bedgz, pattern="bed.gz$") 

# forloop to loop through dataset
for(file.name in file_list) {
	#read in file
  temp.data <- read.table(file.name)
  #save pdf under this title
  pdf(paste0(gsub("\\.bed\\.gz", "", file.name), ".pdf"))
  # make histogram with column 5
  hist(temp.data$V5,
       main=paste("Histogram from", file.name),
       xlab="Coverage",
       ylab="Frequency, bin = 100",
       border="blue",
       col="orange",
       xlim=c(0,80),
       ylim=c(0,15),
       las=1,
       breaks=100)
  dev.off()
}

# remove tempfiles
rm(file.name, temp.data)


#ATAC.A119_empty.1000_window = read.table("ATAC.A119_empty.1000_window.regions.bed.gz",header=F)    
#colnames(ATAC.A119_empty.1000_window) <- c("CHR","START","END","NAME","COVERAGE")
#pdf("ATAC.A119_empty.1000_window.pdf")
#hist(ATAC.A119_empty.1000_window$COVERAGE, 
#     main="Histogram of ATAC data for empty sites in A119 with bins 1kb wide", 
#     xlab="Coverage", 
#     ylab="Frequency, bin = 100", 
#     border="blue", 
#     col="green",
#     xlim=c(0,15),
#     ylim=c(0,60),
#     las=1, 
#     breaks=100)
#dev.off()
########################################################################################################################
#pdf("H3K27me3.A119_empty.500_window.pdf")
#hist(H3K27me3.A119_empty.500_window$COVERAGE, 
#     main="Histogram of H3K27me3 data for empty sites in A119 with bins .5kb wide", 
#     xlab="Coverage", 
#     ylab="Frequency, bin = 100", 
#     border="blue", 
#     col="green",
#     xlim=c(0,80),
#     ylim=c(0,15),
#     las=1, 
#     breaks=100)
#dev.off()