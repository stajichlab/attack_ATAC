#!/usr/bin/env Rscript

library(tools)
library(ggplot2)
library(gridExtra)

# create list of all .csv files in folder
# I WOULD MAKE THIS RELATIVE NOT ABSOLUTE
#infolder <- "/rhome/ysun/bigdata/epigenome/attack_ATAC/"
infolder <- "coverage/processed_empty/H3K56ac"
plotfolder="plots/H3K56ac"

# functions
plot_histo <- function(bedfile) {
  bname = file_path_sans_ext(basename(bedfile),compression=TRUE)
  bname=gsub("\\.region","",bname,perl=TRUE)
  bedin <- read.table(bedfile, col.names = c("Chr", "Start", "End","Name","Coverage"), sep="\t")
  pdfname=sprintf("%s/%shist.pdf",plotfolder,bname)
  pdf(pdfname,width=6)
  hist(bedin$Coverage,
       main=paste("Histogram from", bname),
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

# main code
file_list <- list.files(path=infolder, pattern="bed.gz$",full.names=TRUE)
plots <- lapply(file_list,plot_histo)

# Functions 2
plot_box <- function(bedfile) {
  bname = file_path_sans_ext(basename(bedfile),compression=TRUE)
  bname=gsub("\\.region","",bname,perl=TRUE)
  bedin <- read.table(bedfile, col.names = c("Chr", "Start", "End","Name","Coverage"), sep="\t")
  pdfnamebox= sprintf("%s/%sbox.pdf",plotfolder,bname)
  boxplot = ggplot(bedin, aes(x=factor(Chr),y=Coverage,group=factor(Chr))) +
    geom_boxplot() +
    geom_jitter(color="black", size=0.4, alpha=0.9) +
    theme_bw() +
    theme(
      legend.position="none",
      plot.title = element_text(size=11)
    ) + ggtitle(sprintf("Coverage of %s",bname)) +
    xlab("Chromosome")
  ggnamebox= sprintf("%sbox.pdf", bname)
  ggsave(path = plotfolder, ggnamebox, device = "pdf")
}

# main code
file_list <- list.files(path=infolder, pattern="bed.gz$",full.names=TRUE)
plotsbox <- lapply(file_list,plot_box)





