library(gridExtra)
library(ggpubr)
library(ggplot2)
library(tidyr)
library(tidyverse)
library(pheatmap)
library(RColorBrewer)
statdir="coverage/processed_empty"
bgdir="coverage/processed_background"

# tibble will be
# Isolate Mark Chrom Start End Coverage
ks_set = tibble()
Total = tibble()
all_sites = tibble()
for (mark in c("H3K56ac","H3K27me3","ATAC","H3K36me3") ) {
  for (strain in c("A119","EG4","HEG4","A123") ) {
    for ( size in c(1000) ) {
            for ( rep in c("rep1") ) { # lets just use rep1 for now, not sure they should be combined as you did
                # this needs to be a bunch more simulated values so it would be better
                # if we pulled these file names in from a folder
                for (Simulate in c("Simulate0001_100") ) {
                    #"Simulate0001_165",
                    #               "Simulate0001_208", "Simulate0001_310",
                    #               "Simulate0001_99" )) {
                    if ( mark == "ATAC" ) {
                        simfile = sprintf("%s_%s.%s_background_sites.%d.regions.bed.gz",
                                          mark,strain,Simulate,size)
                    } else {
                        simfile = sprintf("%s_%s.%s.%s_background_sites.%d.regions.bed.gz",
                                          mark,strain,rep,Simulate,size)
                    }
                    simfile = file.path(bgdir,mark,simfile)
                    print(simfile)
                    simin <- read_table2(simfile,col_names=FALSE)
                    colnames(simin) <- c("Chr", "Start", "End","Name","Coverage")
                    simin$Mark <- c(mark)
                    simin$Strain <- c(strain)
                    simin$Type  <- c("Sim")
                }

                if ( mark == "ATAC" ) {
                    filename=sprintf("%s_%s.%s_empty.%d.regions.bed.gz",
                                     mark,strain,strain,size)

                } else {
                    filename=sprintf("%s_%s.%s.%s_empty.%d.regions.bed.gz",
                                     mark,strain,rep,strain,size)
                }
                filename = file.path(statdir,mark,filename)
                print(filename)
                obsin <- read_table2(filename,col_names=FALSE)
                colnames(obsin) <- c("Chr", "Start", "End","Name","Coverage")
                obsin$Mark <- c(mark)
                print(mark)
                obsin$Strain <- c(strain)
                obsin$Type  <- c("Obs")
                all_sites = bind_rows(all_sites,obsin)
               
                combo = bind_rows(simin,obsin)
                Total = bind_rows(Total,combo)
                all_sites = bind_rows(all_sites,obsin)
                #glimpse(combo)
                ks.result <- ks.test(simin$Coverage,obsin$Coverage)
                k = tibble(Strain  = strain,
                           Mark    = mark,
                           Window_size = size,
                           Pvalue = ks.result$p.value,
                           Median_Obs = median(obsin$Coverage),
                           Median_Sim = median(simin$Coverage))
                ks_set = bind_rows(ks_set,k)
            }
        }
  }
}

A119 <- subset(all_sites,all_sites$Strain == "A119")
A119$Position <- sprintf("%s_%s",A119$Chr,A119$Start)
A119sub <- select(A119,Position,Chr,Start,Coverage,Mark)
p <- ggplot(A119sub,aes(x=Mark,y=Position,fill=Coverage)) + geom_tile(colour="white",size=0.25) +  scale_y_discrete(expand=c(0,0)) + 
  theme_grey(base_size=8) + theme(
    #bold font for legend text
    legend.text=element_text(face="bold"),
    #set thickness of axis ticks
    axis.ticks=element_line(size=0.4),
    #remove plot background
    plot.background=element_blank(),
    #remove plot border
    panel.border=element_blank()) 

p

