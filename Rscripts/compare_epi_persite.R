library(gridExtra)
library(ggpubr)
library(ggplot2)
library(tidyr)
library(tidyverse)
library(pheatmap)
library(RColorBrewer)

statdir="coverage/processed_empty"
bgdir="coverage/processed_background"

Simulations =  c("Simulate0001_100")
#                                        #"Simulate0001_165",
#               "Simulate0001_208", "Simulate0001_310",
#               "Simulate0001_99" )) {

parse_bedfile <- function(bedfile) {
    bed <- read_table2(bedfile,
                       col_names=c("Chr", "Start", "End","Name","Coverage"),
                       col_types = cols_only(
                           Chr   = "f",
                           Start = "i",
                           End   = "i",
                           Name  = "_",
                           Coverage = "d"))
    return(bed)
}

# tibble will be
# Isolate Mark Chrom Start End Coverage
#for ( size in c( 500 ) ) {
size=1000
all_sites = tibble()
for (strain in c("A119","EG4","HEG4","A123") ) {
           # this is here because there are no reps for ATAC just one rep
           # first parse ATAC from simulated sites
    for (Simulate in c("Simulate0001_100") ) {
        simfile = sprintf("%s_%s.%s_background_sites.%d.regions.bed.gz",
                          "ATAC",strain,Simulate,size)
        simfile = file.path(bgdir,"ATAC",simfile)
        simin <- parse_bedfile(simfile)
        simin$Strain <- c(strain)
        simin$Type  <- "Sim"
        simin$Mark  <- "ATAC"
        simin$Rep <- "rep1"
        all_sites = bind_rows(all_sites,simin)
    }
           # then parse the ATAC from observed sites for this strain
    filename=sprintf("%s_%s.%s_empty.%d.regions.bed.gz",
                     "ATAC",strain,strain,size)
          
    filename = file.path(statdir,"ATAC",filename)

    obsin <- parse_bedfile(filename)
    obsin$Strain <- c(strain)
    obsin$Type  <- "Obs"
    obsin$Mark  <- "ATAC"
    obsin$Rep   <- "rep1"
    all_sites = bind_rows(all_sites,obsin)
    
    for (mark in c("H3K56ac","H3K27me3","H3K36me3") ){
        for ( rep in c("rep1","rep2") ) { 
            for (Simulate in Simulations ) {
                simfile = sprintf("%s_%s.%s.%s_background_sites.%d.regions.bed.gz",
                                  mark,strain,rep,Simulate,size)
                
                simfile = file.path(bgdir,mark,simfile)
                                        # print(simfile)
                simin <- parse_bedfile(simfile)
                simin$Strain <- c(strain)
                simin$Type  <- c("Sim")
                simin$Mark  <- mark
                simin$Rep <- rep
                all_sites = bind_rows(all_sites,simin)
            }
            
            filename=sprintf("%s_%s.%s.%s_empty.%d.regions.bed.gz",
                             mark,strain,rep,strain,size)
            filename = file.path(statdir,mark,filename)
                                        #print(filename)
            obsin <- parse_bedfile(filename)
            obsin$Strain <- c(strain)
            obsin$Type  <- c("Obs")
            obsin$Mark  <- mark
            obsin$Rep   <- rep
            
            all_sites = bind_rows(all_sites,obsin)
        }
    }
}
marksplot <- all_sites %>% mutate(Position = paste0(Type,"_",Chr, ":", Start)) %>% 
    select(-c(Chr,Start,End,Type)) %>%
    pivot_wider(names_from = c(Mark,Rep), values_from = c(Coverage)) %>%
    as.data.frame() 

#head(marksplot)
#glimpse(marksplot)
pdf("epi_heatmap.pdf",height=20)
for (strain in unique(marksplot$Strain)) {
  s <- subset(marksplot,marksplot$Strain == strain)
  row.names(s) <- s$Position
  s <- s %>% select(-c(Position,Strain))
  pheatmap(s, scale="column",fontsize_row=4,main=sprintf("%s %d windows",strain,size),
         cluster_rows = TRUE,
         cluster_cols = TRUE)
}
