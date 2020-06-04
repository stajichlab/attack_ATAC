library(gridExtra)
library(ggpubr)
library(ggplot2)
library(tidyverse)



statdir="coverage/processed_empty"
bgdir="coverage/processed_background"

# tibble will be
# Isolate Mark Chrom Start End Coverage
ks_set = tibble()
Total = tibble()
for (mark in c("H3K56ac","H3K27me3","ATAC","H3K36me3") ) {
  for (strain in c("A119","EG4","HEG4","A123") ) {
    for ( size in c(500, 1000) ) {
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
                obsin$Strain <- c(strain)
                obsin$Type  <- c("Obs")
                combo = bind_rows(simin,obsin)
                Total = bind_rows(Total,combo)
                #glimpse(combo)
                ks.result <- ks.test(simin$Coverage,obsin$Coverage)
                k = tibble(Compare  = sprintf("%s-%s-%d",strain,mark,size),
                                          Pvalue = ks.result$p.value,
                                          Median_Obs = median(obsin$Coverage),
                                          Median_Sim = median(simin$Coverage))

                ks_set = bind_rows(ks_set,k)
#                boxplot <- ggplot(combo, aes(x=Strain,y=Coverage,color=Strain)) +
#                  geom_boxplot() +
#                  geom_jitter(color="black", size=0.4, alpha=0.9) +
#                  theme_bw() +
#                  theme(
#                    legend.position="none",
#                    plot.title = element_text(size=11)
#                  ) + scale_fill_brewer(palette="Set1") +
#                 ggtitle(sprintf("Coverage of %s %s",mark,strain)) + xlab("Strain") +
#                    ylab(sprintf("%s Coverage",mark))
            }
        }
  }
}
write_csv(ks_set,"ks_compare_Obs-vs-Sim.csv")
outplotfile <- "plots/emptysitemPing_compare_epi.pdf"
boxplot <- ggplot(Total, aes(x=Mark,y=Coverage,color=Mark)) +
  geom_boxplot() +
  geom_jitter(color="black", size=0.4, alpha=0.9) +
  theme_bw() +
  theme(
    legend.position="none",
    plot.title = element_text(size=11)) + scale_fill_brewer(palette="Set1") +
  ylab("Coverage") + facet_grid(cols=vars(Type),rows=vars(Strain))
ggsave(outplotfile,boxplot,width=10,height=8)