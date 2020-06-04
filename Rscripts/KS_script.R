# Check for H3K56ac
setwd("/rhome/ysun/bigdata/epigenome/attack_ATAC/coverage/stats/H3K56ac")

H3K56ac.A119.rep1.A119_empty.500_window = read.table("H3K56ac.A119.rep1.A119_empty.500_window.regions.bed.gz",header=F) 
H3K56ac.A119.rep2.A119_empty.500_window = read.table("H3K56ac.A119.rep2.A119_empty.500_window.regions.bed.gz",header=F) 
H3K56ac_A119.rep1.208_background_sites.500 = read.table("H3K56ac_A119.rep1Simulate0001_208_background_sites.500_window.regions.bed.gz",header=F) 
H3K56ac_A119.rep2.208_background_sites.500 = read.table("H3K56ac_A119.rep2Simulate0001_208_background_sites.500_window.regions.bed.gz",header=F) 

H3K56ac.A119.rep1.A119_empty.500_window <- H3K56ac.A119.rep1.A119_empty.500_window$V5
H3K56ac.A119.rep2.A119_empty.500_window  <- H3K56ac.A119.rep2.A119_empty.500_window$V5
H3K56ac_A119.rep1.208_background_sites.500 <- H3K56ac_A119.rep1.208_background_sites.500$V5
H3K56ac_A119.rep2.208_background_sites.500  <- H3K56ac_A119.rep2.208_background_sites.500$V5

H3K56ac.A119_empty.500_window <- c(H3K56ac.A119.rep1.A119_empty.500_window, H3K56ac.A119.rep2.A119_empty.500_window)
H3K56ac_A119.208_background_sites.500 <- c(H3K56ac_A119.rep1.208_background_sites.500, H3K56ac_A119.rep2.208_background_sites.500)

ks.test(H3K56ac_A119.208_background_sites.500 , H3K56ac.A119_empty.500_window)


H3K56ac.A123.rep1.A123_empty.500 = read.table("H3K56ac.A123.rep1.A123_empty.500_window.regions.bed.gz",header=F) 
H3K56ac.A123.rep2.A123_empty.500 = read.table("H3K56ac.A123.rep2.A123_empty.500_window.regions.bed.gz",header=F) 
H3K56ac_A123.rep1.310_background_sites.500 = read.table("H3K56ac_A123.rep1Simulate0001_310_background_sites.500_window.regions.bed.gz",header=F) 
H3K56ac_A123.rep2.310_background_sites.500 = read.table("H3K56ac_A123.rep2Simulate0001_310_background_sites.500_window.regions.bed.gz",header=F) 

H3K56ac.A123.rep1.A123_empty.500 <- H3K56ac.A123.rep1.A123_empty.500$V5
H3K56ac.A123.rep2.A123_empty.500  <- H3K56ac.A123.rep2.A123_empty.500$V5
H3K56ac_A123.rep1.310_background_sites.500 <- H3K56ac_A123.rep1.310_background_sites.500$V5
H3K56ac_A123.rep2.310_background_sites.500  <- H3K56ac_A123.rep2.310_background_sites.500$V5

H3K56ac.A123_empty.500 <- c(H3K56ac.A123.rep1.A123_empty.500, H3K56ac.A123.rep2.A123_empty.500)
H3K56ac_A123.310_background_sites.500 <- c(H3K56ac_A123.rep1.310_background_sites.500, H3K56ac_A123.rep2.310_background_sites.500)

ks.test(H3K56ac_A123.310_background_sites.500 , H3K56ac.A123_empty.500)


H3K56ac.EG4.rep1.EG4_empty.500 = read.table("H3K56ac.EG4.rep1.EG4_empty.500_window.regions.bed.gz",header=F) 
H3K56ac.EG4.rep2.EG4_empty.500 = read.table("H3K56ac.EG4.rep2.EG4_empty.500_window.regions.bed.gz",header=F) 
H3K56ac_EG4.rep1.165_background_sites.500 = read.table("H3K56ac_EG4.rep1Simulate0001_165_background_sites.500_window.regions.bed.gz",header=F) 
H3K56ac_EG4.rep2.165_background_sites.500 = read.table("H3K56ac_EG4.rep2Simulate0001_165_background_sites.500_window.regions.bed.gz",header=F) 

H3K56ac.EG4.rep1.EG4_empty.500 <- H3K56ac.EG4.rep1.EG4_empty.500$V5
H3K56ac.EG4.rep2.EG4_empty.500  <- H3K56ac.EG4.rep2.EG4_empty.500$V5
H3K56ac_EG4.rep1.165_background_sites.500 <- H3K56ac_EG4.rep1.165_background_sites.500$V5
H3K56ac_EG4.rep2.165_background_sites.500 <- H3K56ac_EG4.rep2.165_background_sites.500$V5

H3K56ac.EG4_empty.500 <- c(H3K56ac.EG4.rep1.EG4_empty.500, H3K56ac.EG4.rep2.EG4_empty.500)
H3K56ac_EG4165_background_sites.500  <- c(H3K56ac_EG4.rep1.165_background_sites.500, H3K56ac_EG4.rep2.165_background_sites.500 )

ks.test(H3K56ac_EG4165_background_sites.500, H3K56ac.EG4_empty.500 )


H3K56ac.HEG4.rep1.HEG4_empty.500 = read.table("H3K56ac.HEG4.rep1.HEG4_empty.500_window.regions.bed.gz",header=F) 
H3K56ac.HEG4.rep2.HEG4_empty.500 = read.table("H3K56ac.HEG4.rep2.HEG4_empty.500_window.regions.bed.gz",header=F) 
H3K56ac_HEG4.rep1.99_background_sites.500 = read.table("H3K56ac_HEG4.rep1Simulate0001_99_background_sites.500_window.regions.bed.gz",header=F) 
H3K56ac_HEG4.rep2.99_background_sites.500 = read.table("H3K56ac_HEG4.rep2Simulate0001_99_background_sites.500_window.regions.bed.gz",header=F) 


H3K56ac.HEG4.rep1.HEG4_empty.500 <- H3K56ac.HEG4.rep1.HEG4_empty.500$V5
H3K56ac.HEG4.rep2.HEG4_empty.500  <- H3K56ac.HEG4.rep2.HEG4_empty.500$V5
H3K56ac_HEG4.rep1.99_background_sites.500 <- H3K56ac_HEG4.rep1.99_background_sites.500$V5
H3K56ac_HEG4.rep2.99_background_sites.500 <- H3K56ac_HEG4.rep2.99_background_sites.500$V5

H3K56ac.HEG4_empty.500 <- c(H3K56ac.HEG4.rep1.HEG4_empty.500, H3K56ac.HEG4.rep2.HEG4_empty.500)
H3K56ac_HEG4.99_background_sites.500 <- c(H3K56ac_HEG4.rep1.99_background_sites.500, H3K56ac_HEG4.rep2.99_background_sites.500)

ks.test(H3K56ac_HEG4.99_background_sites.500, H3K56ac.HEG4_empty.500)
