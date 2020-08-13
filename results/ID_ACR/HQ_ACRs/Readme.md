This is a read me for the modified scripts originally obtained here:https://github.com/schmitzlab/The-prevalence-evolution-and-chromatin-signatures-of-plant-regulatory-elements

Goal: obtain high quality ACRs

Set up: This is a 4 part script. For each of the following parts, I explain the goal of the script and the modifications I made to Bob's original script so it will work for my purposes.

01_calculate_density_mod.sh

Purpose: Identify raw ACRs with MACS2 and calculate the density.
(1) ACRs called with MACS2 were split into 50 bp windows with 25 bp steps
(2) the Tn5 integration frequency in each window was calculated and normalized to the average frequency in the total genome

Modifications: 
1. Input file removed
2. Included the modules
3. MACS2 script was modified to remove input sample file which is just a control and not necessary
4. Since only 1 file is used to generate the .Tn5 file, the name is not edited
5. A dictionary was generated because it was not possible to get total0 to produce the cat $reads |wc -l

# run using
sbatch ./01_calculate_density_mod.sh A119_01 50 25 430000000 /rhome/ysun/bigdata/epigenome/attack_ATAC/results/ID_ACR/HQ_ACRs/A119samplefile.txt
sbatch ./01_calculate_density_mod.sh A123_01 50 25 430000000 /rhome/ysun/bigdata/epigenome/attack_ATAC/results/ID_ACR/HQ_ACRs/A123samplefile.txt
sbatch ./01_calculate_density_mod.sh EG4_01 50 25 430000000 /rhome/ysun/bigdata/epigenome/attack_ATAC/results/ID_ACR/HQ_ACRs/EG4samplefile.txt
sbatch ./01_calculate_density_mod.sh HEG4_01 50 25 430000000 /rhome/ysun/bigdata/epigenome/attack_ATAC/results/ID_ACR/HQ_ACRs/HEG4samplefile.txt
sbatch ./01_calculate_density_mod.sh NP_01 50 25 430000000 /rhome/ysun/bigdata/epigenome/attack_ATAC/results/ID_ACR/HQ_ACRs/NPsamplefile.txt



02_choose_cutoff.mod.sh
Purpose: Choose the density cutoff and filter the ACRs by removing regions similar with Organelle genome(s).
(3) windows passing the integration frequency cut-off were merged together with 150 bp gaps
(4) small regions with only one window were then filtered with ‘length > 50 bp’
(5) regions aligning to the mitochondrial or chloroplast genome from NCBI Organelle Genome Resources were also removed. The sites within ACRs with the highest Tn5 integration frequency were defined as summits.

Modifications: 
1. Added CPU component
2. Added component to genrate a database
3. removed database input from original script
# commented out
export PATH=/rhome/ysun/bigdata/epigenome/attack_ATAC/results/ID_ACR/HQ_ACRs/ncbi-blast-2.8.1
4. redownloaded Oryza_sativa.IRGSP-1.0.30.dna.toplevel.modified.fa and it is now Oryza_sativa.IRGSP-1.0.dna.toplevel.modified.fa
# wget new fa file
ftp://ftp.ensemblgenomes.org/pub/plants/release-47/fasta/oryza_sativa/dna/Oryza_sativa.IRGSP-1.0.dna.toplevel.fa.gz
    # to include mitochondra etc. 
5. Modified .fa file
- add Chr, remove everything after space
# replace headers with proper headers
awk '/^>/ {$0=$1} 1' Oryza_sativa.IRGSP-1.0.dna.toplevel.fa | awk -F "::" '{if($1~">"){gsub(">","Chr");print ">"$1} else {print $0}}' > Oryza_sativa.IRGSP-1.0.dna.toplevel_mod.fa
- reformat .fa file
# reformat fasta file
module load hmmer/3
esl-reformat fasta Oryza_sativa.IRGSP-1.0.dna.toplevel_mod.fa > Oryza_sativa.IRGSP-1.0.dna.toplevel_reformat.fa
rm Oryza_sativa.IRGSP-1.0.dna.toplevel_mod.fa
rm Oryza_sativa.IRGSP-1.0.dna.toplevel_mod.fa.fai
6. replace $input.$density.$gap.bedo with $density.$cutoff.$gap.bedo

# run using
sh 02_choose_cutoff.mod.sh A119.density 5 150 50 1 A119
sh 02_choose_cutoff.mod.sh A123.density 5 150 50 1 A123
sh 02_choose_cutoff.mod.sh EG4.density 5 150 50 1 EG4
sh 02_choose_cutoff.mod.sh HEG4.density 5 150 50 1 HEG4
sh 02_choose_cutoff.mod.sh NP.density 5 150 50 1 NP


03_prepareTag.mod.sh
Purpose: To call the summits, a tag directory containing Tn5 integration sites were generate.

# run using
sbatch 03_prepareTag.mod.sh  A119_03 Oryza_sativa.IRGSP-1.0.chrom_info.txt no 430000000 A119TN5.txt
sbatch 03_prepareTag.mod.sh  A123_03 Oryza_sativa.IRGSP-1.0.chrom_info.txt no 430000000 A123TN5.txt
sbatch 03_prepareTag.mod.sh  EG4_03 Oryza_sativa.IRGSP-1.0.chrom_info.txt no 430000000 EG4TN5.txt
sbatch 03_prepareTag.mod.sh  HEG4_03 Oryza_sativa.IRGSP-1.0.chrom_info.txt no 430000000 HEG4TN5.txt
sbatch 03_prepareTag.mod.sh  NP_03 Oryza_sativa.IRGSP-1.0.chrom_info.txt no 430000000 NPTN5.txt


Modifications: 
- removed $X from '$[total$x]'
- removed * from size'*$2*/ so now it's size'*$2/

04_call_summits.mod.sh
Purpose: Call the summits. The position with highest Tn5 Integration sites were considered as the summit(s).

# run using
sbatch 04_call_summits.mod.sh A119_04 /rhome/ysun/bigdata/epigenome/attack_ATAC/results/ID_ACR/HQ_ACRs/A119cutoff /rhome/ysun/bigdata/epigenome/attack_ATAC/results/ID_ACR/HQ_ACRs/A119.density /rhome/ysun/bigdata/epigenome/attack_ATAC/results/ID_ACR/HQ_ACRs/A119_03 /rhome/ysun/bigdata/epigenome/attack_ATAC/results/ID_ACR/HQ_ACRs/Oryza_sativa.IRGSP-1.0.chrom_info.txt no
sbatch 04_call_summits.mod.sh A123_04 /rhome/ysun/bigdata/epigenome/attack_ATAC/results/ID_ACR/HQ_ACRs/A123cutoff /rhome/ysun/bigdata/epigenome/attack_ATAC/results/ID_ACR/HQ_ACRs/A123.density /rhome/ysun/bigdata/epigenome/attack_ATAC/results/ID_ACR/HQ_ACRs/A123_03 /rhome/ysun/bigdata/epigenome/attack_ATAC/results/ID_ACR/HQ_ACRs/Oryza_sativa.IRGSP-1.0.chrom_info.txt no
sbatch 04_call_summits.mod.sh EG4_04 /rhome/ysun/bigdata/epigenome/attack_ATAC/results/ID_ACR/HQ_ACRs/EG4cutoff /rhome/ysun/bigdata/epigenome/attack_ATAC/results/ID_ACR/HQ_ACRs/EG4.density /rhome/ysun/bigdata/epigenome/attack_ATAC/results/ID_ACR/HQ_ACRs/EG4_03 /rhome/ysun/bigdata/epigenome/attack_ATAC/results/ID_ACR/HQ_ACRs/Oryza_sativa.IRGSP-1.0.chrom_info.txt no
sbatch 04_call_summits.mod.sh HEG4_04 /rhome/ysun/bigdata/epigenome/attack_ATAC/results/ID_ACR/HQ_ACRs/HEG4cutoff /rhome/ysun/bigdata/epigenome/attack_ATAC/results/ID_ACR/HQ_ACRs/HEG4.density /rhome/ysun/bigdata/epigenome/attack_ATAC/results/ID_ACR/HQ_ACRs/HEG4_03 /rhome/ysun/bigdata/epigenome/attack_ATAC/results/ID_ACR/HQ_ACRs/Oryza_sativa.IRGSP-1.0.chrom_info.txt no
sbatch 04_call_summits.mod.sh NP_04 /rhome/ysun/bigdata/epigenome/attack_ATAC/results/ID_ACR/HQ_ACRs/NPcutoff /rhome/ysun/bigdata/epigenome/attack_ATAC/results/ID_ACR/HQ_ACRs/NP.density /rhome/ysun/bigdata/epigenome/attack_ATAC/results/ID_ACR/HQ_ACRs/NP_03 /rhome/ysun/bigdata/epigenome/attack_ATAC/results/ID_ACR/HQ_ACRs/Oryza_sativa.IRGSP-1.0.chrom_info.txt no
