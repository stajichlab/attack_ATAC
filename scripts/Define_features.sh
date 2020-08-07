#!/bin/bash
#SBATCH -p short

# runs script in /rhome/ysun/bigdata/epigenome/attack_ATAC/results/Closest
# Purpose: this takes gff3 file and mPing positions and uses both to define features for mPing positions
# output is a bed file that has mPing features defined
# for output 530 mPing fall into these regions
# some double up because some genes are close together in rice. 

#order data
module load bedtools
module load mosdepth
module unload perl
module load parallel

# pull out column information if it's for CDS (+)
# match with mRNA
# match with (+) or (-)
awk '$3 == "mRNA" {print $1 "\t" $4 "\t" $5 "\t" $3 "\t" $7 "\t" $9}' Oryza_sativa.IRGSP-1.0.47.chr.gff3| awk '$5 == "+" {print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6}' | awk '{gsub(".*=","",$6)}1' > forwardstrand.gtf 
# pull out column information if it's for CDS (+)
#awk '$3 == "mRNA" {print $1 "\t" $4 "\t" $5 "\t" $3 "\t" $7 "\t" $9}' Oryza_sativa.IRGSP-1.0.47.chr.gff3 | awk '$5 == "-" {print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6}' | awk '{gsub(".*:","",$6)}1' | awk '{gsub(";.*","",$6)}1' > reversestrand.gtf
awk '$3 == "mRNA" {print $1 "\t" $4 "\t" $5 "\t" $3 "\t" $7 "\t" $9}' Oryza_sativa.IRGSP-1.0.47.chr.gff3 | awk '$5 == "-" {print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6}' | awk '{gsub(".*=","",$6)}1' > reversestrand.gtf

# assign 5' promoters to genes in file
awk -v proximal=500 -v distal=500 -v upstream=500 -v farupstream=2000 '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $2-proximal "\t" $2 "\t" "five-proximal" "\t" $2-proximal-distal "\t" $2-proximal "\t" "five-distal" "\t" $2-proximal-distal-upstream "\t" $2-proximal-distal "\t" "five-upstream" "\t" $2-proximal-distal-upstream-farupstream "\t" $2-proximal-distal-upstream "\t" "five-farupstream"}' forwardstrand.gtf > annotatedforwardstrand.gtf
awk -v proximal=500 -v distal=500 -v downstream=500 -v farupstream=2000 '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $3 "\t" $3+proximal "\t" "five-proximal" "\t" $3+proximal "\t" $3+proximal+distal "\t" "five-distal" "\t" $3+proximal+distal "\t" $3+proximal+distal+downstream "\t" "five-upstream" "\t" $2-proximal-distal-upstream-farupstream "\t" $2-proximal-distal-upstream "\t" "five-farupstream"}' reversestrand.gtf > annotatedreversestrand.gtf

#forwardstrand.gtf
cat forwardstrand.gtf | awk '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6}' >> intermediatefinalForwardstrand.gtf
cat annotatedforwardstrand.gtf | awk '{print $1 "\t" $7 "\t" $8 "\t" $9 "\t" $5 "\t" $6}'>> intermediatefinalForwardstrand.gtf
cat annotatedforwardstrand.gtf | awk '{print $1 "\t" $10 "\t" $11 "\t" $12 "\t" $5 "\t" $6}' >> intermediatefinalForwardstrand.gtf
cat annotatedforwardstrand.gtf | awk '{print $1 "\t" $13 "\t" $14 "\t" $15 "\t" $5 "\t" $6}' >> intermediatefinalForwardstrand.gtf
cat annotatedforwardstrand.gtf | awk '{print $1 "\t" $16 "\t" $17 "\t" $18 "\t" $5 "\t" $6}' >> intermediatefinalForwardstrand.gtf
sort -k1,1 -k2,2n intermediatefinalForwardstrand.gtf| tail -n +2 > finalForwardstrand.gtf
# mannually removed this from dataset
11	-1451	549	five-farupstream	+	Os11t0100150-00
12	-818	1182	five-farupstream	+	Os12t0100100-01
12	-773	1227	five-farupstream	+	Os12t0100100-02
12	-740	1260	five-farupstream	+	Os12t0100100-03
2	-3108	-1108	five-farupstream	+	Os02t0100050-00
2	-1108	-608	five-upstream	+	Os02t0100050-00
2	-608	-108	five-distal	+	Os02t0100050-00
2	-108	392	five-proximal	+	Os02t0100050-00
#reversestrand.gtf
cat reversestrand.gtf | awk '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6}' >> intermediatefinalReversestrand.gtf
cat annotatedreversestrand.gtf | awk '{print $1 "\t" $7 "\t" $8 "\t" $9 "\t" $5 "\t" $6}'>> intermediatefinalReversestrand.gtf
cat annotatedreversestrand.gtf | awk '{print $1 "\t" $10 "\t" $11 "\t" $12 "\t" $5 "\t" $6}' >> intermediatefinalReversestrand.gtf
cat annotatedreversestrand.gtf | awk '{print $1 "\t" $13 "\t" $14 "\t" $15 "\t" $5 "\t" $6}' >> intermediatefinalReversestrand.gtf
cat annotatedreversestrand.gtf | awk '{print $1 "\t" $16 "\t" $17 "\t" $18 "\t" $5 "\t" $6}' >> intermediatefinalReversestrand.gtf
sort -k1,1 -k2,2n intermediatefinalReversestrand.gtf > finalReversestrand.gtf
# manually removed this from dataset
11	-1820	180	five-farupstream	-	Os11t0100100-01
Pt	-2918	-918	five-farupstream	-	transcript-psbA

# assign 3' UTR to genes in file
awk -v proximal=499 -v distal=499 -v upstream=499 '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $2-proximal "\t" $2 "\t" "three-proximal" "\t" $2-proximal-distal "\t" $2-proximal "\t" "three-distal" "\t" $2-proximal-distal-upstream "\t" $2-proximal-distal "\t" "three-downstream"}' reversestrand.gtf > annotatedforwardstrand3prime.gtf
awk -v proximal=499 -v distal=499 -v downstream=499 '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $3 "\t" $3+proximal "\t" "three-proximal" "\t" $3+proximal "\t" $3+proximal+distal "\t" "three-distal" "\t" $3+proximal+distal "\t" $3+proximal+distal+downstream "\t" "three-downstream"}' forwardstrand.gtf > annotatedreversestrand3prime.gtf
#forwardstrand.gtf
cat annotatedforwardstrand3prime.gtf | awk '{print $1 "\t" $7 "\t" $8 "\t" $9 "\t" $5 "\t" $6}'>> intermediatefinalForwardstrand3prime.gtf
cat annotatedforwardstrand3prime.gtf | awk '{print $1 "\t" $10 "\t" $11 "\t" $12 "\t" $5 "\t" $6}' >> intermediatefinalForwardstrand3prime.gtf
cat annotatedforwardstrand3prime.gtf | awk '{print $1 "\t" $13 "\t" $14 "\t" $15 "\t" $5 "\t" $6}' >> intermediatefinalForwardstrand3prime.gtf
sort -k1,1 -k2,2n intermediatefinalForwardstrand3prime.gtf| tail -n +2 > finalForwardstrand3prime.gtf
# mannually removed this from dataset
11	-323	178	three-downstream	-	Os11t0100100-01
Pt	-1421	-920	three-downstream	-	transcript-psbA
Pt	-920	-419	three-distal	-	transcript-psbA
Pt	-419	82	three-proximal	-	transcript-psbA
#reversestrand.gtf
cat annotatedreversestrand3prime.gtf | awk '{print $1 "\t" $7 "\t" $8 "\t" $9 "\t" $5 "\t" $6}'>> intermediatefinalReversestrand3prime.gtf
cat annotatedreversestrand3prime.gtf | awk '{print $1 "\t" $10 "\t" $11 "\t" $12 "\t" $5 "\t" $6}' >> intermediatefinalReversestrand3prime.gtf
cat annotatedreversestrand3prime.gtf | awk '{print $1 "\t" $13 "\t" $14 "\t" $15 "\t" $5 "\t" $6}' >> intermediatefinalReversestrand3prime.gtf
sort -k1,1 -k2,2n intermediatefinalReversestrand3prime.gtf > finalReversestrand3prime.gtf


# run bedtools

# for forward strand, there are 392 cases
bedtools intersect -a finalForwardstrand.gtf -b all_empty_sites.bed -wo | awk '{print $7 "\t" $8 "\t" $9 "\t" $10 "\t" $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6}' >> presortmPingPromoter.bed
# for reverse strand, there are 359 cases
bedtools intersect -a finalReversestrand.gtf -b all_empty_sites.bed -wo | awk '{print $7 "\t" $8 "\t" $9 "\t" $10 "\t" $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6}' >> presortmPingPromoter.bed

# do the 3prime
bedtools intersect -a finalForwardstrand3prime.gtf -b all_empty_sites.bed -wo | awk '{print $7 "\t" $8 "\t" $9 "\t" $10 "\t" $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6}' >> presortmPingPromoter.bed
# do the 3prime
bedtools intersect -a finalReversestrand3prime.gtf -b all_empty_sites.bed -wo | awk '{print $7 "\t" $8 "\t" $9 "\t" $10 "\t" $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6}' >> presortmPingPromoter.bed

sort -k1,1 -k2,2n presortmPingPromoter.bed > mPingPromoter.bed

### There are 926 positions in this file, covering 530 mPing sites
