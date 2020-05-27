#!/opt/Python/2.7.3/bin/python
import sys
from collections import defaultdict
from numpy import *
import re
import os, gzip
import argparse
from Bio import SeqIO
import random

random.seed(a=None)

def usage():
    test="name"
    message='''
python mPing_sim_9merTSDV2_Random.py --input pictogram/somatic.tsd.matrix --output simulateV2_Random_TSD9mer_somaticMat --replicate 1 --size 100 --use_freq 0

--input: frequency matrix from known mPing insertion
--output output directory where multiple simulated replicates can be stored
--replicate: replicate id, could be 1-100 or a-z
--size: number of mPing sites to simulate
--use_freq: 1 is to use frequency matrix to random pick mPing sites; 0 is to only random pick site from genome sequence.

    '''
    print message

def blackchroms(infile):
    data = defaultdict(str)
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2:
                unit = re.split(r'\t',line)
                data[unit[1] + '_' + unit[2]] = unit[0] 
    return data

def convert_position(start, end, chroms):
    for p in sorted(chroms.keys()):
        p1 = re.split(r'_', p)   
        if start >= int(p1[0]) and end <= int(p1[1]):
            start1 = start - int(p1[0])
            end1   = end - int(p1[0])
            return chroms[p], start1, end1


def blacklist(infile):
    data = defaultdict(int)
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2: 
                unit = re.split(r'\t',line)
                data[unit[0]] = 1
    return data



def fasta(fastafile):
    fastaid = defaultdict(str)
    fastalen = defaultdict(int)
    fh = gzip.open(fastafile,"r")
    for record in SeqIO.parse(fh,"fasta"):
        #print record.id
        seq = str(record.seq.upper())
        fastaid[record.id] = seq
        fastalen[record.id] = len(seq) 
    return fastaid, fastalen

def complement(s): 
    basecomplement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'} 
    letters = list(s) 
    letters = [basecomplement[base] for base in letters] 
    return ''.join(letters)

def revcom(s):
    return complement(s[::-1])


'''
read frquency matirx
>A	C	G	T
0.110	0.335	0.139	0.416
0.247	0.213	0.138	0.402
0.406	0.370	0.144	0.081
0.000	0.000	0.001	0.999
0.478	0.012	0.010	0.500
0.999	0.001	0.000	0.000
0.092	0.149	0.357	0.403
0.397	0.147	0.205	0.252
0.386	0.150	0.346	0.118
'''
def readmatrix(infile):
    data = defaultdict(lambda : defaultdict(lambda: float))
    count = 0
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2 and line.startswith(r'0'): 
                unit = re.split(r'\t',line)
                data['A'][count] = float(unit[0])
                data['C'][count] = float(unit[1])
                data['G'][count] = float(unit[2])
                data['T'][count] = float(unit[3])
                count += 1
    return data

def validdna(dna):
    flag = 1
    bases = ['A', 'T', 'C', 'G'] 
    for b in dna:
        if b not in bases:
            flag = 0
    return flag


def insertion_in_genome(fastaseq, fastalen, matrix, use_freq, blacks, chroms):
    data = []
    mping = 0
    while (mping == 0):
        #chrn = random.randint(1,12)
        #chri = 'Chr' + str(chrn)
        # generate random site
        # make whole genome in one chromosome; it is easier to generate random.
        chri = 'Chr1'
        seq = fastaseq[chri]
        rpos = random.randint(1,fastalen[chri])
        # skip 10 bp at chromosome end, which will cause problems in calculate 9mer frequency
        if blacks.has_key(rpos):
            continue
        # get 9mer of random site
        s    = rpos-1
        e    = s + 9 
        sitep = seq[s:e]
        # skip sites if contain non-ATCG bases
        if not validdna(sitep):
            continue
        # calculate probablity of mPing insertion based on sequencee of plus or minus strand
        siten = revcom(sitep)
        probp = 1.00 #plus strand
        probn = 1.00 #minus strand
        for i in range(0,9):
            basep = sitep[i]
            basen = siten[i]
            probp = probp*matrix[basep][i]
            probn = probn*matrix[basen][i]
        prob = max(probp, probn)
        # use_freq==1 : use max probability to detertmine plus or minus strand the new insertion is going to be
        # use_freq==0 : random plus and minus strand
        site = ''
        strand = ''
        if int(use_freq) == 1:
            site = sitep if probp == prob else siten
            strand = '+' if probp == prob else '-'
        else:
            rn_strand = random.random()
            site = sitep if rn_strand < 0.5 else siten
            strand = '+' if rn_strand < 0.5 else '-'
        # determine whether put this insertion in simuation gff file
        # if use frequency matrix, the chance of outputing this insertion gff depend on prob
        if int(use_freq) == 1:
            rn = random.random()
            if rn <= prob:
                mping = 1
                start = s + 3
                end   = s + 5
                chr1, start1, end1 = convert_position(start, end, chroms)
                data = [chr1, str(start1), str(end1), strand]
        # if not use frequency matrix, output every insertio in gff
        else:
            mping = 1
            start = s + 3
            end   = s + 5
            chr1, start1, end1 = convert_position(start, end, chroms)
            data = [chr1, str(start1), str(end1), strand]
            #print '%s\t%s\t%s' %(chr1, str(start1), str(end1))
    return data


def simulate(fastaseq, fastalen, matrix, outdir, replicate, size, use_freq, blacks, chroms):
    data = defaultdict(list)
    cmd0 = 'mkdir %s' %(outdir)
    os.system(cmd0)
    if 1:
        sufix = '%04d' %(int(replicate))
        filename = '%s/Simulate' %(outdir) + sufix + '.gff'
        tempgff  = '%s/Temp' %(outdir) + sufix + '.gff'
        ofile = open (tempgff, 'w')
        # simulate n=size random mPing insertions
        for i in range(1, int(size)+1):
            print 'Simulating %s mPing site' %(i)
            sys.stdout.flush()
            mping = insertion_in_genome(fastaseq, fastalen, matrix, use_freq, blacks, chroms) 
            print >> ofile, '%s\tSimulate\tmPing\t%s\t%s\t.\t%s\t.\tID=mPing_%s' %(mping[0], mping[1], mping[2], mping[3], i)
        ofile.close()
        # sort gff file by coordinate 
        cmd1 = 'sort -k1,1 -k4,4n %s > %s' %(tempgff, filename)
        cmd2 = 'rm %s' %(tempgff)
        os.system(cmd1)
        os.system(cmd2)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-o', '--output')
    parser.add_argument('-r', '--replicate')
    parser.add_argument('-s', '--size')
    parser.add_argument('-i', '--input')
    parser.add_argument('--use_freq')
    parser.add_argument('-v', dest='verbose', action='store_true')
    args = parser.parse_args()

    if args.output is None:
        args.output = 'simulate_TSD9mer'

    if args.replicate is None:
        args.replicate = '0'

    if args.size is None:
        args.size = 50

    if args.use_freq is None:
        args.use_freq = 0

    try:
        len(args.input) > 0
    except:
        usage()
        sys.exit(2)
   
 
    #ref = '/rhome/cjinfeng/BigData/00.RD/seqlib/MSU_r7.fa'
    ref = 'MSU7.Simulation.Genome.fa.gz'
    fastaseq, fastalen = fasta(ref)
    matrix = readmatrix(args.input)
    blacks = blacklist('MSU7.Simulation.Blacklist')
    chroms = blackchroms('MSU7.Simulation.Chr')
    simulate(fastaseq, fastalen, matrix, args.output, args.replicate, args.size, args.use_freq, blacks, chroms)
    

if __name__ == '__main__':
    main()

