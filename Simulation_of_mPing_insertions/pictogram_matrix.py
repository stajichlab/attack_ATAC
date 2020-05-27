from Bio.Seq import Seq
from Bio import motifs
from Bio import SeqIO
from Bio.Alphabet import IUPAC
from collections import defaultdict
import numpy as np
import sys

matrix = defaultdict(lambda : defaultdict(lambda : int()))
#for seq_record in SeqIO.parse("mping_flank_simulateV2_Chromatin0.99_Random.tsd.fa", 'fasta'):
#../pictogram/strain.tsd.fa
for seq_record in SeqIO.parse(sys.argv[1], 'fasta'):
    string = str(seq_record.seq)
    for i in range(0, len(string)):
        matrix[i][string[i]] += 1


print '>A\tC\tG\tT'
for rank in sorted(matrix.keys(), key=int):
    frq = []
    for base in ['A', 'C', 'G', 'T']:
        frq.append('%0.3f' %(float(matrix[rank][base])/np.sum(matrix[rank].values())))
        #print np.sum(matrix[rank].values())
    print '\t'.join(map(str, frq))
