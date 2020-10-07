
import sys
import numpy as np

from Bio import SeqIO, AlignIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

from Bio.Align.Applications import ClustalwCommandline

param = int(sys.argv[1])

spike = list(SeqIO.parse('s_gene.fasta', 'fasta'))[0]

records = list(SeqIO.parse('group_90_spike.fasta', 'fasta'))

result = []

count = 0

for i in range(16):
    tmp = records[(param * 320 + i * 20):(param * 320 + i * 20 + 20)]
    tmp.append(spike)
    SeqIO.write(tmp, 'tmp_' + str(param) +'.fasta', 'fasta')

    align_cmd = ClustalwCommandline('clustalw2', infile = 'tmp_' + str(param) +'.fasta')
    align_stdout, align_stderr = align_cmd()

    align = AlignIO.read('tmp_' + str(param) +'.aln', 'clustal')
    align_annot = align.column_annotations['clustal_consensus']

    for t in align:
        if t.id == 'NC_045512.2_21563-25384':
            break

    tmp = list(t.seq)
    align_seq = list(t.seq)
    for i in range(len(align_annot)):
        if align_annot[i] != '*':
            tmp[i] = ' '
            align_seq[i] = '-'
    result.append(SeqRecord(Seq(''.join(align_seq)), 'annotated-{:02d}-{:06d}'.format(param, count), '', ''))
    count += 1
    
SeqIO.write(result, 'annotated-{:02d}.fasta'.format(param), 'fasta')