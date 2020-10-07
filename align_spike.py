
import numpy as np

from Bio import SeqIO, AlignIO
from Bio.SeqRecord import SeqRecord

from Bio.Align.Applications import ClustalwCommandline

spike_records = list(SeqIO.parse('s_gene.fasta', 'fasta'))
std_seq = spike_records[0]
records = list(SeqIO.parse('all.fasta', 'fasta'))

counter = 0
group_90 = []
group_80 = []

for record in records[7500:]:
#     if len(record) < 26000:
#         counter += 1
#         continue
#     else:
#         len_list.append(len(record))
#         mito_frag = record.seq
    mito_frags=[]
    mito_frags.append(std_seq)
    record.id = 'covid-19-gene-seq-{:06d}'.format(counter)
    mito_frags.append(record) # SeqRecord(mito_frag, record.id, '', ''))
    
    SeqIO.write(mito_frags, 'tmp_3.fasta', 'fasta')
    align_cmd = ClustalwCommandline('clustalw2', infile = 'tmp_3.fasta')
    align_stdout, align_stderr = align_cmd()
    counter += 1
    
    score = 0
    t = align_stdout.splitlines()
    for s in t:
        if s.startswith('Sequences ('):
            score = int(s.split()[-1]) 
            
    if score == 100:
        continue
    align = AlignIO.read('tmp_3.aln', 'clustal')
    for a in align:
        if a.id == record.id:
            break
    annot = align.column_annotations['clustal_consensus']
    if score > 90:        
        group_90.append(SeqRecord(a.seq[annot.find('*') : (annot.rfind('*') + 1)], record.id, '', ''))
    else:
        group_80.append(SeqRecord(a.seq[annot.find('*') : (annot.rfind('*') + 1)], record.id, '', ''))
    
    if counter % 100 == 0:
        SeqIO.write(group_90, 'group_90_spike_3.fasta', 'fasta')
        SeqIO.write(group_80, 'group_80_spike_3.fasta', 'fasta')