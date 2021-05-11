import pandas as pd
from Bio import SeqIO

# read raw data
y2h_union = pd.read_csv('data/raw/yuri/Y2H_union.txt',header=None,names=['gene1','gene2'], sep='\t')
y2h_union.to_csv('data/interim/yuri/Y2H_union.csv',index=False)
ccsb_yi2 = pd.read_csv('data/raw/yuri/CCSB-YI2_preliminary_release.tsv',header=None,names=['gene1','gene2'], sep='\ ')
ccsb_yi2.to_csv('data/interim/yuri/CCSB-YI2_preliminary_release.csv',index=False)

y2h_combined = pd.concat([y2h_union,ccsb_yi2])
y2h_combined.to_csv('data/interim/yuri/yuri_combined.csv',index=False)

#yi2_space = SeqIO.parse('data/raw/yuri/YI2_space.fa','fasta')
#orf_names = [a.name.split('|')[0] for a in yi2_space]


