#!/usr/bin/env python
# coding: utf-8

# # Analyzing a Genome

import os
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from Bio import Entrez
import pandas as pd


# ## Parse out length

ruddii = SeqIO.read("ruddii.fasta", "fasta")
length = len(ruddii)
print(length)


# ## Parse out GC Content

for seq_record in SeqIO.parse("ruddii.fasta","fasta"):
    g = seq_record.count('G')
    c = seq_record.count('C')
    gc = ((g + c)/len(ruddii.seq))*100
print(gc)


# ## Parse out ATG forward and backward

ruddii = SeqIO.read("ruddii.fasta", "fasta")
atg_forward = ruddii.seq.count('ATG')
print(atg_forward) 

rev_ruddii = ruddii.reverse_complement(id="TESTING")
atg_reverse = rev_ruddii.seq.count('ATG')
print(atg_reverse)


# ## Create data frame

df = pd.DataFrame([length, gc, atg_forward, atg_reverse], index = ['Length_of_Genome', 'GC_content', 'atg_forward', 'atg_reverse'], columns = ['value'])
print(df)

df.to_csv('ruddi.csv')





