#!/usr/bin/env python
# coding: utf-8

# # Fasta Parse

import os
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from Bio import Entrez
import pandas as pd


# ## Parse out first 10 AAs

fasta1 = SeqIO.read("fasta1.fasta", "fasta")
first10_1 = fasta1.seq[0:9]

fasta2 = SeqIO.read("fasta2.fasta", "fasta")
first10_2 = fasta2.seq[0:9]

fasta3 = SeqIO.read("fasta3.fasta", "fasta")
first10_3 = fasta3.seq[0:9]

fasta4 = SeqIO.read("fasta4.fasta", "fasta")
first10_4 = fasta4.seq[0:9]

first10s = [first10_1, first10_2, first10_3, first10_4]
print(first10s)


# ## Parse out lengths


for seq_record in SeqIO.parse("fasta1.fasta", "fasta"):
    len1 = len(seq_record)

for seq_record in SeqIO.parse("fasta2.fasta", "fasta"):
    len2 = len(seq_record)

for seq_record in SeqIO.parse("fasta3.fasta", "fasta"):
    len3 = len(seq_record)

for seq_record in SeqIO.parse("fasta4.fasta", "fasta"):
    len4 = len(seq_record)

lengths = [len1, len2, len3, len4]
print(lengths)


# ## Parse out number of C's

num_c1 = fasta1.seq.count('C')

num_c2 = fasta2.seq.count('C')

num_c3 = fasta3.seq.count('C')

num_c4 = fasta4.seq.count('C')

number_cs = [num_c1, num_c2, num_c3, num_c4]
print(number_cs)


# ## Write into csv using pandas

ids = [fasta1.id, fasta2.id, fasta3.id, fasta4.id]

df = pd.DataFrame(list(zip(*[ids, first10s, lengths, number_cs])), columns = ['ID', 'First_10_AAs', 'length', 'number_cs'])

df.to_csv('protein_info.csv', index = False)

