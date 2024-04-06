#!/usr/bin/env python
# coding: utf-8

# # Exercise 1: Genbank Parse

# ### Import Modules from BioPython in pandas_practice environment

import os
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from Bio import Entrez
import pandas as pd


# ### Open genbank files for reading

gb1 = "NZ_BHVZ01000001.gb"
gb_rec1 = SeqIO.read(open(gb1,"r"), "genbank")

gb2 = "NZ_CAJTFZ010000019.gb"
gb_rec2 = SeqIO.read(open(gb2,"r"), "genbank")

gb3 = "NZ_CALPCP010000001.gb"
gb_rec3 = SeqIO.read(open(gb3,"r"), "genbank")

gb4 = "NZ_CALPCY010000130.gb"
gb_rec4 = SeqIO.read(open(gb4,"r"), "genbank")

gb5 = "NZ_SRYA01000017.gb"
gb_rec5 = SeqIO.read(open(gb5,"r"), "genbank")

all_recs = [gb_rec1, gb_rec2, gb_rec3, gb_rec4, gb_rec5]


# ### Parsing out the source, taxonomy lists, number of features, and ids for all records


for rec in all_recs:
    source = rec.annotations['source']
    split = source.split(" ")
    species = split[1]
    taxonomy = rec.annotations['taxonomy']
    family = taxonomy[-2]
    genus = taxonomy[-1]
    num_features = len(rec.features)
    accession = rec.id
    print([accession, family, genus, species, num_features, source])


# ### Take lists and make a dataframe, convert to csv file (use pandas). 

data = [['Anaerotignaceae', 'Bacteroidaceae', 'Eubacteriales', 'Lachnospiraceae', 'Lachnospiraceae'], ['Anaerotignum', 'Bacteroides', 'Flintibacter', 'Otoolea', 'Petralouisia'], ['faecicola', 'acidifaciens', 'muris', 'muris', 'muris'], [1, 143, 903, 3, 161], ['Anaerotignum faecicola', 'Bacteroides acidifaciens', 'Flintibacter muris',  'Otoolea muris', 'Petralouisia muris']]

df = pd.DataFrame(data, columns=['family', 'genus', 'species', 'num_features', 'source'], index=['NZ_BHVZ01000001.1', 'NZ_CAJTFZ010000019.1', 'NZ_CALPCP010000001.1', 'NZ_CALPCY010000130.1', 'NZ_SRYA01000017.1'])
print(df)

df.to_csv('genbank_parse.csv')

