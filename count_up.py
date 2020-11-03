#!/usr/bin/env python3

#%%
# this is a python script template
# this next line will download the file using curl

gff="Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.37.gff3.gz"
fasta="Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.dna.chromosome.Chromosome.fa.gz"

import os,gzip,itertools,csv,re

# this is code which will parse FASTA files
# define what a header looks like in FASTA format
def isheader(line):
    return line[0] == '>'

def aspairs(f):
    seq_id = ''
    sequence = ''
    for header,group in itertools.groupby(f, isheader):
        if header:
            line = next(group)
            seq_id = line[1:].split()[0]
        else:
            sequence = ''.join(line.strip() for line in group)
            yield seq_id, sequence



if not os.path.exists(gff):
    os.system("curl -O ftp://ftp.ensemblgenomes.org/pub/bacteria/release-45/gff3/bacteria_0_collection/escherichia_coli_str_k_12_substr_mg1655/Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.37.gff3.gz")

if not os.path.exists(fasta):
    os.system("curl -O ftp://ftp.ensemblgenomes.org/pub/bacteria/release-45/fasta/bacteria_0_collection/escherichia_coli_str_k_12_substr_mg1655/dna/Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.dna.chromosome.Chromosome.fa.gz")

i = 0
# [0: sequence, 1: source, 2: feature, 3: start, 4: end, 5: score, 6: strand, 7: phase, 8: attributes]
number_of_genes = 0
with gzip.open(gff,"rt") as fh:
    # now add code to process this
    gff = csv.reader(fh,delimiter="\t")
    for row in gff:
        if not row[0].startswith("#"):
            if row[2] == "gene":
                number_of_genes += 1
                # if i is 30:
                    # break
                # print(row[2])
                i += 1

# Bash test: 
# zcat < Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.37.gff3.gz | cut -f 3,3 | grep -v ### | grep -v #! | grep -v ## | grep gene | sort | uniq -c
print('Total number of genes: ' +  str(number_of_genes))
assert number_of_genes == 4142, 'number_of_genes = ' + str(number_of_genes) + ', expected = 4142'

# %%
