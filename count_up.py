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
    for header, group in itertools.groupby(f, isheader):
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
ENABLE_TESTS = True #False
# [0: sequence, 1: source, 2: feature, 3: start, 4: end, 5: score, 6: strand, 7: phase, 8: attributes]
number_of_genes = 0
total_length_of_genes = 0
with gzip.open(gff,"rt") as fh:
    # now add code to process this
    gff = csv.reader(fh,delimiter="\t")
    for row in gff:
        if not row[0].startswith("#"):
            if not ENABLE_TESTS and i == 30:
                break
            if row[2] == "gene":
                number_of_genes += 1
                total_length_of_genes += int(row[4]) - int(row[3])
            i += 1

with gzip.open(fasta,"rt") as fh:
    seqs = dict(aspairs(fh))
    first_codon = {}
    last_codon  = {}
    sequence_count = 0
    strand_count = {}
    i = 0

    for seqname in seqs:
        total_length_of_genome = len(seqs[seqname])

answers = [
    '2. Total number of genes: {:,}'.format(number_of_genes),
    '3. Total length of genes: {:,}'.format(total_length_of_genes),
    '4. Total length of genome: {:,}'.format(total_length_of_genome),
    '5. Perecent of genome which is coding: {:.2%}'.format(total_length_of_genes / total_length_of_genome)
]

for answer in answers:
    print(answer)

if ENABLE_TESTS:
    # Bash test: zcat < Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.37.gff3.gz | cut -f 3,3 | grep -v ### | grep -v #! | grep -v ## | grep gene | sort | uniq -c
    assert answers[0] == '2. Total number of genes: 4,142', 'Answer 1 = ' + answers[0] + ', expected = 2. Total number of genes: 4,142'
    # Bash test: zcat < Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.37.gff3.gz | awk '!/^#+/' |  awk '$3 ~ /^gene$/ {print $0}' | awk '{sum += $5-$4} END {print sum}'
    assert answers[1] == '3. Total length of genes: 3,944,686', 'Answer 1 = ' + answers[1] + ', expected = 3. Total length of genes: 3,944,686'
    # Bash test: zcat < Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.dna.chromosome.Chromosome.fa.gz | tail -n 77361 | tr -d '\n' | wc -c
    assert answers[2] == '4. Total length of genome: 4,641,652', 'Answer 1 = ' + answers[2] + ', expected = 4. Total length of genome: 4,641,652'
    assert answers[3] == '5. Perecent of genome which is coding: 84.98%', 'Answer 1 = ' + answers[3] + ', expected = 5. Perecent of genome which is coding: 84.98%'

# %%
