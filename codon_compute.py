#!/usr/bin/env python3
#%%
import os, gzip, itertools

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

url1="ftp://ftp.ensemblgenomes.org/pub/bacteria/release-45/fasta/bacteria_0_collection/salmonella_enterica_subsp_enterica_serovar_typhimurium_str_lt2/cds/Salmonella_enterica_subsp_enterica_serovar_typhimurium_str_lt2.ASM694v2.cds.all.fa.gz"
url2="ftp://ftp.ensemblgenomes.org/pub/bacteria/release-45/fasta/bacteria_0_collection/mycobacterium_tuberculosis_h37rv/cds/Mycobacterium_tuberculosis_h37rv.ASM19595v2.cds.all.fa.gz"
file1="Salmonella_enterica_subsp_enterica_serovar_typhimurium_str_lt2.ASM694v2.cds.all.fa.gz"
file2="Mycobacterium_tuberculosis_h37rv.ASM19595v2.cds.all.fa.gz"

if not os.path.exists(file1):
    os.system("curl -O %s"%(url1))

if not os.path.exists(file2):
    os.system("curl -O %s"%(url2))

from collections import defaultdict 

def get_answers_to_questions(files, species_labels, nucleotides):
    codons_seen_in_each_file = defaultdict(dict) # { AUG: { SP1: int, SP2: int, ...} }
    species_to_codons = defaultdict(dict) # { SP1: [ AUG, ATC, ...], SP2: [ AUG, GGT, ...] }
    total_number_of_genes_in_each_file = defaultdict(dict) # { SP1: int, SP2: int, ...}
    total_length_of_genes_in_each_file = defaultdict(dict) # { SP1: int, SP2: int, ...}
    g_plus_c_count_in_each_file = defaultdict(dict) # { SP1: int, SP2: int, ...}
    total_codons_in_each_file = defaultdict(dict) # { SP1: int, SP2: int, ...}
    species_analyzed = 0

    for file in files:
        with gzip.open(file,"rt") as fh:
            seqs = aspairs(fh)

            for seq in seqs:
                seqname  = seq[0]
                seqstring = seq[1]
                spn = species_labels[species_analyzed]

                if spn not in total_number_of_genes_in_each_file:
                    total_number_of_genes_in_each_file[
                        spn
                    ] = 1
                else:
                    total_number_of_genes_in_each_file[
                        spn
                    ] += 1                  

                total_length_of_sequence = len(seqstring)

                if spn not in total_length_of_genes_in_each_file:
                    total_length_of_genes_in_each_file[
                        spn
                    ] = total_length_of_sequence
                else:
                    total_length_of_genes_in_each_file[
                        spn
                    ] += total_length_of_sequence
                
                for i in range(0, total_length_of_sequence, 3):
                    codon = seqstring[i:i+3]

                    # piggy back off each codon to check to see if G or C are in there, increase counts
                    for nucleotide in nucleotides:
                        if nucleotide in codon:
                            if spn not in g_plus_c_count_in_each_file:
                                g_plus_c_count_in_each_file[
                                    spn
                                ] = 1
                            else:
                                g_plus_c_count_in_each_file[
                                    spn
                                ] += 1

                    if len(codon) == 3:
                        if spn not in species_to_codons:
                            species_to_codons[spn] = [codon]
                        else:
                            species_to_codons[spn].append(codon)

                        if codon not in codons_seen_in_each_file:
                            # add this codon to our hash table and initialize with +1 for the current species we're on and 0 for the others
                            codons_seen_in_each_file[
                                codon
                            ] = defaultdict(dict)
                            for species in species_labels:
                                codons_seen_in_each_file[
                                    codon
                                ][
                                    species 
                                ] = 1 if species == spn else 0
                        else:
                            codons_seen_in_each_file[
                                codon
                            ][
                                spn 
                            ] += 1


        species_analyzed += 1

    total_number_of_genes_in_each_file_display_txt = []
    for species in total_number_of_genes_in_each_file:
        total_number_of_genes_in_each_file_display_txt.append(
            '{}: {:,}'.format(species, total_number_of_genes_in_each_file[species])
        )

    total_length_of_genes_in_each_file_display_txt = []
    g_plus_c_frequency_in_each_file = {}
    for species in total_length_of_genes_in_each_file:
        total_length_of_genes_in_each_file_display_txt.append(
            '{}: {:,}'.format(species, total_length_of_genes_in_each_file[species])
        )
        g_plus_c_frequency_in_each_file[species] = g_plus_c_count_in_each_file[species] / total_length_of_genes_in_each_file[species]
    
    g_plus_c_frequency_in_each_file_display_txt = []
    for species in g_plus_c_frequency_in_each_file:
        g_plus_c_frequency_in_each_file_display_txt.append(
            '{}: {:.2%}'.format(species, g_plus_c_frequency_in_each_file[species])
        )
    
    codons_seen_in_each_file_display_txt = []
    for species in species_labels:
        codons_seen_in_each_file_display_txt.append(
            '{}: {:,}'.format(species, len(species_to_codons[species]))
        )
    
    codons_display_table  = [
        '{:^20}{:^30}{:^18}'.format('Codon', 'Frequency in Sp1', 'Frequency in Sp2')
    ]
    for codon in codons_seen_in_each_file:
        codons_display_row = ['{:^20}'.format(codon)]
        for species in species_labels:
            codons_display_row.append('{:^20.2%}'.format(codons_seen_in_each_file[codon][species]/len(species_to_codons[species])))
        codons_display_table.append('\t'.join(codons_display_row))

    return '\n\n'.join([    
        '1. The total number of genes in each species: ' + '; '.join(total_number_of_genes_in_each_file_display_txt),
        '2. Total length of these gene sequences for each file: ' + '; '.join(total_length_of_genes_in_each_file_display_txt),
        '3. The G+C percentage for the whole dataset (eg the frequency of G + the frequency of C): ' + '; '.join(g_plus_c_frequency_in_each_file_display_txt),
        '4. Total number codons in each genome.: ' + '; '.join(codons_seen_in_each_file_display_txt),
        '5. Print out table with three columns: Codon, Frequency in Sp1, Frequency in Sp2:',
        '\n'.join(codons_display_table)
    ])

answers = get_answers_to_questions([file1, file2], ['SP1', 'SP2'], ['G', 'C'])

print(answers)

# %%
