#goal is too return all sequences in fasta format that are found to be in a prodigal protein

#input 1 is a folder of prodigal files
#has both txt and fasta files, the txt file's title look like
#        ({seq_id}.fa_gene.coords.txt)

#then within that text file it looks like
#.       Line 1: Defenition ....
#.       Line 2: Features ....
#.           CDS   {Start}..{End} 
#.                  or
#.           CDS   complement({Start}..{End}) 
#plan to take every line after line 2 that .startsWith('CDS') and .strip()[2:]
#then .split('..')

#input 2 is a file of fasta files of all the small genomic proteing looking like 
#       (>{seq_id}:{start}:{end}:{direction})
#.      (Sequence)

#one function takes all the txt files and turns them into a list of tuples that look like 
#.       (id, start, end, direction)

#one function takes the fasta file
#   reads and turns the sequence into (id, start, end, direction, sequence)
#   only compare the start and end if the ids and direction match
#   compare that start.txt < start.fasta and end.txt > end.fasta
#   then write it as a new super_negative.fasta

import os
import sys

def create_prodigal_list(txt):
    prod_list = []
    direction = '+'
    with open(txt) as f:
        next(f)
        next(f)
        for line in f:
            line = str(line).strip().replace(" ", "")
            if line.startswith('CDS'):
                line = line[3:]
                if line.startswith('complement'):
                    direction = '-'
                    line = line[11:-1]
                line = line.split('..')
                start = line[0]
                if start.startswith('<'):
                    start = 0
                end = line[1]
                if end.startswith('>'):
                    end = int(end[1:]) + 1
                prod_list.append((start, end, direction))
    return prod_list

def create_fasta_list(fasta):
    fasta_list = []
    with open(fasta) as f:
        next(f)
        for line in f:
            line = str(line).strip().replace(" ", "")
            if line.startswith('>'):
                line=line[1:].split(':')
                start = line[0]
                end = line[1]
                direction = line[2]
            else:
                sequence = line
                fasta_list.append((start, end, direction, sequence))
    return fasta_list
            

def fasta_vs_txt(fasta, txt):
    super_negative_proteins = []
    prodigal_stats = create_prodigal_list(txt)
    fasta_stats = create_fasta_list(fasta)
    for fstart, fend, fdirection, fsequence in fasta_stats:
        for pstart, pend, pdirection in prodigal_stats:
            if fdirection == pdirection and int(fstart)>int(pstart) and int(fend)<int(pend):
                super_negative_proteins.append(fsequence)
    return super_negative_proteins

def list_to_fasta(sequences, outpath):
    with open(outpath, 'w') as f:
        for i, seq in enumerate(sequences, start=1):
            f.write(f">seq{i}\n{seq}\n")

directory_path = '/wistar/auslander/Timothy/microproteins_timka/contained_within_prodigal/all_genome_chunks/'
prodigal_path = '/wistar/auslander/microproteins/prodigal_translations/'


    
fasta_name = sys.argv[1]
fasta_name = fasta_name.removesuffix('.fa')

fasta_id = fasta_name.removeprefix('all_genomes_protein_with_loc_').removesuffix('.fasta')

txt_path = os.path.join(prodigal_path, f'{fasta_id}.fa_gene.coords.txt')
fa_path = os.path.join(directory_path, f'all_genomes_protein_with_loc_{fasta_name}.fasta')
    
print(f"Comparing {fa_path} vs {txt_path}")
    
negative_proteins = fasta_vs_txt(fa_path, txt_path)

output_file = f'/wistar/auslander/Timothy/microproteins_timka/contained_within_prodigal/negative_protein_chunks/negative_proteins_{fasta_id}.fa'
list_to_fasta(negative_proteins, output_file)
