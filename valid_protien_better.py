import re
import random

######################### 6-frame translation #########################

def handle_non_ATGC(sequence):
    """Replace non-ATGC with random base."""
    ret = re.sub('[^ATGC]', lambda x: random.choice(['A','C','G','T']), sequence)
    assert len(ret) == len(sequence)
    return ret

def translate_frameshifted(sequence, gcode):
    """Translate sequence using codon table."""
    translate =  ''.join([gcode.get(sequence[3*i:3*i+3], 'X') for i in range(len(sequence)//3)])
    return translate

def reverse_complement(sequence, bpairs):
    """Reverse complement of a nucleotide sequence."""
    reversed_sequence = (sequence[::-1])
    rc = ''.join([bpairs.get(reversed_sequence[i]) for i in range(len(sequence))])
    return rc

def six_frame_trans(seq, gcode, bpairs):
    """Return 6-frame translations."""
    x1 = translate_frameshifted(seq, gcode)
    x2 = translate_frameshifted(seq[1:], gcode)
    x3 = translate_frameshifted(seq[2:], gcode)
    rc = reverse_complement(seq, bpairs)
    x4 = translate_frameshifted(rc, gcode)
    x5 = translate_frameshifted(rc[1:], gcode)
    x6 = translate_frameshifted(rc[2:], gcode)
    x = [x1, x2, x3, x4, x5, x6]
    return x

def get_short_prot(seq):
    """Return short proteins 15-70aa starting with M."""
    proteins = []
    for prot in seq.split('*'):
        if 'M' in prot:
            prot = prot[prot.find('M'):]
            if 15 <= len(prot) <= 70:
                proteins.append(prot)
    return proteins

######################### Fasta reading #########################

def read_in_fasta(fasta_path):
    print("staring read in")
    sequences = []
    with open(fasta_path, 'r') as f:
        seq = ''
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if seq:
                    sequences.append(seq)
                    seq = ''
            else:
                seq += line
        if seq:
            sequences.append(seq)
    print("finished reading")
    return sequences

######################### Main processing #########################

def dna_to_protein(sequences, gcode, bpairs):
    proteins=[]
    for i, sequence in enumerate(sequences):
        seq = handle_non_ATGC(sequence.strip())
        frames = six_frame_trans(seq, gcode, bpairs)
        for translated in frames:
            proteins.extend(get_short_prot(translated))
        print(f"Sequence {i} done", flush=True)
    return proteins

def list_to_fasta(sequences, outpath):
    with open(outpath, 'w') as f:
        for i, seq in enumerate(sequences, start=1):
            f.write(f">seq{i}\n{seq}\n")

######################### Example usage #########################

# Standard genetic code table
gcode = {
    'TTT':'F','TTC':'F','TTA':'L','TTG':'L',
    'CTT':'L','CTC':'L','CTA':'L','CTG':'L',
    'ATT':'I','ATC':'I','ATA':'I','ATG':'M',
    'GTT':'V','GTC':'V','GTA':'V','GTG':'V',
    'TCT':'S','TCC':'S','TCA':'S','TCG':'S',
    'CCT':'P','CCC':'P','CCA':'P','CCG':'P',
    'ACT':'T','ACC':'T','ACA':'T','ACG':'T',
    'GCT':'A','GCC':'A','GCA':'A','GCG':'A',
    'TAT':'Y','TAC':'Y','TAA':'*','TAG':'*',
    'CAT':'H','CAC':'H','CAA':'Q','CAG':'Q',
    'AAT':'N','AAC':'N','AAA':'K','AAG':'K',
    'GAT':'D','GAC':'D','GAA':'E','GAG':'E',
    'TGT':'C','TGC':'C','TGA':'*','TGG':'W',
    'CGT':'R','CGC':'R','CGA':'R','CGG':'R',
    'AGT':'S','AGC':'S','AGA':'R','AGG':'R',
    'GGT':'G','GGC':'G','GGA':'G','GGG':'G'
}

bpairs = {'A':'T','T':'A','C':'G','G':'C'}

fasta_path = '/wistar/auslander/microproteins/all_genomes.fa'
outpath = '/wistar/auslander/Timothy/microproteins_timka/all_genomes_protein.fasta'

seqs = read_in_fasta(fasta_path)
proteins = dna_to_protein(seqs, gcode, bpairs)
list_to_fasta(proteins, outpath)
