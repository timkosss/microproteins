from Bio.Seq import Seq
import pandas as pd

fasta_path = '/wistar/auslander/microproteins/all_genomes.fa'
#fasta_path = '/wistar/auslander/Timothy/microproteins_timka/small_genome.fa'
df = pd.read_csv('/wistar/auslander/microproteins/TableS2c.csv')

def read_in_fasta(fasta_path):
    print("Reading Fasta")
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
    print("Done Reading Fasta")
    return sequences

def dna_to_protein(sequences):
    proteins=[]
    for i, sequence in enumerate(sequences):
        seq = Seq(sequence)
        revseq = seq.reverse_complement()
        for nuc_seq in [seq, revseq]:
            for frame in range(3):
                translated = str(nuc_seq[frame:].translate(to_stop=False))
                proteins.extend(get_short_prot(translated))
        print(f"Sequence {i} is done", flush=True)
    return proteins

def get_short_prot(seq):
    proteins = []
    for prot in seq.split('*'):
        if 'M' in prot:  
            prot = prot[prot.find('M'):]
            if 15 <= len(prot) <= 70:
                proteins.append(prot)
    return proteins

def pd_to_fasta(df):
    with open("/wistar/auslander/Timothy/microproteins_timka/output.fasta", "w") as f:
        for i, seq in enumerate(df['Sequence'], start=1):
            f.write(f">seq{i}\n{seq}\n")

def list_to_fasta(sequences):
    with open("/wistar/auslander/Timothy/microproteins_timka/output2.fasta", "w") as f:
        for i, seq in enumerate(sequences, start=1):
            f.write(f">seq{i}\n{seq}\n")
               
list_to_fasta(dna_to_protein(read_in_fasta(fasta_path)))

