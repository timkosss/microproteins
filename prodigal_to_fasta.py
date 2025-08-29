import os

directory_path = '/wistar/auslander/microproteins/prodigal_translations/'
output_file = '/wistar/auslander/Timothy/microproteins_timka/complete_prodigal.fa'
min_len = 15
max_len = 70

def filter_and_write_fasta(fasta_path, output):
    filtered_count = 0
    total_count = 0
    with open(fasta_path, 'r') as f:
        seq = ''
        header = ''
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if seq and min_len <= len(seq) <= max_len:
                    output.write(f"{header}\n{seq}\n")
                    filtered_count+=1
                header = line
                seq = ''
                total_count+=1
                if total_count % 1000 ==0:
                    print(f"{total_count} sequences read, {filtered_count} filtered so far...", flush=True)
            else:
                seq += line
        if seq and min_len <= len(seq) <= max_len:
            output.write(f"{header}\n{seq}\n")
            filtered_count+=1


with open(output_file, 'w') as fasta_out:
    print("Started Writing ")
    for file in os.listdir(directory_path):
        filename = os.fsdecode(file)
        if filename.endswith(".fa"):
            print(f"Filtering {filename}")
            fa_path = os.path.join(directory_path, filename)
            filter_and_write_fasta(fa_path, fasta_out)
            print(f"FILTERED {filename}")

print(f"All filtered sequences written to {output_file}")
