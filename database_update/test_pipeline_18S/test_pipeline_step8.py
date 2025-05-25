import os
from pathlib import Path

input_dir = Path("./final_results")
output_file = Path("fin_file.fasta")

if not input_dir.exists():
    print(f"Error: dir {input_dir} does not exist!")
    exit(1)

with open(output_file, 'w') as out_f:
    for fasta_file in input_dir.glob("*.fasta"):
        with open(fasta_file, 'r') as in_f:
            lines = in_f.readlines()
            header = lines[0].strip()
            
            # join lines with seqs to one line
            sequence = ''.join(line.strip() for line in lines[1:])
            sequence = sequence.upper()
            
            out_f.write(f"{header}\n{sequence}\n")

print(f"All seqs saved in {output_file}")
