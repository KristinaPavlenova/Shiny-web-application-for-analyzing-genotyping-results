import os
from pathlib import Path

input_file = Path("./filtered_sequences.fasta")
output_file = Path("./final_file.fasta")

if not input_file.exists():
    print(f"Error: file {input_file} does not exist!")
    exit(1)

with open(output_file, 'w') as out_f:
        with open(input_file, 'r') as in_f:
            sequence = ""
            for line in in_f:
                if line[0] == ">":
                    out_f.write(f"{sequence}\n")
                    sequence = ""
                    out_f.write(line)
                else:
                    sequence += line.strip()
            out_f.write(sequence)
            
print(f"All seqs saved in {output_file}")
