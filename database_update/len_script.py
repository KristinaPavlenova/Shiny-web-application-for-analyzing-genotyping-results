import glob
import os
from Bio import SeqIO




fasta_files = glob.glob("mafft_results_18S/*_trimmed_best.fasta")
os.makedirs('top1_check', exist_ok=True)

for fasta_file in fasta_files:
    name = fasta_file.split('_trimmed_best.fasta')[0].split('/')[1]
    with open(fasta_file, "r") as f:
        record = next(SeqIO.parse(f, "fasta"))
        seq_id = record.id
        
        clean_seq = str(record.seq).replace('-', '')
        
        if (len(clean_seq) >= 2260 and len(clean_seq) <= 2274) and len(clean_seq) == len(record.seq) and 'n' not in clean_seq:
             continue

        table_path = 'top1_check/seqs_len_18S_after.txt'
        with open(table_path, 'a') as table_file:
                table_file.write(
                f"{name:<40}\t"
                f"{seq_id:<50}\t"
                f"{len(record.seq):<10}\t"
                f"{len(clean_seq):<10}\t"
                f"{'n' not in clean_seq}\n")