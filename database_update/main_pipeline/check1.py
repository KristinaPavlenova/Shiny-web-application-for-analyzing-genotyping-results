import glob
import os, os.path
from Bio import SeqIO
from config import gene, idxs, assembly

os.makedirs(os.path.join(
    "main_pipeline_results", "top1_check"), exist_ok=True)
table_path = os.path.join("main_pipeline_results",
    "top1_check", f"seqs_len_exceptions_{assembly}_{gene}.txt"
)  # файл для результатов проверки
dir_for_mafft_results = os.path.join("main_pipeline_results",
    f"mafft_results_{assembly}_{gene}"  # директория для результатов mafft-выравнивания
)
fastanames = glob.glob("*_trimmed_best.fasta", root_dir=dir_for_mafft_results)
fastanames = [fastanames[i] for i in idxs]

if not os.path.exists(table_path):
    with open(table_path, "w") as table_file:
        table_file.write(
            f"{'Name':<40}\t"
            f"{'Seq_id':<50}\t"
            f"{'Seq_length':<10}\t"
            f"{'Cleaned_seq_length':<10}\t"
            f"{'Without_NNN'}\n"
        )
        table_file.write("-" * 120 + "\n")

len_lower_limit = 1531 if gene == "COI" else 2260
len_upper_limit = 1538 if gene == "COI" else 2274

for fasta in fastanames:
    file = os.path.join(dir_for_mafft_results, fasta)
    name = os.path.basename(fasta).split("_trimmed_best.fasta")[0]
    with open(file, "r") as f:
        record = next(SeqIO.parse(f, "fasta"))
        seq_id = record.id
        clean_seq = str(record.seq).replace("-", "")
        if (
            (len(clean_seq) >= len_lower_limit and len(clean_seq) <= len_upper_limit)
            and len(clean_seq) == len(record.seq)
            and "n" not in clean_seq
        ):
            continue
        with open(table_path, "a") as table_file:
            table_file.write(
                f"{name:<40}\t"
                f"{seq_id:<50}\t"
                f"{len(record.seq):<10}\t"
                f"{len(clean_seq):<10}\t"
                f"{'n' not in clean_seq}\n"
            )
print(f"Check1 done")
