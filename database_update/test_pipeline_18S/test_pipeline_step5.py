import os
from Bio import SeqIO

folder_path = "./5_hits_seqs"
folder_res = "./test_1_best_hit"
os.makedirs(folder_res, exist_ok=True)


def get_best_hit(file_path):
    with open(file_path, "r") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            # лучший хит —первый контиг
            return f">{record.id}\n{record.seq}"
    return None


def process_all_files(folder_path):
    files_in_folder = os.listdir(folder_path)
    fasta_files = [file for file in files_in_folder if file.endswith(".fasta")]
    for fasta_file in fasta_files:
        file_path = os.path.join(folder_path, fasta_file)
        best_hit = get_best_hit(file_path)  # первый контиг
        if best_hit:
            new_file_name = f"{os.path.splitext(fasta_file)[0][:-5]}_best.fasta"
    new_file_path = os.path.join(folder_res, new_file_name)
    with open(new_file_path, "w") as new_handle:
        new_handle.write(best_hit)


process_all_files(folder_path)
