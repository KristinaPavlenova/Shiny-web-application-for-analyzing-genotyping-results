import os
from Bio import SeqIO

best_hits_dir = "./best_hits"
fasta_dir = "/media/tertiary/transcriptome_assemblies/rnaspades_reassemblies"
output_dir = "./5_hits_seqs_TEST"
os.makedirs(output_dir, exist_ok=True)

hit_files = [
    f for f in os.listdir(best_hits_dir) if f.endswith("_vs_18S_reference_top5.txt")
]

for hit_file in hit_files:
    species_name = hit_file[:-26]
    fasta_filename = species_name + ".fasta"
    fasta_path = os.path.join(fasta_dir, fasta_filename)
    if not os.path.exists(fasta_path):
        print(f"Файл {fasta_path} не найден, пропускаем...")
        continue

    best_hits_path = os.path.join(best_hits_dir, hit_file)
    with open(best_hits_path, "r") as f:
        best_hits = [
            line.strip().split("\t")[1] for line in f.readlines()
        ]  # contig names
    output_fasta_path = os.path.join(output_dir, f"{species_name}_18S_top5.fasta")
    # извлечение последовательности
    with open(output_fasta_path, "w") as out_f:
        for record in SeqIO.parse(fasta_path, "fasta"):
            if record.id in best_hits:
                out_f.write(
                    f">Naumenko_18S_{species_name}_{record.id}\n{str(record.seq)}\n"
                )
print(f"Сохранено в файле {output_fasta_path}")
