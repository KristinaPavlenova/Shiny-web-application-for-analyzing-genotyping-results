import os
from Bio import SeqIO

INPUT_DIR = "./processed_results"
REFERENCE = "18S_reference.fasta"
OUTPUT_DIR = "./final_results"
BLASTDB_DIR = "./blast_reference_db"
MIN_IDENTITY = 90.0
os.makedirs(OUTPUT_DIR, exist_ok=True)
os.makedirs(BLASTDB_DIR, exist_ok=True)

# BLAST бд из референса
print("Создаём BLAST базу данных...")
base_name = os.path.splitext(os.path.basename(REFERENCE))[0]
db_command = (
    f"makeblastdb -in {REFERENCE} -blastdb_version 5 "
    f"-dbtype nucl -hash_index -out {BLASTDB_DIR}/{base_name}"
)
os.system(db_command)

print("\nПроверка последовательности...")
for filename in os.listdir(INPUT_DIR):
    if not filename.endswith("_processed.fasta"):
        continue
    input_file = os.path.join(INPUT_DIR, filename)
    output_file = os.path.join(OUTPUT_DIR, filename.replace("_processed", "_final"))
    # BLAST
    blast_output = "temp_blast_result.txt"
    os.system(
        f"blastn -query {input_file} -db {BLASTDB_DIR}/{base_name} "
        f"-out {blast_output} -outfmt '6 pident' -max_target_seqs 1"
    )

    try:
        with open(blast_output) as f:
            identity = float(f.readline().split()[0])
        print(f"{filename}: идентичность {identity:.1f}%", end=" ")

        if identity >= MIN_IDENTITY:
            with open(input_file) as src, open(output_file, "w") as dst:
                dst.write(src.read())
            print("Сохранено")
        else:
            print("Слишком низкая идентичность")
    except Exception as e:
        print(f"{filename}: ошибка при проверке ({str(e)})")
    if os.path.exists(blast_output):
        os.remove(blast_output)
print("\nПроверка закончена!")
