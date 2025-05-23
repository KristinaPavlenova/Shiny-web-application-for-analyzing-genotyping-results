import os
from Bio import SeqIO

INPUT_DIR = "./test_1_best_hit"
REFERENCE = "18S_reference.fasta"
OUTPUT_DIR = "./processed_results"
MAFFT = "mafft"
TRIMAL = "trimal"
os.makedirs(OUTPUT_DIR, exist_ok=True)
ref_length = len(next(SeqIO.parse(REFERENCE, "fasta")).seq)
print(f"Reference length: {ref_length} bp")

for filename in os.listdir(INPUT_DIR):
    if not filename.endswith("_18S_best.fasta"):
        continue
    species_name = filename.replace("_18S_best.fasta", "")
    input_file = os.path.join(INPUT_DIR, filename)
    output_prefix = os.path.join(OUTPUT_DIR, species_name)
    print(f"\nОбработка {species_name}...")
    combined_file = f"{output_prefix}_combined.fasta"  # temp union file
    with open(combined_file, "w") as f:
        # 1st ref then specie seq
        ref_seq = next(SeqIO.parse(REFERENCE, "fasta"))
        query_seq = next(SeqIO.parse(input_file, "fasta"))
        SeqIO.write([ref_seq, query_seq], f, "fasta")
    # mafft
    aligned_file = f"{output_prefix}_aligned.fasta"
    exit_code = os.system(f"{MAFFT} --auto {combined_file} > {aligned_file}")

    if exit_code != 0:
        print(f"ERROR: MAFFT failed for {filename}")
        os.remove(combined_file)
        continue
    # trimal
    trimmed_file = f"{output_prefix}_trimmed.fasta"
    exit_code = os.system(f"{TRIMAL} -in {aligned_file} -out {trimmed_file} -nogaps")

    if exit_code != 0:
        print(f"ERROR: trimal failed for {filename}")
        os.remove(combined_file)
        continue
    if not os.path.exists(trimmed_file):
        print(f"ERROR: Trimmed file not created for {filename}")
        os.remove(combined_file)
        continue
    # проверка длины
    try:
        trimmed_seqs = list(SeqIO.parse(trimmed_file, "fasta"))
        if len(trimmed_seqs) != 2:
            print(
                f"ERROR: Expected 2 sequences after trimming, got {len(trimmed_seqs)}"
            )
            os.remove(combined_file)
            continue
        trimmed_length = len(trimmed_seqs[1].seq)

        if abs(trimmed_length - ref_length) > 0.1 * ref_length:  # +-10%
            print(f"ERROR: {filename} - bad length after trimming")
            print(f"Trimmed length: {trimmed_length}, Reference: {ref_length}")
        else:
            with open(f"{output_prefix}_processed.fasta", "w") as out:
                SeqIO.write(trimmed_seqs[1], out, "fasta")
                print(f"Saved processed sequence (length: {trimmed_length})")
    except Exception as e:
        print(f"ERROR processing {filename}: {str(e)}")
    os.remove(combined_file)
print("\nProcessing complete!")
