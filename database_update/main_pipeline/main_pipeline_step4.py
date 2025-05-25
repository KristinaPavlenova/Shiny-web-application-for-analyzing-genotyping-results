import glob
import os, os.path
from Bio import SeqIO, AlignIO
import pandas as pd
import logging
import re
from config import input_dir, log_file, assembly, idxs, assembly, gene

fastanames = glob.glob("*.fasta", root_dir=input_dir)
logging.basicConfig(filename=log_file, level=logging.ERROR)
dir_for_mafft_results = os.path.join(
    "main_pipeline_results",
    f"mafft_results_{assembly}_{gene}"  # директория для результатов mafft-выравнивания
)
fastanames = [fastanames[i] for i in idxs]

for fasta in fastanames:
    file = os.path.join(input_dir, fasta)
    name = ".".join((fasta).split(".")[:-1])

    # Выбираем лучшую последовательность: 1) по идентичности к референсу; 2) по покрытию (для rnaspades)
    try:
        mafft_result = os.path.join(dir_for_mafft_results, name + "_mafft.fasta")
        best_to_ref = os.path.join(dir_for_mafft_results, name + "_best_to_ref.fasta")
        with open(best_to_ref, "w") as out_file:
            alignment = AlignIO.read(mafft_result, "fasta")
            reference_sequence = alignment[0]
            results = []
            for record in alignment:
                matches = sum(
                    1
                    for a, b in zip(reference_sequence.seq, record.seq)
                    if a == b and a != "-"
                )
                identity = matches / len(reference_sequence.seq) * 100

                if assembly == "rnaspades":
                    match = re.search(r"cov[_\-]?(\d+\.?\d*)", record.id)
                    coverage = float(match.group(1)) if match else 0.0
                    results.append((record.id, identity, coverage))
                else:
                    results.append((record.id, identity))
            if assembly == "rnaspades":
                df = pd.DataFrame(
                    results, columns=["Sequence ID", "Avg Identity (%)", "Coverage"]
                )
                df_sorted = df.sort_values(
                    by=["Avg Identity (%)", "Coverage"],
                    ascending=False,
                    ignore_index=True,
                )
            else:
                df = pd.DataFrame(results, columns=["Sequence ID", "Avg Identity (%)"])
                df_sorted = df.sort_values(
                    by=["Avg Identity (%)"], ascending=False, ignore_index=True
                )
            for record in SeqIO.parse(mafft_result, "fasta"):
                if record.id in list(df_sorted["Sequence ID"][:2]):
                    SeqIO.write(record, out_file, "fasta")
    except Exception as e:
        logging.error(f"Ошибка при выборе топ1 результатов MAFFT для {name}: {e}")
        continue
