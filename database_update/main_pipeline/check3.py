from Bio.Align import PairwiseAligner
import pandas as pd
import os, os.path
from modules.seqs_process import (
    collect_trimmed_seqs,
    load_sequences,
    align_and_score,
    align_wo_gaps,
    save_alignment_text,
)
from config import gene, dir_final_fasta

dir_with_mafft_results_rnaspades = os.path.join("main_pipeline_results",
    f"mafft_results_rnaspades_{gene}"
)  # директория с готовыми fasta-файлами rnaspades
dir_with_mafft_results_trinity = os.path.join("main_pipeline_results",
    f"mafft_results_trinity_{gene}"  # директория с готовыми fasta-файлами trinity
)
dir_results = os.path.join("main_pipeline_results",f"comparison_{gene}")  # директория для сохранения результатов выравнивания последовательностей из разных сборок

os.makedirs(dir_final_fasta, exist_ok=True)
fasta_rnaspades = collect_trimmed_seqs(
    dir_with_mafft_results_rnaspades, "rnaspades", gene, dir_final_fasta
)
fasta_trinity = collect_trimmed_seqs(
    dir_with_mafft_results_trinity, "trinity", gene, dir_final_fasta
)
os.makedirs(dir_results, exist_ok=True)

seqs_rnaspades = load_sequences(fasta_rnaspades)
seqs_trinity = load_sequences(fasta_trinity)
aligner = PairwiseAligner()
aligner.mode = "global"
aligner.match_score = 1.0
aligner.mismatch_score = -1.0
aligner.open_gap_score = -2.0
aligner.extend_gap_score = -0.5
common_ids = set(seqs_rnaspades.keys()) & set(seqs_trinity.keys())
ids_unique_trinity = set(seqs_trinity.keys()) - set(seqs_rnaspades.keys())
ids_unique_rnaspades = set(seqs_rnaspades.keys()) - set(seqs_trinity.keys())
results = []

for seq_id in common_ids:
    result = align_and_score(seqs_rnaspades[seq_id], seqs_trinity[seq_id], aligner)
    if len(seqs_rnaspades[seq_id]) != len(seqs_rnaspades[seq_id].strip("-")) or len(
        seqs_trinity[seq_id]
    ) != len(seqs_trinity[seq_id].strip("-")):
        result_wo_gaps = align_wo_gaps(
            seqs_rnaspades[seq_id], seqs_trinity[seq_id], aligner
        )

    alignment_file_path = os.path.join(dir_results, f"{seq_id}_alignment.txt")
    save_alignment_text(result["Alignment_Object"], alignment_file_path)

    results.append(
        {
            "Sequence_ID": seq_id,
            "Len_Seq_rnaspades": len(seqs_rnaspades[seq_id]),
            "Len_Seq_trinity": len(seqs_trinity[seq_id]),
            "Identity": round(result["Identity_Percentage"], 3),
            "Alignment_Len": result["Alignment_Length"],
            "Len_Seq_wo_gaps_rnaspades": len(seqs_rnaspades[seq_id].strip("-")),
            "Len_Seq_wo_gaps_trinity": len(seqs_trinity[seq_id].strip("-")),
            "Identity_wo_gaps": (
                round(result_wo_gaps["Identity_Percentage"], 2)
                if len(seqs_rnaspades[seq_id]) != len(seqs_rnaspades[seq_id].strip("-"))
                or len(seqs_trinity[seq_id]) != len(seqs_trinity[seq_id].strip("-"))
                else "n/a"
            ),
        }
    )

for seq_id in ids_unique_rnaspades:
    results.append(
        {
            "Sequence_ID": seq_id,
            "Len_Seq_rnaspades": len(seqs_rnaspades[seq_id]),
            "Len_Seq_trinity": "n/a",
            "Identity": "n/a",
            "Alignment_Len": "n/a",
            "Len_Seq_wo_gaps_rnaspades": len(seqs_rnaspades[seq_id].strip("-")),
            "Len_Seq_wo_gaps_trinity": "n/a",
            "Identity_wo_gaps": "n/a",
        }
    )

for seq_id in ids_unique_trinity:
    results.append(
        {
            "Sequence_ID": seq_id,
            "Len_Seq_rnaspades": "n/a",
            "Len_Seq_trinity": len(seqs_trinity[seq_id]),
            "Identity": "n/a",
            "Alignment_Len": "n/a",
            "Len_Seq_wo_gaps_rnaspades": "n/a",
            "Len_Seq_wo_gaps_trinity": len(seqs_trinity[seq_id].strip("-")),
            "Identity_wo_gaps": "n/a",
        }
    )

df = pd.DataFrame(results)
df = df.sort_values("Sequence_ID", ignore_index=True)
df.to_csv(os.path.join("main_pipeline_results", f"comparison_summary_{gene}.txt"), sep="\t", index=False)
print(f"Check3 done")
