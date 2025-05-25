import glob
import os, os.path
from Bio import SeqIO, AlignIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import logging
from config import input_dir, log_file, idxs, assembly, gene

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

    # Обрезаем лучшую последовательность по границам референса
    try:
        best_to_ref = os.path.join(dir_for_mafft_results, name + "_best_to_ref.fasta")
        trimmed_best = os.path.join(dir_for_mafft_results, name + "_trimmed_best.fasta")
        with open(trimmed_best, "w") as out_file:
            alignment = AlignIO.read(best_to_ref, "fasta")
            reference_sequence = alignment[0].seq
            start = next(
                (i for i, char in enumerate(reference_sequence) if char != "-"), -1
            )
            end = (
                len(reference_sequence) - next(
                    (i for i, char in enumerate(reference_sequence[::-1]) if char != "-"), - 1,
                ) - 1
            )
            trimmed_seq = alignment[-1].seq[start : end + 1]
            trimmed_ref = reference_sequence[start : end + 1]
            cleaned_seq = "".join(
                base
                for base, ref_base in zip(trimmed_seq, trimmed_ref)
                if not (base == "-" and ref_base == "-")
            )
            trimmed_record = SeqRecord(
                Seq(cleaned_seq), id=alignment[-1].id, description=""
            )
            SeqIO.write(trimmed_record, out_file, "fasta")
        print(f"Step5 for {name} done")
    except Exception as e:
        logging.error(f"Ошибка при обрезке последовательности {name}: {e}")
        continue
