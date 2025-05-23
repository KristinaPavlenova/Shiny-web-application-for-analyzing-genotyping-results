import subprocess
import glob
import shutil
import os, os.path
from Bio import SeqIO, AlignIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import pandas as pd
import logging
import re

input_dir = "input_fasta"  # указать директорию входных fasta-файлов
dir_for_blast_results = (
    "blast_results_COI"  # указать директорию для результатов blast-выравнивания
)
dir_for_mafft_results = (
    "mafft_results_COI"  # указать директорию для результатов mafft-выравнивания
)
os.makedirs("logs", exist_ok=True)
log_file = os.path.join("logs", "main_script_errors.log")  # указать файл для логов
reference = os.path.join(
    "references", "Eulimnogammarus_verrucosus_COI.fasta"
)  # указать fasta-файл с референсом
assembly = "rnaspades"  # указать вариант сборки (rnaspades или trinity)
gene = "COI"  # указать ген для поиска (COI или 18S)
top_n = 5  # указать количество лучших результатов после blast-выравнивания
fastanames = glob.glob("*.fasta", root_dir=input_dir)
idxs = list(
    range(len(fastanames))
)  # указать индексы (номер по порядку) входных fasta-файлов для пайплайна

logging.basicConfig(filename=log_file, level=logging.ERROR)
os.makedirs(dir_for_blast_results, exist_ok=True)
os.makedirs(dir_for_mafft_results, exist_ok=True)
fastanames = [fastanames[i] for i in idxs]

for fasta in fastanames:
    file = os.path.join(input_dir, fasta)
    name = ".".join((fasta).split(".")[:-1])

    # Бластим референс на результаты сборки
    subprocess.run(
        [
            "makeblastdb",
            "-in",
            file,
            "-title",
            file,
            "-dbtype",
            "nucl",
            "-out",
            os.path.join("databases", fasta),
        ]
    )
    subprocess.run(
        [
            "blastn",
            "-db",
            "databases/" + fasta,
            "-query",
            reference,
            "-evalue",
            "1e-3",
            "-word_size",
            "11",
            "-outfmt",
            "6 sseqid evalue slen",
            "-out",
            os.path.join(dir_for_blast_results, fasta),
        ]
    )
    shutil.rmtree("databases")

    # Топ хитов бласта: 1) сначала фильтруем по длине (+ покрытию для rnaspades); 2) далее сортируем по e-value (+ покрытию для rnaspades)
    blast_top_hits = os.path.join(dir_for_blast_results, name + f"_top{top_n}.txt")
    with open(blast_top_hits, "w") as out_file:
        len_filter = 1400 if gene == "COI" else 2000
        if assembly == "rnaspades":
            if gene == "COI":
                awk_extract = subprocess.Popen(
                    [
                        "awk",
                        r"""{match($1, /length_([0-9]+)/, len); match($1, /cov_([0-9.]+)/, cov);
                                                if (len[1] >= 1400 && cov[1] >= 10) {print $0 "\t" cov[1];}}""",
                        os.path.join(dir_for_blast_results, fasta),
                    ],
                    stdout=subprocess.PIPE,
                )
            else:
                awk_extract = subprocess.Popen(
                    [
                        "awk",
                        r"""{match($1, /length_([0-9]+)/, len); match($1, /cov_([0-9.]+)/, cov);
                                                if (len[1] >= 2000 && cov[1] >= 10) {print $0 "\t" cov[1];}}""",
                        os.path.join(dir_for_blast_results, fasta),
                    ],
                    stdout=subprocess.PIPE,
                )
            sort_e = subprocess.Popen(
                ["sort", "-k2,2g", "-k4,4gr", "-s"],
                stdin=awk_extract.stdout,
                stdout=subprocess.PIPE,
            )
        else:
            awk_extract = subprocess.Popen(
                [
                    "awk",
                    f"$3 >= {len_filter}",
                    os.path.join(dir_for_blast_results, fasta),
                ],
                stdout=subprocess.PIPE,
            )
            sort_e = subprocess.Popen(
                ["sort", "-k2,2g", "-s"],
                stdin=awk_extract.stdout,
                stdout=subprocess.PIPE,
            )
        cut_id = subprocess.Popen(
            ["cut", "-f", "1"], stdin=sort_e.stdout, stdout=subprocess.PIPE
        )
        subprocess.run(["head", "-n", f"{top_n}"], stdin=cut_id.stdout, stdout=out_file)
    with open(blast_top_hits, "r") as in_file:
        required_ids = [line.strip() for line in in_file.readlines()]
    # Создаем fasta референс + Топ хитов
    top_ref = os.path.join(dir_for_blast_results, name + f"_top{top_n}" + "_ref.fasta")
    with open(top_ref, "w") as out_file:
        for record in SeqIO.parse(reference, "fasta"):
            SeqIO.write(record, out_file, "fasta")
        for record in SeqIO.parse(file, "fasta"):
            if record.id in required_ids:
                SeqIO.write(record, out_file, "fasta")

    # Проводим множественное выравнивание
    try:
        mafft_result = os.path.join(dir_for_mafft_results, name + "_mafft.fasta")
        with open(mafft_result, "w") as out_file:
            subprocess.run(
                [
                    "mafft",
                    "--maxiterate",
                    "1000",
                    "--globalpair",
                    "--adjustdirection",
                    top_ref,
                ],
                stdout=out_file,
            )
    except subprocess.CalledProcessError as e:
        logging.error(f"Ошибка при запуске MAFFT для {name}: {e}")
        continue

    # Выбираем лучшую последовательность: 1) по идентичности к референсу; 2) по покрытию (для rnaspades)
    try:
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
        logging.error(f"Ошибка при записи результатов MAFFT для {name}: {e}")
        continue

    # Обрезаем лучшую последовательность по границам референса
    try:
        trimmed_best = os.path.join(dir_for_mafft_results, name + "_trimmed_best.fasta")
        with open(trimmed_best, "w") as out_file:
            alignment = AlignIO.read(best_to_ref, "fasta")
            reference_sequence = alignment[0].seq
            start = next(
                (i for i, char in enumerate(reference_sequence) if char != "-"), -1
            )
            end = (
                len(reference_sequence)
                - next(
                    (
                        i
                        for i, char in enumerate(reference_sequence[::-1])
                        if char != "-"
                    ),
                    -1,
                )
                - 1
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
    except Exception as e:
        logging.error(f"Ошибка при обрезке последовательности {name}: {e}")
        continue
