import subprocess
import glob
import os, os.path
from Bio import SeqIO
import logging
from config import input_dir, log_file, reference, assembly, gene, top_n, idxs

fastanames = glob.glob("*.fasta", root_dir=input_dir)
logging.basicConfig(filename=log_file, level=logging.ERROR)
dir_for_blast_results = os.path.join(
    "main_pipeline_results",
    f"blast_results_{assembly}_{gene}"  # директория для результатов blast-выравнивания
)
fastanames = [fastanames[i] for i in idxs]

for fasta in fastanames:
    file = os.path.join(input_dir, fasta)
    name = ".".join((fasta).split(".")[:-1])

    # Топ хитов бласта: 1) сначала фильтруем по длине (+ покрытию для rnaspades); 2) далее сортируем по e-value (+ покрытию для rnaspades)
    try:
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
    except Exception as e:
        logging.error(f"Ошибка при получении топ хитов blast-выравнивания для {name}: {e}")
        continue
