import subprocess
import glob
import shutil
import os, os.path
import logging
from config import input_dir, log_file, reference, assembly, gene, idxs

os.makedirs("main_pipeline_results", exist_ok=True)
fastanames = glob.glob("*.fasta", root_dir=input_dir)
logging.basicConfig(filename=log_file, level=logging.ERROR)
dir_for_blast_results = os.path.join(
    "main_pipeline_results",
    f"blast_results_{assembly}_{gene}"  # директория для результатов blast-выравнивания
)
os.makedirs(dir_for_blast_results, exist_ok=True)
fastanames = [fastanames[i] for i in idxs]

for fasta in fastanames:
    file = os.path.join(input_dir, fasta)
    name = ".".join((fasta).split(".")[:-1])

    # Бластим референс на результаты сборки
    try:
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
        print(f"Step1 'BLAST alignment of reference to BLAST database with transcriptome results' for {name} done")
    except Exception as e:
        logging.error(f"Ошибка при blast-выравнивании для {name}: {e}")
        continue

