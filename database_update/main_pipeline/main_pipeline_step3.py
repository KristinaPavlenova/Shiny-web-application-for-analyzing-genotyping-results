import subprocess
import glob
import os, os.path
import logging
from config import input_dir, log_file, top_n, idxs, assembly, gene

fastanames = glob.glob("*.fasta", root_dir=input_dir)
logging.basicConfig(filename=log_file, level=logging.ERROR)
dir_for_blast_results = os.path.join(
    "main_pipeline_results",
    f"blast_results_{assembly}_{gene}"  # директория для результатов blast-выравнивания
)
dir_for_mafft_results = os.path.join(
    "main_pipeline_results",
    f"mafft_results_{assembly}_{gene}"  # директория для результатов mafft-выравнивания
)
os.makedirs(dir_for_mafft_results, exist_ok=True)
fastanames = [fastanames[i] for i in idxs]

for fasta in fastanames:
    file = os.path.join(input_dir, fasta)
    name = ".".join((fasta).split(".")[:-1])

    # Проводим множественное выравнивание
    try:
        top_ref = os.path.join(dir_for_blast_results, name + f"_top{top_n}" + "_ref.fasta")
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
