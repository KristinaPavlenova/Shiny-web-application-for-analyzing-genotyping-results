import time
import requests
import logging
import glob
import os, os.path
from Bio import SeqIO
import re
from modules.seqs_process import extract_hits
from config import assembly, gene, idxs, log_file

os.makedirs(os.path.join(
    "main_pipeline_results", "top1_check"), exist_ok=True)
table_path = os.path.join("main_pipeline_results",
    "top1_check", f"blast_summary_{assembly}_{gene}.txt"
)  # указать файл для результатов проверки
dir_for_mafft_results = os.path.join("main_pipeline_results",
    f"mafft_results_{assembly}_{gene}"  # директория для результатов mafft-выравнивания
)
fastanames = glob.glob("*_trimmed_best.fasta", root_dir=dir_for_mafft_results)
fastanames = [fastanames[i] for i in idxs]
logging.basicConfig(
    filename=log_file,
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
)

if not os.path.exists(table_path):
    with open(table_path, "w") as table_file:
        table_file.write(
            f"{'Name':<40}{'Seq_ID':<30}{'Type':<15}"
            f"{'Subject_ID':<20}{'Organism':<40}"
            f"{'E-value':<10}{'Identity':<10}{'Align_len':<10}\n"
        )
        table_file.write("-" * 185 + "\n")

len_lower_limit = 1450 if gene == "COI" else 2200
len_upper_limit = 1700 if gene == "COI" else 2300

for fasta in fastanames:
    file = os.path.join(dir_for_mafft_results, fasta)
    name = os.path.basename(fasta).split("_trimmed_best.fasta")[0]

    with open(file, "r") as f:
        record = next(SeqIO.parse(f, "fasta"))
        seq_id = record.id
        clean_seq = str(record.seq).replace("-", "")
        if len(clean_seq) < len_lower_limit or len(clean_seq) > len_upper_limit:
            logging.warning(f"[{name}] Пропущена: длина {len(clean_seq)} нт.")
            continue

        sequence = f">query\n{str(record.seq).strip('-')}"
        logging.info(f"{name}: отправка BLAST запроса...")

        params = {
            "CMD": "Put",
            "PROGRAM": "blastn",
            "DATABASE": "core_nt",
            "QUERY": sequence,
            "MEGABLAST": "on",
        }
        try:
            response = requests.post(
                "https://blast.ncbi.nlm.nih.gov/blast/Blast.cgi", data=params
            )
            response.raise_for_status()
        except Exception as e:
            logging.error(f"{name}: ошибка при отправке запроса: {e}")
            continue

        rid_search = re.search(r"RID = ([\w-]+)", response.text)
        rid = rid_search.group(1) if rid_search else None
        if not rid:
            logging.error(f"{name}: RID не найден.")
            continue

        while True:
            try:
                check = requests.get(
                    "https://blast.ncbi.nlm.nih.gov/Blast.cgi",
                    params={"CMD": "Get", "RID": rid, "FORMAT_OBJECT": "SearchInfo"},
                )
                check.raise_for_status()
                if "Status=READY" in check.text:
                    break
                elif "Status=FAILED" in check.text or "Status=UNKNOWN" in check.text:
                    logging.error(f"{name}: ошибка статуса BLAST.")
                    break
                time.sleep(5)
            except Exception as e:
                logging.error(f"{name}: ошибка при проверке статуса: {e}")
                break

        try:
            result = requests.get(
                "https://blast.ncbi.nlm.nih.gov/Blast.cgi",
                params={"CMD": "Get", "RID": rid, "FORMAT_TYPE": "Text"},
            )
            result.raise_for_status()
        except Exception as e:
            logging.error(f"{name}: ошибка при получении результата: {e}")
            continue

        hits = extract_hits(result.text)
        if not hits:
            logging.warning(f"{name}: хиты не найдены.")
            continue

        best_by_identity = max(hits, key=lambda x: x["identity"])
        best_by_evalue = min(hits, key=lambda x: x["evalue"])

        with open(table_path, "a") as table_file:
            for hit_type, hit in [
                ("Best_Identity", best_by_identity),
                ("Best_Evalue", best_by_evalue),
            ]:
                table_file.write(
                    f"{name:<40}{seq_id:<30}\t"
                    f"{hit_type:<15}\t"
                    f"{hit['subject_id']:<20}\t"
                    f"{hit['organism'][:40]:<40}\t"
                    f"{hit['evalue']:<10.1e}\t"
                    f"{str(hit['identity']) + '%':<10}\t"
                    f"{hit['align_len']:<10}\n"
                )
print(f"Check2 'The best hits by e-value and identity of trimmed top1 sequences BLAST to NCBI database' done")