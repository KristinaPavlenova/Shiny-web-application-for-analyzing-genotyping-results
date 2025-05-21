import time
import requests
import logging
import glob
import os
from Bio import SeqIO
import re

logging.basicConfig(filename='blast_request.log', level=logging.INFO,
                    format='%(asctime)s - %(levelname)s - %(message)s')

os.makedirs('top1_check', exist_ok=True)
table_path = 'top1_check/blast_summary.txt'

if not os.path.exists(table_path):
    with open(table_path, 'w') as table_file:
        table_file.write(f"{'Name':<40}{'Seq_ID':<30}{'Type':<15}"
                         f"{'Subject_ID':<20}{'Organism':<40}"
                         f"{'E-value':<10}{'Identity':<10}{'Align_len':<10}\n")
        table_file.write('-' * 185 + '\n')


def extract_hits(blast_text):
    hits = []
    entries = re.split(r"^>", blast_text, flags=re.MULTILINE)[1:]
    for entry in entries:
        lines = entry.strip().splitlines()
        if not lines:
            continue
        header = lines[0]
        match_header = re.match(r"([\w|.]+) (.+)", header)
        subject_id = match_header.group(1) if match_header else 'N/A'
        organism = match_header.group(2).strip() if match_header else 'N/A'

        identities_match = re.search(r'Identities\s*=\s*(\d+)\s*/\s*(\d+)', entry)
        evalue_match = re.search(r'Expect *= *([\deE.-]+)', entry)

        if identities_match and evalue_match:
            matches = int(identities_match.group(1))
            align_len = int(identities_match.group(2))
            identity = round((matches / align_len) * 100, 2)
            evalue = float(evalue_match.group(1).replace('e', 'E'))
            hits.append({
                "subject_id": subject_id,
                "organism": organism,
                "identity": identity,
                "align_len": align_len,
                "evalue": evalue
            })
    return hits


fasta_files = glob.glob("mafft_results_new1/*_trimmed_best.fasta")

for fasta_file in fasta_files:
    name = os.path.basename(fasta_file).split('_trimmed_best.fasta')[0]

    with open(fasta_file, "r") as f:
        record = next(SeqIO.parse(f, "fasta"))
        seq_id = record.id
        clean_seq = str(record.seq).replace('-', '')
        if len(clean_seq) < 1450 or len(clean_seq) > 1700:
            logging.warning(f"[{name}] Пропущена: длина {len(clean_seq)} нт.")
            continue

        sequence = f">query\n{str(record.seq).strip('-')}"
        logging.info(f"{name}: отправка BLAST запроса...")

        params = {'CMD': 'Put', 'PROGRAM': 'blastn', 'DATABASE': 'core_nt', 'QUERY': sequence, 'MEGABLAST': 'on'}
        try:
            response = requests.post('https://blast.ncbi.nlm.nih.gov/blast/Blast.cgi', data=params)
            response.raise_for_status()
        except Exception as e:
            logging.error(f"{name}: ошибка при отправке запроса: {e}")
            continue

        rid_search = re.search(r'RID = ([\w-]+)', response.text)
        rid = rid_search.group(1) if rid_search else None
        if not rid:
            logging.error(f"{name}: RID не найден.")
            continue

        while True:
            try:
                check = requests.get('https://blast.ncbi.nlm.nih.gov/Blast.cgi',
                                     params={'CMD': 'Get', 'RID': rid, 'FORMAT_OBJECT': 'SearchInfo'})
                check.raise_for_status()
                if 'Status=READY' in check.text:
                    break
                elif 'Status=FAILED' in check.text or 'Status=UNKNOWN' in check.text:
                    logging.error(f"{name}: ошибка статуса BLAST.")
                    break
                time.sleep(5)
            except Exception as e:
                logging.error(f"{name}: ошибка при проверке статуса: {e}")
                break

        try:
            result = requests.get('https://blast.ncbi.nlm.nih.gov/Blast.cgi',
                                  params={'CMD': 'Get', 'RID': rid, 'FORMAT_TYPE': 'Text'})
            result.raise_for_status()
        except Exception as e:
            logging.error(f"{name}: ошибка при получении результата: {e}")
            continue

        hits = extract_hits(result.text)
        if not hits:
            logging.warning(f"{name}: хиты не найдены.")
            continue

        best_by_identity = max(hits, key=lambda x: x['identity'])
        best_by_evalue = min(hits, key=lambda x: x['evalue'])

        with open(table_path, 'a') as table_file:
            for hit_type, hit in [('Best_Identity', best_by_identity), ('Best_Evalue', best_by_evalue)]:
                table_file.write(
                    f"{name:<40}{seq_id:<30}{hit_type:<15}"
                    f"{hit['subject_id']:<20}{hit['organism'][:40]:<40}"
                    f"{hit['evalue']:<10.1e}{str(hit['identity']) + '%':<10}{hit['align_len']:<10}\n"
                )
