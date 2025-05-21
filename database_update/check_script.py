import time
import requests
import logging
import glob
import os
from Bio import SeqIO
import re
from bs4 import BeautifulSoup

logging.basicConfig(filename='blast_request_new1.log', level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

os.makedirs('top1_check', exist_ok=True)
table_path = 'top1_check/blast_summary_new1.txt'
if not os.path.exists(table_path):
    with open(table_path, 'w') as table_file:
        table_file.write(f"{'Name':<45}{'Sequence_ID':<50}{'Subject_ID':<15}{'Organism':<50}{'E-value':<10}{'Identity':<10}{'Align_len':<10}\n")
        table_file.write('-' * 190 + '\n')

fasta_files = glob.glob("mafft_results_new1/*_trimmed_best.fasta")


for fasta_file in fasta_files[26:]:
    name = fasta_file.split('_trimmed_best.fasta')[0].split('/')[1]
    with open(fasta_file, "r") as f:
        record = next(SeqIO.parse(f, "fasta"))
        seq_id = record.id
        
        clean_seq = str(record.seq).replace('-', '')
        if len(clean_seq) < 1450:
            logging.error(f"[{name}] Ошибка: длина последовательности < 1450 ({len(clean_seq)} нуклеотидов). Запрос не отправлен.")
            continue
        if len(clean_seq) > 1700:
            logging.error(f"[{name}] Ошибка: длина последовательности > 1700 ({len(clean_seq)} нуклеотидов). Запрос не отправлен.")
            continue
        
        sequence = f">query\n{str(record.seq).strip('-')}"
      
        logging.info(f"Отправка запроса на BLAST для последовательности {name} из файла {fasta_file}...")
        params = {'CMD': 'Put', 'PROGRAM': 'blastn', 'DATABASE': 'core_nt', 'QUERY': sequence}

        try:
            response = requests.post('https://blast.ncbi.nlm.nih.gov/blast/Blast.cgi', data=params)
            response.raise_for_status()
        except Exception as e:
            logging.error(f"{name}, ошибка при отправке запроса: {e}")
            raise
        # logging.info(f"{name}, Ответ от сервера: {response.text}")
        
        soup = BeautifulSoup(response.content, 'html.parser')
        rid_input = soup.find('input', {'name': 'RID'})
        rid = rid_input.get('value')

        if not rid: 
        # or not rtoe:
            logging.error("RID или RTOE не найдены в ответе BLAST")
            raise ValueError(f"{name}, RID или RTOE не найдены")

        logging.info(f"{name}, RID получен: {rid}")
        # time.sleep(rtoe + 5)

        while True:
            check_params = {'CMD': 'Get', 'RID': rid, 'FORMAT_OBJECT': 'SearchInfo'}
            try:
                status_response = requests.get('https://blast.ncbi.nlm.nih.gov/Blast.cgi', params=check_params)
                status_response.raise_for_status()
            except Exception as e:
                logging.error(f"{name}, ошибка при проверке статуса: {e}")
                raise

            status_text = status_response.text
            if 'Status=READY' in status_text:
                logging.info(f"{name}, результаты готовы.")
                break
            elif 'Status=FAILED' in status_text:
                logging.error(f"{name}, Ошибка: запрос не удался.")
                raise Exception(f"{name}, Ошибка: запрос не удался.")
            elif 'Status=UNKNOWN' in status_text:
                logging.error(f"{name}, Ошибка: RID не найден.")
                raise Exception(f"{name}, Ошибка: RID не найден.")
            else:
                logging.info(f"{name}, Ожидание завершения анализа...")
                time.sleep(5)

        result_params = {'CMD': 'Get', 'RID': rid, 'FORMAT_TYPE': 'Text', 'FORMAT_OBJECT': 'Alignment', 'HITLIST_SIZE': 1}

        try:
            result_response = requests.get('https://blast.ncbi.nlm.nih.gov/Blast.cgi', params=result_params)
            result_response.raise_for_status()
        except Exception as e:
            logging.error(f"{name}, Ошибка при получении результата: {e}")
            raise

        match = re.search(r'^> *([\w|.]+) (.+)', result_response.text, re.MULTILINE)
        subject_id = match.group(1) if match else 'N/A'
        organism = match.group(2).strip() if match else 'N/A'

        evalue = re.search(r'Expect *= *([\deE.-]+)', result_response.text)
        identity = re.search(r'Identities = *\d+/\d+ *\((\d+%)\)', result_response.text)
        align_len = re.search(r'Identities = *(\d+)/(\d+)', result_response.text)

        with open(table_path, 'a') as table_file:
            table_file.write(
                f"{name:<45}\t"
                f"{seq_id:<50}\t"
                f"{subject_id[:15] if subject_id else 'N/A':<15}\t"
                f"{organism[:50] if organism else 'N/A':<50}\t"
                f"{evalue.group(1) if evalue else 'N/A':<10}\t"
                f"{identity.group(1) if identity else 'N/A':<10}\t"
                f"{align_len.group(2) if align_len else 'N/A':<10}\n"
            )
