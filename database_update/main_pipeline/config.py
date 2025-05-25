import os, os.path

input_dir = os.path.join("example", "input_fasta_rnaspades")  # указать директорию входных fasta-файлов
log_file = os.path.join("main_pipeline_results", "main_pipeline.log")  # указать файл для логов
reference = os.path.join(
    "example", "references", "Eulimnogammarus_verrucosus_COI.fasta"
)  # указать fasta-файл с референсом
assembly = "rnaspades"  # указать вариант сборки (rnaspades или trinity)
gene = "COI"  # указать ген для поиска (COI или 18S)
top_n = 5  # указать количество лучших результатов после blast-выравнивания
idxs = list(
    range(1)
)  # указать индексы (номер по порядку) входных fasta-файлов для пайплайна
dir_final_fasta = os.path.join("main_pipeline_results", "final_seqs")  # указать директорию для сохранения fasta-файла с обрезанными последовательностями для всех видов
