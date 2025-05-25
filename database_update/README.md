# Getting to know your amphipods: Updates to the speCOIdent app for analysis of DNA barcoding results
This repo was created to update database of the speCOIdent app. Updating was performed by adding to database cytochrome oxidase I (COI) sequence data from the mitochondrial genomes of Baikal amphipods published in [NCBI](https://www.ncbi.nlm.nih.gov) and COI and 18S ribosomal subunit sequence data from two published amphipod transcriptome [assemblies](https://doi.org/10.1186/s12862-021-01806-9), Trinity and rnaSPAdes.

The use of the pipelines presented below for extracting sequences from transcriptome assemblies is also possible for other transcriptome results.

## Update steps:

### 1. COI from NCBI

For 14 species COI sequences from the mitochondrial genome published in NCBI were added to the database with geographic coordinates: *Crypturopus tuberculatus*, *Gmelinoides fasciatus*, *Linevichella vortex*, *Brachyuropus grewingkii*, *Eulimnogammarus cyaneus*, *Eulimnogammarus verrucosus*, *Eulimnogammarus vittatus*, *Garjajewia cabanisii*, *Acanthogammarus victorii*, *Pallaseopsis kessleri*, *Macrohectopus branickii*, *Echiuropus macronychus*, *Baikalogammarus pullus*, *Crypturopus inflatus*.


### 2. COI and 18S from transcriptome assemblies

Two pipelines were developed: main pipeline for COI and 18S assembly results processing, and additional test pipeline for 18S assembly results.

#### 2.1 Main Pipeline and Checks

This pipeline consists of **5 sequential steps** for processing COI and 18S sequences and 3 checks of obtained genes (the last one for comparison sequences from different assemblies). All processes are run using Python scripts. To read files of pipeline see folder `Shiny-web-application-for-analyzing-genotyping-results/database_update/main_pipeline`. File `LabNotebook.pdf` contains descriptions about main pipeline results. Example of main pipeline usage is presented in section below.

**Requirements:**
- Working directory must contain:
  - Reference file: `./references/reference.fasta`
  - Folder with transcriptome assemblies: `./input_fasta`
  - All pipeline scripts (`config.py`, `main_pipeline_step1.py` to `main_pipeline_step5.py`, `check1.py` to `check3.py`, `modules/seqs_process.py`)
- Packages and tools are available in file `main_pipeline_requirements.txt`

**Pipeline Steps:**

The file `config.py` must specify the path to the directory with the input data (`input_dir`) and to the reference fasta (`reference`), the path to file for recording logs (`log_file`), the variant of assembly (`assembly = "rnaspades"` or `assembly = "trinity"`), the gene (`gene = "COI"` or `gene = "18S"`), the number of best hits for step 2 (`top_n`, default `top_n = 5`), the indexes (sequence number) of the input fasta files for cases when not all files in the directory need to be processed (`idxs`) and the path to directory for saving the fasta file with trimmed sequences (`dir_final_fasta`).

1. **BLAST alignment of reference to BLAST database with transcriptome results**
```bash
python main_pipeline_step1.py
```

2. **Hits filtration by length (+ coverage for rnaspades) and top hits selection by e-value (+ coverage for rnaspades)**
```bash
python main_pipeline_step2.py
```

3. **Global multiple alignment of top hits to reference**
```bash
python main_pipeline_step3.py
```

4. **Top1 selection by reference identity (+ coverage for rnaspades)**
```bash
python main_pipeline_step4.py
```

5. **Top1 trimming by reference gene boundaries**
```bash
python main_pipeline_step5.py
```

6. **Check 1 - trimmed top1 sequences filtration by length, ‘NNN’ and gaps content**
```bash
python check1.py
```

During this check all exceptions are collected in the file `./main_pipeline_results/top1_check/seqs_len_exceptions_{assembly}_{gene}.txt`. The following length restrictions apply to genes:

- COI: 1531 <= length <= 1538
- 18S: 2260 <= length <= 2274

After this check pipeline results for exceptions were reprocessed (only some steps) using the following approaches:

- choosing a suitable by lenght hit for trimming
- moving one nucleotide in multiple alignment results
- multiple alignment of the best coverage hit
- multiple alignment without the worst coverage hits
- multiple alignment without the hits containing ‘NNN’
- `top_n = 20` for step 2

Detailed instructions are contained in the file `LabNotebook.pdf`.

7. **Check 2 - the best hits by e-value and identity of trimmed top1 sequences BLAST to NCBI database**
```bash
python check2.py
```

During this check results are collected in the file `./main_pipeline_results/top1_check/blast_summary_{assembly}_{gene}.txt`.

8. **Check 3 - comparison of trimmed top1 sequences from different assemblies by identity and length**
```bash
python check3.py
```

During this check results of alignment of gene sequences obtained for different assemblies are collected in the directory `./main_pipeline_results/comparison_{gene}`, results of comparison - in the file  `./main_pipeline_results/comparison_summary_{gene}.txt`


**Output:**
- Final files containing the found sequences: `./main_pipeline_results/final_seqs/{assembly}_{gene}.fasta` 
- Results:
  - Processed 65 sequences of each assembly and each gene for the Baikal amphipod species
  - Some COI sequences showed high identity to another species or low identity to the same during check 3
  - Many 18S sequences showed high similarity to reference during check 3


**Post-processing:**
- Analysis of the results of check 3 for assembly selection for sequences that had different lengths and low identity when compared (see `Shiny-web-application-for-analyzing-genotyping-results/database_update/main_pipeline/LabNotebook.pdf`)
- Gaps removing in sequences in final fasta files 
- Condensing information to the file `Shiny-web-application-for-analyzing-genotyping-results/database_update/seqs_to_database.csv` for sequences that can be added from transcriptome assemblies to app database.
- Updating database


**Notes:**
- All obtained sequences passed all checks and COI sequences that have not passed check 3 are collected in the directory `Shiny-web-application-for-analyzing-genotyping-results/database_update/final_seqs` 
- Files with results of step 2 are collected in the directory `Shiny-web-application-for-analyzing-genotyping-results/database_update/top_hits`
- Files with checks 1-3 results are collected in the directory `Shiny-web-application-for-analyzing-genotyping-results/database_update/top1_check`


#### 2.2 Test pipeline 

This pipeline consists of **11 sequential steps** for processing 18S sequences. The first two steps are executed in the terminal, while steps 3-11 are run using Python scripts. To read files of pipeline see folder `Shiny-web-application-for-analyzing-genotyping-results/database_update/test_pipeline_18S`. Pipeline is based on the idea that the first e-value hit belongs to the best contig. File `test_pipeline_18S_LabNotebook.pdf` contains descriptions about pipeline results. Example of test pipeline usage and output files are available in `./example`.

**Requirements:**
- Working directory must contain:
  - Reference file: `18S_reference.fasta`
  - Folder with transcriptome assemblies: `./fastas`
  - All pipeline scripts (`test_pipeline_step1.txt` to `test_pipeline_step11.py`)
- Packages and tools are available in file `test_pipeline_requirements.txt`

**Pipeline Steps:**
1. **Make BLAST databases**  
   Execute commands from `test_pipeline_step1.txt`
2. **Align to reference**  
   Execute commands from `test_pipeline_step2.txt`
3. **Extract top 5 BLAST hits by e-value**  
```bash
python test_pipeline_step3.py
```
4. **Create FASTA files for top 5 sequences** (backup if top sequence fails)  
```bash
python test_pipeline_step4.py
```
5. **Select the best sequence from the top 5**  
```bash
python test_pipeline_step5.py
```
6. **MAFFT + TrimAl alignment + BLAST identity check (≥90% identity)**  
```bash
python test_pipeline_step6.py
```
7. **Local BLAST identity verification**  
```bash
python test_pipeline_step7.py
```
8. **Merge sequences passing checks**  
```bash
python test_pipeline_step83.py
```
9. **Final global alignment check (PairwiseAligner)**  
```bash
python test_pipeline_step9.py
```
10. **Merge all validated sequences**  
```bash
python test_pipeline_step10.py
```
11. **UGENE-based identity check (global PairwiseAligner)**  
```bash
python test_pipeline_step11.py
```

**Output:**
- Final file: `final_file.fasta` (headers matched to `./geographic_coordinates` table)
- Results:
  - Processed 25 sequences (rnaspades) and 18 (trinity)
  - After filtering, retained **26 sequences**
  - Many showed high similarity to reference via NCBI BLAST

**Post-processing:**
- Performed multiple contig alignment (UGENE) of 26 best hits + reference
- Trimmed at reference boundaries
- Identity comparison files in `./final_seqs`:
  - `union_nogaps.fasta`
  - `18S_ugene.fasta`
- Identity % calculated via `test_pipeline_step11.py`

**Notes:**
- Adjust file paths in scripts if needed
- Ensure dependencies are installed (see `test_pipeline_requirements.txt`)


## Example of main pipeline usage

Download ```main_pipeline``` directory to your computer

🚀 Execute the following commands in the terminal:

```bash
python main_pipeline_step1.py
python main_pipeline_step2.py
python main_pipeline_step3.py
python main_pipeline_step4.py
python main_pipeline_step5.py
python check1.py
python check2.py
```

📝 Edit the following lines in `config.py`:

```python
# Line 3:
input_dir = os.path.join("example", "input_fasta_trinity")

# Line 8:
assembly = "trinity"
```

🚀 Execute the following commands in the terminal:

```bash
python main_pipeline_step1.py
python main_pipeline_step2.py
python main_pipeline_step3.py
python main_pipeline_step4.py
python main_pipeline_step5.py
python check1.py
python check2.py
python check3.py
```
You must obtain the files, as in ```main_pipeline_results_example``` directory

## System requirements
- Python 3.x

## 👥 Authors

The pipelines was developed by:

- **Kristina Pavlenova**
  📧 Email: [pavlenova.kristina@gmail.com]  
  🏢 Affiliation: 
  - BIOCAD, 198515, St. Petersburg, Russia
  - Bioinformatics Institute, Kantemirovskaya st. 2A, 197342, St. Petersburg, Russia

- **Valeria Afanasyeva**
  📧 Email: [digitalvaleria0.0@gmail.com]  
  🏢 Affiliation: 
  - I.M. Sechenov First Moscow State Medical University, 119048, Moscow, Russia
  - Bioinformatics Institute, Kantemirovskaya st. 2A, 197342, St. Petersburg, Russia

- **Polina Drozdova** 
  🏢 Affiliation: 
  - Irkutsk State University, 664003, Irkutsk, Russia
