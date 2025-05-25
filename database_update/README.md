# Getting to know your amphipods: Updates to the speCOIdent app for analysis of DNA barcoding results
This repo was created to update database of the speCOIdent app. 
## Required packages:


## Update steps:
### 1. COI from NCBI

### 2. COI and 18S from transcriptome

#### 2.1 Main Pipeline and Checks
здесь нужно расписать, что пайплайн делает

#### 2.2 Test pipeline 

This pipeline consists of **11 sequential steps** for processing 18S sequences. The first two steps are executed in the terminal, while steps 3-11 are run using Python scripts. To read files of pipeline see folder `Shiny-web-application-for-analyzing-genotyping-results/database_update/test_pipeline_18S`. Pipeline is based on the idea that the first e-value hit belongs to the best contig. File "test_pipeline_18S_LabNotebook.pdf" contains descriptions about pipeline results. Example of test pipeline usage and output files are available in `./example`.

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
   `python test_pipeline_step3.py`
4. **Create FASTA files for top 5 sequences** (backup if top sequence fails)  
   `python test_pipeline_step4.py`
5. **Select the best sequence from the top 5**  
   `python test_pipeline_step5.py`
6. **MAFFT + TrimAl alignment + BLAST identity check (≥90% identity)**  
   `python test_pipeline_step6.py`
7. **Local BLAST identity verification**  
   `python test_pipeline_step7.py`
8. **Merge sequences passing checks**  
   `python test_pipeline_step8.py`
9. **Final global alignment check (PairwiseAligner)**  
   `python test_pipeline_step9.py`
10. **Merge all validated sequences**  
    `python test_pipeline_step10.py`
11. **UGENE-based identity check (global PairwiseAligner)**  
    `python test_pipeline_step11.py`

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