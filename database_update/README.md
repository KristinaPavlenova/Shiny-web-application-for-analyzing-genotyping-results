# Getting to know your amphipods: Updates to the speCOIdent app for analysis of DNA barcoding results
This repo was created to update database of the speCOIdent app. 
## Required packages:


## Update steps:
### 1. COI from NCBI

### 2. COI and 18S from transcriptome

#### Main Pipeline and Checks
здесь нужно расписать, что пайплайн делает

##### Example

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



#### Test pipeline


## Example of main pipeline usage
