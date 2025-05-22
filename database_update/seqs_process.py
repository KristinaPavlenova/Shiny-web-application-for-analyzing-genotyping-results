import re
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Align import PairwiseAligner, Alignment
import glob
import os, os.path
from typing import List, Dict, Any


def extract_hits(blast_text: str) -> List[Dict[str, Any]]:
    hits = []
    entries = re.split(r"^>", blast_text, flags=re.MULTILINE)[1:]
    for entry in entries:
        lines = entry.strip().splitlines()
        if not lines:
            continue
        header = lines[0]
        match_header = re.match(r"([\w|.]+) (.+)", header)
        subject_id = match_header.group(1) if match_header else "N/A"
        organism = match_header.group(2).strip() if match_header else "N/A"

        identities_match = re.search(r"Identities\s*=\s*(\d+)\s*/\s*(\d+)", entry)
        evalue_match = re.search(r"Expect *= *([\deE.-]+)", entry)

        if identities_match and evalue_match:
            matches = int(identities_match.group(1))
            align_len = int(identities_match.group(2))
            identity = round((matches / align_len) * 100, 2)
            evalue = float(evalue_match.group(1).replace("e", "E"))
            hits.append(
                {
                    "subject_id": subject_id,
                    "organism": organism,
                    "identity": identity,
                    "align_len": align_len,
                    "evalue": evalue,
                }
            )
    return hits


def load_sequences(fasta_file: str) -> Dict[str, Seq]:
    sequences = {}
    for record in SeqIO.parse(fasta_file, "fasta"):
        sequences[record.id] = record.seq
    return sequences


def align_and_score(seq1: Seq, seq2: Seq, aligner: PairwiseAligner) -> Dict[str, Any]:
    alignment = aligner.align(seq1, seq2)[0]
    aln_seq1 = alignment.aligned[0]
    aln_seq2 = alignment.aligned[1]
    matches = 0
    for (start1, end1), (start2, end2) in zip(aln_seq1, aln_seq2):
        matches += sum(
            seq1[start1:end1][i] == seq2[start2:end2][i]
            for i in range(min(end1 - start1, end2 - start2))
        )
    identity = matches / max(len(seq1), len(seq2)) * 100
    alignment_length = alignment.shape[1]
    return {
        "Alignment_Score": alignment.score,
        "Identity_Percentage": identity,
        "Alignment_Length": alignment_length,
        "Alignment_Object": alignment,
    }


def align_wo_gaps(seq1: Seq, seq2: Seq, aligner: PairwiseAligner) -> Dict[str, Any]:
    trimmed_seq1 = str(seq1).strip("-")
    trimmed_seq2 = str(seq2).strip("-")
    alignment = aligner.align(trimmed_seq1, trimmed_seq2)[0]
    aln_seq1 = alignment.aligned[0]
    aln_seq2 = alignment.aligned[1]
    matches = 0
    for (start1, end1), (start2, end2) in zip(aln_seq1, aln_seq2):
        for i in range(min(end1 - start1, end2 - start2)):
            if trimmed_seq1[start1 + i] == trimmed_seq2[start2 + i]:
                matches += 1
    identity = matches / min(len(trimmed_seq1), len(trimmed_seq2)) * 100
    alignment_length = sum(end - start for start, end in aln_seq1)
    return {
        "Alignment_Score": alignment.score,
        "Identity_Percentage": identity,
        "Alignment_Length": alignment_length,
        "Alignment_Object": alignment,
    }


def save_alignment_text(alignment: Alignment, output_path: str) -> None:
    alignment_str = alignment.format()
    with open(output_path, "w") as f:
        f.write(alignment_str)


def collect_trimmed_seqs(dir: str, assembly: str, gene: str) -> None:
    fastanames = glob.glob("*_trimmed_best.fasta", root_dir=dir)
    for fasta in fastanames:
        file = os.path.join(dir, fasta)
        name = os.path.basename(fasta).split("_trimmed_best.fasta")[0]
        with open(file, "r") as f:
            record = next(SeqIO.parse(f, "fasta"))
            record.id = "_".join(name.split("_")[:-1])
            with open(f"{assembly}_{gene}.fasta", "a") as out_file:
                SeqIO.write(record, out_file, "fasta")
