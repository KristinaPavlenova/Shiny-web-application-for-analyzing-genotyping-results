from Bio import SeqIO
from Bio.Align import PairwiseAligner


def load_sequence(fasta_file):
    record = next(SeqIO.parse(fasta_file, "fasta"))
    return record.seq


def compute_identity(seq1, seq2):
    aligner = PairwiseAligner()
    aligner.mode = "global"
    alignment = aligner.align(seq1, seq2)[0]

    aligned_seq1 = alignment.aligned[0]
    aligned_seq2 = alignment.aligned[1]

    matches = 0
    for (start1, end1), (start2, end2) in zip(aligned_seq1, aligned_seq2):
        matches += sum(
            seq1[start1 + i] == seq2[start2 + i]
            for i in range(min(end1 - start1, end2 - start2))
        )

    total_length = max(len(seq1), len(seq2))
    identity = matches / total_length * 100
    return identity


# Пример использования
fasta1 = "18S_reference.fasta"
fasta2 = "fin_file.fasta"

seq1 = load_sequence(fasta1)
seq2 = load_sequence(fasta2)

identity = compute_identity(seq1, seq2)
print(f"Процент идентичности: {identity:.2f}%")
