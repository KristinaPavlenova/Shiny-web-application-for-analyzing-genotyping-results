import os
from Bio import SeqIO
from Bio.Align import PairwiseAligner


def load_single_sequence(fasta_file):
    return next(SeqIO.parse(fasta_file, "fasta")).seq


def calculate_identity(seq1, seq2):
    aligner = PairwiseAligner()
    aligner.mode = "global"
    alignment = aligner.align(seq1, seq2)[0]
    aln_seq1 = alignment.aligned[0]
    aln_seq2 = alignment.aligned[1]

    matches = 0
    for (start1, end1), (start2, end2) in zip(aln_seq1, aln_seq2):
        for i in range(min(end1 - start1, end2 - start2)):
            if seq1[start1 + i] == seq2[start2 + i]:
                matches += 1

    identity = matches / max(len(seq1), len(seq2)) * 100  # % of identity
    return identity, alignment


def process_fasta(input_file, reference_seq, min_identity=90.0):
    sequences_to_save = []

    for seq_record in SeqIO.parse(input_file, "fasta"):
        seq = seq_record.seq
        identity, alignment = calculate_identity(seq, reference_seq)

        print(f"{seq_record.id}: identity {identity:.1f}%")

        if identity >= min_identity:
            # print(f"  - saved (identity > {min_identity}%)")
            sequences_to_save.append(seq_record)
        else:
            print(f"  - too low identity (< {min_identity}%)")

    if sequences_to_save:
        with open("filtered_sequences.fasta", "w") as output_file:
            SeqIO.write(sequences_to_save, output_file, "fasta")
        print(
            f"\nSaved {len(sequences_to_save)} seqs in file 'filtered_sequences.fasta'."
        )
    else:
        print("\nNo seqs passed checking.")


def main():
    input_fasta = "fin_file.fasta"
    reference_fasta = "18S_reference.fasta"
    reference_seq = load_single_sequence(reference_fasta)
    process_fasta(input_fasta, reference_seq, min_identity=90.0)


if __name__ == "__main__":
    main()
