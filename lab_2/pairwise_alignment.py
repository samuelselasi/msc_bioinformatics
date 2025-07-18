from Bio import pairwise2
from Bio.pairwise2 import format_alignment
from Bio import SeqIO

# Load actual sequences
seq1 = str(SeqIO.read("sequence_1.fasta", "fasta").seq)
seq2 = str(SeqIO.read("sequence_2.fasta", "fasta").seq)

# Perform global alignment (Needleman-Wunsch)
alignments = pairwise2.align.globalxx(seq1, seq2)  # match=1, mismatch=0

# Format the best alignment
best_alignment = format_alignment(*alignments[0])

# Save to file
with open("alignment_output.txt", "w") as f:
    f.write(best_alignment)

# Optional: also print to console
print("Alignment saved to alignment_output.txt")
print(best_alignment)
