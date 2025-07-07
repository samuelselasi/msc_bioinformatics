from Bio import pairwise2
from Bio.pairwise2 import format_alignment
from Bio import SeqIO

# Load actual sequences
seq1 = str(SeqIO.read("sequence_1.fasta", "fasta").seq)
seq2 = str(SeqIO.read("sequence_2.fasta", "fasta").seq)

# Perform global alignment (Needleman-Wunsch)
alignments = pairwise2.align.globalxx(seq1, seq2)  # xx: match=1, mismatch=0

# Print the best alignment
print(format_alignment(*alignments[0]))
