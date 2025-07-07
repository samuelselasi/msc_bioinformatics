from Bio import pairwise2
from Bio.pairwise2 import format_alignment
from Bio import SeqIO
from Bio.Align import substitution_matrices

# Load sequences (same human and mouse example)
seq1 = str(SeqIO.read("human.fasta", "fasta").seq)
seq2 = str(SeqIO.read("mouse.fasta", "fasta").seq)

# Load BLOSUM62 matrix
matrix = substitution_matrices.load("BLOSUM62")

# Gap penalties
gap_open = -10
gap_extend = -0.5

# Local alignment using Smith-Waterman algorithm
alignments = pairwise2.align.localds(seq1, seq2, matrix, gap_open, gap_extend)

# Display the best local alignment
print(format_alignment(*alignments[0]))
