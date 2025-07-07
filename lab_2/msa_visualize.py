from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
from Bio.Align.AlignInfo import SummaryInfo

# Read the alignment using the modern API
alignment = AlignIO.read("all_sequences.aln", "clustal")

# Generate consensus using SummaryInfo
summary = SummaryInfo(alignment)
consensus = summary.dumb_consensus(ambiguous='X')  # 'X' for uncertain bases

print(alignment)
print(f"\nAlignment length: {alignment.get_alignment_length()}")
print(f"\nConsensus Sequence:\n{consensus}")
