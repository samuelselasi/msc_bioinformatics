from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
from Bio.Align.AlignInfo import SummaryInfo

# Read the alignment using the modern API
alignment = AlignIO.read("all_sequences.aln", "clustal")

# Generate consensus using SummaryInfo
summary = SummaryInfo(alignment)
consensus = summary.dumb_consensus(ambiguous='X')  # 'X' for uncertain bases

# Prepare the output text
output_text = []
output_text.append(str(alignment))
output_text.append(f"\nAlignment length: {alignment.get_alignment_length()}")
output_text.append(f"\nConsensus Sequence:\n{consensus}")

# Save to file
with open("msa_summary.txt", "w") as f:
    f.write("\n".join(output_text))

# Optional: also print to console
print("MSA summary saved to msa_summary.txt")
print("\n".join(output_text))
