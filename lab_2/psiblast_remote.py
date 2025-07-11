from Bio.Blast import NCBIWWW

# PSI-BLAST parameters for nucleotide
result_handle = NCBIWWW.qblast(
    program="blastn", 
    database="nt", 
    sequence=open("human.fasta").read(),
    entrez_query="hemoglobin",
    expect=0.01,
    hitlist_size=10
)

with open("psiblast_result.xml", "w") as out_handle:
    out_handle.write(result_handle.read())

print("PSI-BLAST search completed. Results saved to psiblast_result.xml")
