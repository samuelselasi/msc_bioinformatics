from Bio import Entrez, SeqIO

Entrez.email = "samelsekasi@gmail.com"

ids = [
    ("NM_000518.5", "human.fasta"),   # Human HBB
    ("NM_008220.2", "mouse.fasta"),   # Mouse HBB
    ("NM_013096.2", "rat.fasta")      # Rat HBB
]

for gene_id, filename in ids:
    handle = Entrez.efetch(db="nucleotide", id=gene_id, rettype="fasta", retmode="text")
    record = SeqIO.read(handle, "fasta")
    SeqIO.write(record, filename, "fasta")
