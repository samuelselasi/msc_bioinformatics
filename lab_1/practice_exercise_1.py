from Bio import Entrez, SeqIO, AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML

Entrez.email = 'samelsekasi@gmail.com'
handle = Entrez.efetch(db = 'nuccore', id = 'AB000824.1', rettype = 'fasta')
# print(handle.read())

gene = SeqIO.read(handle, 'fasta')
handle.close()

with open("./datasets/p1_gene_sequence.fasta", "w") as output_file:
    SeqIO.write(gene, output_file, "fasta")

# print(f"Downloaded sequence:\n{gene.seq}")

dna_seq = gene.seq
mrna_seq = dna_seq.transcribe()
# print(f"mRNA sequence:\{mrna_seq}")

for frame in range(3):
    protein = mrna_seq[frame:].translate(to_stop=False)
    # print(f"\nProtein translation (frame {frame+1}):\n{protein}")

forward_primer = dna_seq[:30]
reverse_primer = dna_seq[-30:].reverse_complement()

# print(f"Forward Primer (5'->3'): {forward_primer}")
# print(f"Reverse Primer (5'->3'): {reverse_primer}")


handle = Entrez.efetch(db="nuccore", id="AB000824.1", rettype="gb")
gb_gene = SeqIO.read(handle, "genbank")
handle.close()

protein_id = None
for feature in gb_gene.features:
    if feature.type == "CDS" and "protein_id" in feature.qualifiers:
        protein_id = feature.qualifiers["protein_id"][0]
        break

# print(f"Protein ID: {protein_id}")


handle = Entrez.efetch(db="protein", id=protein_id, rettype="fasta")
protein_record = SeqIO.read(handle, "fasta")
handle.close()

with open("./datasets/p1_protein_sequence.fasta", "w") as output_file:
    SeqIO.write(protein_record, output_file, "fasta")

# print(f"Protein sequence:\n{protein_record.seq}")

domain_seq = protein_record.seq[47:230]
# print(f"Catalytic domain sequence:\n{domain_seq}")

domain_record = SeqRecord(
    domain_seq,
    id="Catalytic_Domain",
    description="Glycoside hydrolase family 37 (IPR001661)"
)

with open("./datasets/catalytic_domain.fasta", "w") as f:
    SeqIO.write(domain_record, f, "fasta")

domain_seq = """YQDDKQFVDMPLSIAPEQVLQTFTELSRDHNHSIPREQLQAFVHEHFQAKGQELQPWTPADWKDSPQFLQKISDAKLRAWAGQLHQLWKKLGKKMKPEVLSHPERFSLIYSEHPFIVPGGRFVEFYYWDSYWVMEGLLLSEMAETVKGMLQNFLDLVKTYGHVPNGGRVYYLQRSQPPLLTLM"""

fasta_sequence = f">domain_seq\n{domain_seq}"

result_handle = NCBIWWW.qblast("blastp", "nr", fasta_sequence)

with open("./datasets/blast_result.xml", "w") as out_handle:
    out_handle.write(result_handle.read())
result_handle.close()


with open("./datasets/blast_result.xml") as result_handle:
    blast_record = NCBIXML.read(result_handle)

homologs = []
for alignment in blast_record.alignments[:5]:
    homologs.append((alignment.hit_id, alignment.hit_def))

'''print("Top 5 Homologs:")
for h in homologs:
    print(h)'''

for i, (acc_id, description) in enumerate(homologs):
    organism = description.split("[")[-1].replace("]", "")
    print(f"{organism}: {acc_id}")


ids = [h[0].split('|')[1] for h in homologs]
handle = Entrez.efetch(db="protein", id=ids, rettype="fasta", retmode="text")
records = list(SeqIO.parse(handle, "fasta"))
handle.close()

SeqIO.write(records, "./datasets/homologs.fasta", "fasta")


alignment = AlignIO.read("./datasets/aligned.fasta", "clustal")

print(f"Number of sequences: {len(alignment)}")
print(f"Alignment length: {alignment.get_alignment_length()}")
for record in alignment:
    print(record.id)
