from Bio import Entrez
from Bio import SeqIO

Entrez.email= 'samelsekasi@gmail.com'
handle=Entrez.efetch(db='nucleotide', id=' KC668274.1', rettype='fasta', retmode='text')
gene= SeqIO.read(handle,'fasta')
print(gene)
print(gene.id)
print(gene.description)
print(gene.seq)
print(len(gene.seq))
print(gene.features)
