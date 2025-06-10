from Bio import Entrez
from Bio import SeqIO

Entrez.email = 'samelsekasi@gmail.com'
handle = Entrez.efetch(db = 'nuccore', id = '529158032', rettype = 'fasta')
print(handle.read())
