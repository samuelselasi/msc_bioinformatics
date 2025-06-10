from Bio import Entrez
from Bio import SeqIO

Entrez.email= 'samelsekasi@gmail.com'
handle=Entrez.efetch(db='nucleotide', id=' KC668274.1', rettype='gb')
print(handle.read())
