from Bio import Entrez
from Bio import SeqIO
from Bio.Seq import Seq
for sequence in SeqIO.parse("./datasets/all_hxks.fasta",'fasta'):
  j=sequence
  break
print(j.id)
print(j.seq)
print(type(j))
