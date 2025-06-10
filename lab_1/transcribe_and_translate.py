from Bio import Entrez
from Bio import SeqIO
from Bio.Seq import Seq

for sequence in SeqIO.parse("./datasets/all_hxks.fasta",'fasta'):
  j=sequence
  break

print(j)
print(j.seq.transcribe())
print(j.translate(to_stop=True))
print(j[1:].translate(to_stop=True))
print(j[2:].translate(to_stop=True))
