from Bio import Entrez
from Bio import SeqIO


ids = ['NM_000188.3', 'NM_001322365.2', 'NM_033496.3', 'NM_001322364.2', 'NM_033498.3', 'NM_001358263.1', 'NM_001322367.1', 'NM_001322366.1', 'NM_033500.2']
Entrez.email= 'samelsekasi@gmail.com'
handle=Entrez.efetch(db='nucleotide', id=ids, rettype='fasta', retmode='text')
gene= SeqIO.parse(handle,'fasta')
all_genes= [i for i in gene]
len(all_genes)

for i in all_genes:
  print(i.description)
  print(len(i.seq))

out=open('./datasets/all_hxks.fasta', 'w')
for i in all_genes:
  SeqIO.write(i,out,'fasta')
out.close()
