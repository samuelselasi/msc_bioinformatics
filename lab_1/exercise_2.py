from Bio import Entrez
from Bio import SeqIO


ids = ['KC668274.1','KC668273.1','KC668272.1','KC668271.1', 'KC668270.1','KC668269.1']
Entrez.email= 'samelsekasi@gmail.com'
handle=Entrez.efetch(db='nucleotide', id=ids, rettype='fasta', retmode='text')
gene= SeqIO.parse(handle,'fasta')
all_genes= [i for i in gene]
len(all_genes)

for i in all_genes:
  print(i.description)
  print(len(i.seq))

out=open('./datasets/all_lux.fasta', 'w')
for i in all_genes:
  SeqIO.write(i,out,'fasta')
out.close()
