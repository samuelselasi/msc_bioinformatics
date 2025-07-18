# Lab 2
## Sequence Alignment with Biopython

### ‍Project Overview

This assignment explores **pairwise and multiple sequence alignment (MSA)** techniques using **Biopython** and **ClustalW**. Real-world gene sequences from the hemoglobin beta (HBB) gene family were retrieved from **NCBI** for three species: **human**, **mouse**, and **rat**. The project involves:

- Pairwise alignment using Biopython’s `pairwise2`
- MSA using ClustalW
- Consensus sequence generation
- Visualization of alignment results
- Use of substitution matrices (BLOSUM62)
- Both global (Needleman-Wunsch) and local (Smith-Waterman) alignment

---

### Files and Structure

```bash
.
├── fetch_sequences.py             # Script to fetch FASTA sequences from NCBI
├── pairwise_alignment.py          # Performs pairwise global alignment (globalxx)
├── pairwise_alignment_blosum.py   # Performs global alignment using BLOSUM62
├── local_alignment.py             # Performs local alignment (Smith-Waterman)
├── msa_visualize.py               # Loads MSA results and prints consensus
├── human.fasta                    # Human HBB gene (NM_000518.5)
├── mouse.fasta                    # Mouse HBB gene (NM_008220.2)
├── rat.fasta                      # Rat HBB gene (NM_013096.2)
├── all_sequences.fasta            # Combined FASTA for MSA
├── all_sequences.aln              # ClustalW alignment result
├── all_sequences.dnd              # ClustalW guide tree
├── psiblast_remote.py             # Runs remote PSI-BLAST via NCBI
├── psiblast_result.xml            # PSI-BLAST results (XML format)
├── 1A3N.pdb                       # Hemoglobin structure file (from RCSB)
├── structure_parser.py            # Parses PDB file and prints structure hierarchy
└── README.md                      # This documentation file
```

---

### Biological Sequences

The sequences used in this project were retrieved using NCBI accession numbers:

- **Human**: `NM_000518.5`
- **Mouse**: `NM_008220.2`
- **Rat**: `NM_013096.2`

They were stored in FASTA format for both pairwise and multiple sequence alignment tasks.

---

### Setup Instructions (Ubuntu)

1. **Install Python and Dependencies**

```bash
sudo apt update
sudo apt install python3 python3-pip clustalw
pip3 install biopython
```

2. **Run Scripts**

- Fetch Sequences:  
  ```bash
  python3 fetch_sequences.py
  ```

- Pairwise Alignment (simple global):
  ```bash
  python3 pairwise_alignment.py
  ```

- Pairwise Alignment using BLOSUM62:
  ```bash
  python3 pairwise_alignment_blosum.py
  ```
## Output
```
----ACA-TTTGCTTCTGACACAA--CTGT-GT-TC--ACTA-GCAACCTCA-AA-CAGACAC-CATGGTGCATC-TGACTCCTGA-G--GAGAAGT-CTGCC-GT-TACTGC-CCTGTGGGGCAA-GGTGAACGTG--GATGAAGTTGGTGGTGAGGCCCTGGGCAGGCTGCTGGTG-GTCTACCCTTGGACCCAGA-GGTT-CTTTGAGTC--CTTTGGG-GATC-TG-TCCACTCCTGA--TGCTG-T-T-ATGGGC-AAC--CCTAAG-GTGAAGGCTC-ATGGCAAGAAAGTGC-T--CGGTGCCTTTAGT--GATGGCCTGGC--TCACCT-GGACAA-CCTCAAGGGCACCTTTGCCA-CACTG-AGTGAGCTGC-ACTGTGACAAGCTGCAC-GTGGATCCTGAGAACTTCAGGCTCCTGGGCAACGTGC-TGG-TC-TG-T-GTGCTGGC-CCATC-ACTT-TGG-CAAAGA-ATT-CACCCCACCAG-TGCA--GGCTGCCTATC-AGAAAG-TGGTGGCTGGT-GTGGCTAA--TGCCCTGGCC-CACAAGTATC-ACTAAGCTCGCT-TTCTTGCTGTCCA-ATT-TCTATT--A-AA-GGTTCCT-T-TGTTCCC-TA-AGTCCAACTA--CTA---AACT-G-GGGGATAT-T-ATGAAGGG-CCTT-GAG-CATCTGGA--TTCTGCCTAAT--AA-AAA--AC-ATTTATT-T-TCATTGCAA--
    ||  ||||||||||     |  |||| || |   |||  ||||||||| || ||||||  ||||||||| | |||   |||| |  ||||||  ||| | || | |||  ||||||||| || |||||||  |  ||||||||||||||||||||||||||||||||||||||  ||||||||||||||||||  || | |||||| |   |||| || || | |  | | || |||   | ||| | | |||||  ||   || ||  |||||||| | ||||||||||||||  |  |  ||||||||    |||||||| |   ||| || |||| | |||||||||||||||||||| | ||  |||||||| | ||||||||||||||||  ||||||||||||||||||||||||||||||||  |   | | || || | |||||||  ||| | ||   ||| | |||  ||| || ||| || | ||||  |||||||| || || ||| ||||||||||  ||||||    |||||||| | |||||||| | ||||||| | |  || || |||  |  ||| |||| |  | || || |  | | || |||| || ||   ||  |  ||    || | | ||||| |  | |||||  | |||| | | |||||  |  || |   | ||  || |||  |  |||||   | ||| | |    
GTTTAC-GTTTGCTTCTG-----ATTCTGTTGTGT-TGACT-TGCAACCTCAGAAACAGACA-TCATGGTGCA-CCTGA---CTGATGCTGAGAAG-GCTG-CTGTCT-CTG-GCCTGTGGGG-AAAGGTGAAC--GCCGATGAAGTTGGTGGTGAGGCCCTGGGCAGGCTGCTGGT-TGTCTACCCTTGGACCCAG-CGG-TACTTTGA-T-AGCTTT-GGAGA-CCT-AT-C-CT-CTG-CCT-CTGCTATCATGGG-TAA-TGCC-AA-AGTGAAGGC-CCATGGCAAGAAAGTG-ATAAC--TGCCTTTA--ACGATGGCCT-G-AATCA-CTTGGAC-AGCCTCAAGGGCACCTTTGCCAGC-CT-CAGTGAGCT-CCACTGTGACAAGCTGCA-TGTGGATCCTGAGAACTTCAGGCTCCTGGGCAA--T--AT-GATCGTGATTGTGCTGG-GCCA-CCAC--CTGGGC-AAG-GATTTCA-CCC-CC-GCTGCACAGGCTGCCT-TCCAG-AAGGTGGTGGCTGG-AGTGGCT--GCTGCCCTGG-CTCACAAGTA-CCACTAAGC-C-C-CTT-TT-CTG--C-TATTGTCTA-TGCACAAAGG-T--TATATG-TCCCCTAGAG---AA--AAACT-GTCAA-TTGTGGGGA-A-ATGATGAA--GACCTTTG-GGCATCT--AGCTT-T---T-ATCTAATAAATGA-TATTTA--CTGTCA-T-C--CC
  Score=517
```

- Local Alignment using Smith-Waterman:
  ```bash
  python3 local_alignment.py
  ```

- Combine Sequences for MSA:
  ```bash
  cat human.fasta mouse.fasta rat.fasta > all_sequences.fasta
  ```

- Run ClustalW:
  ```bash
  clustalw -INFILE=all_sequences.fasta
  ```

- View Alignment and Consensus:
  ```bash
  python3 msa_visualize.py
  ```

---

### Output Summary

- **Pairwise Alignment Scores**:
  - Human vs Mouse: 79
  - Human vs Rat: 49
  - Mouse vs Rat: 40

- **MSA Score**: 6802
- **Alignment Length**: *[Printed by script]*
- **Consensus Sequence**: *[Printed by script]*

---

### Key Concepts

- **Global Alignment**: Needleman-Wunsch algorithm (`pairwise2.globalxx`)
- **Local Alignment**: Smith-Waterman algorithm (`pairwise2.localds`)
- **Substitution Matrices**: BLOSUM62 used via `substitution_matrices.load("BLOSUM62")`
- **Multiple Sequence Alignment**: Performed via ClustalW
- **Consensus Sequence**: Derived using `AlignInfo.SummaryInfo`

---

### References

- [Biopython Documentation](https://biopython.org/wiki/Documentation)
- [NCBI Nucleotide](https://www.ncbi.nlm.nih.gov/nuccore)
- ClustalW2: Thompson et al., Nucleic Acids Res. 1994

---
