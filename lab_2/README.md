# Lab 2
## Sequence Alignment with Biopython

### ‍Project Overview

This assignment explores **pairwise and multiple sequence alignment (MSA)** techniques using **Biopython** and **ClustalW**. Real-world gene sequences from the hemoglobin beta (HBB) gene family were retrieved from **NCBI** for three species: **human**, **mouse**, and **rat**. The project involves:

- Pairwise alignment using Biopython’s `pairwise2`
- MSA using ClustalW
- Consensus sequence generation
- Visualization of alignment results

---

### Files and Structure

```bash
.
├── fetch_sequences.py          # Script to fetch FASTA sequences from NCBI
├── pairwise_alignment.py       # Performs pairwise global alignment
├── msa_visualize.py            # Loads MSA results and prints consensus
├── human.fasta                 # Human HBB gene (NM_000518.5)
├── mouse.fasta                 # Mouse HBB gene (NM_008220.2)
├── rat.fasta                   # Rat HBB gene (NM_013096.2)
├── all_sequences.fasta         # Combined FASTA for MSA
├── all_sequences.aln           # ClustalW alignment result
├── all_sequences.dnd           # ClustalW guide tree
└── README.md                   # This documentation file
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

- Pairwise Alignment:
  ```bash
  python3 pairwise_alignment.py
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
- **Substitution Matrices**: Defaults or BLOSUM62 can be used
- **Multiple Sequence Alignment**: Performed via ClustalW
- **Consensus Sequence**: Derived using `AlignInfo.SummaryInfo`

---

### References

- [Biopython Documentation](https://biopython.org/wiki/Documentation)
- [NCBI Nucleotide](https://www.ncbi.nlm.nih.gov/nuccore)
- ClustalW2: Thompson et al., Nucleic Acids Res. 1994

---

