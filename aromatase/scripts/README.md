# Scripts

Command line scripts used by the Makefile.

## Data Build

- `fetch_chembl.py` → `data/raw/bioactivty_data.csv`
- `build_labels.py` → labels classes and writes `data/processed/bioactivty_data_preprocessed_class.csv` and `bioactivity_data_no_intermediate.csv`
- `make_pic50.py` → computes pIC50 and writes `data/processed/bioactivity_data_3class_pIC50.csv` (drops invalid IC50 rows)
- `morgan_fp.py` → computes RDKit Morgan fingerprints and writes `data/processed/bioactivity_data_descriptors_morgan.csv`

---

## Models
- `train_cls.py` → trains logistic regression and random forest on 2-class data (active vs inactive)
- `train_reg.py` → trains ridge and random forest to predict pIC50

Outputs go to `results/`:
- `results/models/*.joblib`
- `results/figures/*.png`
- `results/tables/*.csv`
- `results/*_metrics.json`

---

## Inference
- `predict.py` → predict class or pIC50 from SMILES
  - Single SMILES: `--smiles "CCO"`
  - CSV with `smiles` column: `--in file.csv`
  - Text file with one SMILES per line: `--in file.txt`
  - Choose model: `--model-name cls_random_forest` or `--model-name reg_random_forest`
  - Optional: `--as-ic50-nm` converts predicted pIC50 to IC50 nM

---

## Protein
- `fetch_pdb.py` → downloads PDB file (default 5JL9) into `data/raw/pdb/`
