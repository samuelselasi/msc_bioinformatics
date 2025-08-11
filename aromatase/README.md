# Aromatase Bioactivity Pipeline

End-to-end reproducible pipeline to:
1. Fetch aromatase IC50 data from ChEMBL.
2. Clean, label (active, intermediate, inactive), and compute pIC50.
3. Generate RDKit Morgan fingerprints.
4. Train baseline models for classification and regression.
5. Run single or batch inference from SMILES.

---

## Quick start

```
python3 -m venv .venv && source .venv/bin/activate
pip install -r env/requirements.txt

# Build data and descriptors
make data-all

# Train models
make train-cls
make train-reg

# Predict from a SMILES
make predict-cls SMILES="CCO" MODEL=cls_random_forest
make predict-reg SMILES="CCO" MODEL=reg_random_forest
```

## Project Layout

```
data/
  raw/                 # Raw ChEMBL CSV
  processed/           # Labeled, pIC50, and descriptors
  external/            # Optional third-party files (SDF)
env/                   # Environment files
results/
  figures/             # Plots from training
  models/              # Saved joblib models
  tables/              # Reports and feature ranks
  predictions/         # Inference outputs
scripts/               # All pipeline and training scripts
notebooks/             # Optional exploration
```

---

## Reproducibility

* Python virtual environment is recommended.
* `env/requirements.txt` pins NumPy 1.26.4 to match RDKit wheels.
* Makefile uses file targets so only stale steps rebuild.

---

# Data Sources & License

* ChEMBL bioactivity data is fetched at build time. Check ChEMBL terms if you plan to redistribute.
* Optional natural products SDF can be placed in data/external/. See data/external/README.md.

---

## Makefile Summary

Run `make help` for a quick list of targets.

# env

Environment and dependency management.

## Install

```
python3 -m venv .venv && source .venv/bin/activate
pip install -r env/requirements.txt
```

`requirements.txt` pins:
* `numpy==1.26.4` to match `rdkit-pypi==2022.9.5`
* `scikit-learn==1.3.2` and `scipy==1.10.1` for compatibility

---

## Alternative (Conda)

```
conda create -n aromatase -c conda-forge python=3.10 rdkit pandas numpy scikit-learn chembl_webresource_client biopython
conda activate aromatase
```
