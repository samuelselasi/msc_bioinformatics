PY := python
VENV := .venv
PIP := $(VENV)/bin/pip

# data artifacts
RAW := data/raw/bioactivty_data.csv
P3  := data/processed/bioactivty_data_preprocessed_class.csv
PPI := data/processed/bioactivity_data_3class_pIC50.csv
FPS := data/processed/bioactivity_data_descriptors_morgan.csv

# default PDB ID
PDB ?= 5JL9

.PHONY: help install clean data-all data-raw data-processed fps pdb train-cls train-reg predict-cls predict-reg

help:
	@echo "make install        # create venv and install deps"
	@echo "make data-all       # fetch chembl -> labels -> pIC50 -> morgan fps"
	@echo "make data-raw       # build just the raw ChEMBL CSV"
	@echo "make data-processed # build labeled + pIC50 CSVs"
	@echo "make fps            # build Morgan fingerprints only"
	@echo "make pdb PDB=5JL9   # download PDB (default 5JL9)"
	@echo "make train-cls      # train classifier (active vs inactive)"
	@echo "make train-reg      # train regressor (pIC50)"
	@echo "make predict-cls SMILES='CCO' MODEL=cls_random_forest"
	@echo "make predict-reg SMILES='CCO' MODEL=reg_random_forest"

install:
	python3 -m venv $(VENV)
	. $(VENV)/bin/activate && $(PIP) install -r env/requirements.txt

# file-based rules so we only rebuild when inputs change
$(RAW): scripts/fetch_chembl.py
	$(VENV)/bin/python scripts/fetch_chembl.py

$(P3): $(RAW) scripts/build_labels.py
	$(VENV)/bin/python scripts/build_labels.py

$(PPI): $(P3) scripts/make_pic50.py
	$(VENV)/bin/python scripts/make_pic50.py

$(FPS): $(PPI) scripts/morgan_fp.py
	$(VENV)/bin/python scripts/morgan_fp.py

# user-friendly aggregations
data-raw: $(RAW)

data-processed: $(PPI)

fps: $(FPS)

data-all: $(FPS)

pdb:
	$(VENV)/bin/python scripts/fetch_pdb.py $(PDB)

# training
train-cls: $(FPS) scripts/train_cls.py
	$(VENV)/bin/python scripts/train_cls.py

train-reg: $(FPS) scripts/train_reg.py
	$(VENV)/bin/python scripts/train_reg.py

# inference helpers (single-SMILES convenience)
# Example: make predict-cls SMILES="CCO" MODEL=cls_random_forest
predict-cls: $(FPS) scripts/predict.py
	@if [ -z "$(SMILES)" ]; then echo 'Usage: make predict-cls SMILES="CCO" MODEL=cls_random_forest'; exit 1; fi
	$(VENV)/bin/python scripts/predict.py --task cls --smiles "$(SMILES)" --model-name $${MODEL:-cls_random_forest}

# Example: make predict-reg SMILES="CCO" MODEL=reg_random_forest
predict-reg: $(FPS) scripts/predict.py
	@if [ -z "$(SMILES)" ]; then echo 'Usage: make predict-reg SMILES="CCO" MODEL=reg_random_forest'; exit 1; fi
	$(VENV)/bin/python scripts/predict.py --task reg --smiles "$(SMILES)" --model-name $${MODEL:-reg_random_forest}

clean:
	rm -rf data/processed/*.csv data/raw/*.csv data/raw/pdb results/*
