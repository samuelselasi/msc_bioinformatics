# Raw Data

Raw aromatase IC50 data fetched from ChEMBL.

- `bioactivty_data.csv`
  - `molecule_chembl_id`  ChEMBL molecule ID
  - `canonical_smiles`    SMILES string
  - `standard_value`      IC50 in nM

Rebuild with:
```bash
make data-raw
