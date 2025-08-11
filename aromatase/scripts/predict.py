#!/usr/bin/env python3
import argparse, sys
from pathlib import Path
import numpy as np
import pandas as pd
from joblib import load
from rdkit import Chem
from rdkit.Chem import AllChem

ROOT = Path(".")
RESULTS = ROOT / "results"
MODELS = RESULTS / "models"
PRED = RESULTS / "predictions"
for d in [RESULTS, MODELS, PRED]:
    d.mkdir(parents=True, exist_ok=True)

def read_smiles(args) -> pd.DataFrame:
    # 1) direct string
    if args.smiles:
        return pd.DataFrame([{"id": "mol0", "smiles": args.smiles}])
    # 2) file input
    if not args.infile:
        print("Error: provide --smiles or --in <file>", file=sys.stderr)
        sys.exit(2)
    p = Path(args.infile)
    if not p.exists():
        print(f"Error: file not found: {p}", file=sys.stderr)
        sys.exit(2)
    if p.suffix.lower() == ".csv":
        df = pd.read_csv(p)
        for col in ["smiles", "canonical_smiles", "SMILES"]:
            if col in df.columns:
                df = df.rename(columns={col: "smiles"})
                if "id" not in df.columns:
                    df["id"] = df.index.astype(str)
                return df[["id", "smiles"]]
        print(f"Error: CSV must have a 'smiles' or 'canonical_smiles' column. Got: {list(df.columns)}", file=sys.stderr)
        sys.exit(2)
    else:
        # assume newline-delimited txt
        rows = []
        with open(p, "r") as fh:
            for i, line in enumerate(fh):
                s = line.strip()
                if s:
                    rows.append({"id": f"mol{i}", "smiles": s})
        return pd.DataFrame(rows)

def morgan_df(smiles, radius=2, n_bits=2048):
    fps, ok = [], []
    for s in smiles:
        m = Chem.MolFromSmiles(s)
        if m:
            fp = AllChem.GetMorganFingerprintAsBitVect(m, radius, nBits=n_bits)
            fps.append(list(fp))
            ok.append(True)
        else:
            fps.append([0]*n_bits)
            ok.append(False)
    X = pd.DataFrame(fps, columns=[f"FP_{i}" for i in range(n_bits)])
    return X, pd.Series(ok, name="parsed_ok")

def load_model(task, model_name=None, model_path=None):
    if model_path:
        return load(model_path), Path(model_path).stem
    # sensible defaults
    default = {"cls": "cls_random_forest.joblib", "reg": "reg_random_forest.joblib"}[task]
    fname = model_name + ".joblib" if model_name else default
    path = MODELS / fname
    if not path.exists():
        print(f"Error: model not found at {path}", file=sys.stderr)
        sys.exit(2)
    return load(path), path.stem

def predict_cls(model, X, thr=0.5):
    prob = model.predict_proba(X)[:, 1]
    lab = (prob >= thr).astype(int)
    lab_txt = np.where(lab == 1, "active", "inactive")
    return prob, lab_txt

def predict_reg(model, X, as_ic50_nm=False):
    pic50 = model.predict(X)
    if as_ic50_nm:
        # IC50[nM] = 10**(9 - pIC50)
        ic50 = np.power(10.0, 9.0 - pic50)
        return pic50, ic50
    return pic50, None

def main():
    ap = argparse.ArgumentParser(description="Predict activity class or pIC50 from SMILES using saved models.")
    ap.add_argument("--task", choices=["cls", "reg"], required=True, help="Prediction task: classification or regression")
    ap.add_argument("--model-name", help="Model name without .joblib (e.g., cls_random_forest, cls_logreg, reg_random_forest, reg_ridge)")
    ap.add_argument("--model-path", help="Explicit path to a .joblib model")
    ap.add_argument("--smiles", help="Single SMILES string")
    ap.add_argument("--in", dest="infile", help="Input file (.csv with 'smiles' column or newline .txt)")
    ap.add_argument("--threshold", type=float, default=0.5, help="Classification threshold (default 0.5)")
    ap.add_argument("--as-ic50-nm", action="store_true", help="For regression: also output IC50 in nM")
    ap.add_argument("--out", help="Output CSV path (default results/predictions/pred_<task>_<model>.csv)")
    args = ap.parse_args()

    df_in = read_smiles(args)
    X, ok = morgan_df(df_in["smiles"])
    model, tag = load_model(args.task, args.model_name, args.model_path)

    if args.task == "cls":
        prob, lab = predict_cls(model, X, args.threshold)
        df_out = df_in.copy()
        df_out["parsed_ok"] = ok
        df_out["prob_active"] = prob
        df_out["pred_class"] = lab
        out = args.out or (PRED / f"pred_cls_{tag}.csv")
    else:
        pic50, ic50 = predict_reg(model, X, args.as_ic50_nm)
        df_out = df_in.copy()
        df_out["parsed_ok"] = ok
        df_out["pred_pIC50"] = pic50
        if ic50 is not None:
            df_out["pred_IC50_nM"] = ic50
        out = args.out or (PRED / f"pred_reg_{tag}.csv")

    df_out.to_csv(out, index=False)
    print(f"Wrote {out}  rows={len(df_out)}")
    # also show a small preview
    print(df_out.head(10).to_string(index=False))

if __name__ == "__main__":
    main()
