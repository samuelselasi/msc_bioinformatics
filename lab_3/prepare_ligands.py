from rdkit import Chem
from rdkit.Chem import AllChem
import os
import subprocess

# Make output folder
os.makedirs("ligands_pdbqt", exist_ok=True)

# Paste your SMILES here
smiles_list = [
    "C[C@H](CCC[C@H](C)C(=O)O)[C@H]1CC[C@@H]2[C@@]1(CC[C@H]3C2=CC[C@@H]4[C@@]3(CCC(=O)C4)C)C",
    "C[C@H](CCC[C@H](C)C(=O)O)[C@H]1CC[C@@H]2[C@@]1(CC[C@H]3[C@H]2CCC4=CC(=O)CC[C@]34C)C",
    "C[C@H](CCC[C@@H](C)C(=O)O)[C@H]1CC[C@@H]2[C@@]1(CC[C@H]3[C@H]2CCC4=CC(=O)CC[C@]34C)C",
    "C[C@H](CCC[C@H](C)C(=O)O)[C@H]1CC[C@@H]2[C@@]1(CC[C@H]3C2=CC[C@@H]4[C@@]3(C=CC(=O)C4)C)C",
    "C[C@H](CCCC(C)C(=O)O)[C@H]1CC[C@@H]2[C@@]1(CC[C@H]3[C@H]2CCC4=CC(=O)CC[C@]34C)C",
    "C[C@H](CCC[C@@H](C)C(=O)O)[C@H]1CC[C@@H]2[C@@]1(CC[C@H]3C2=CC[C@@H]4[C@@]3(CCC(=O)C4)C)C",
    "C[C@H](CCCC(C)C(=O)O)[C@H]1CC[C@@H]2[C@@]1(CC[C@H]3C2=CC[C@@H]4[C@@]3(CCC(=O)C4)C)C",
    "C[C@H](CCC[C@H](C)C(=O)O)[C@H]1CC[C@@H]2[C@@]1(CC[C@H]3C2=CC[C@@H]4[C@@]3(CC[C@@H](C4)O)C)C",
    "C[C@H](CCC[C@H](C)C(=O)O)[C@H]1CC[C@@H]2[C@@]1(CC[C@H]3[C@H]2CC[C@@H]4[C@@]3(CCC(=O)C4)C)C",
    "C[C@H](CCCC(C)C(=O)O)[C@H]1CC[C@@H]2[C@@]1(CC[C@H]3[C@H]2CC[C@@H]4[C@@]3(CCC(=O)C4)C)C"
]

for i, smi in enumerate(smiles_list):
    try:
        mol = Chem.MolFromSmiles(smi)
        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol)
        AllChem.MMFFOptimizeMolecule(mol)

        sdf_file = f"ligand_{i+1}.sdf"
        pdbqt_file = f"ligands_pdbqt/ligand_{i+1}.pdbqt"

        writer = Chem.SDWriter(sdf_file)
        writer.write(mol)
        writer.close()

        # Convert SDF to PDBQT using Open Babel
        subprocess.run(["obabel", sdf_file, "-O", pdbqt_file, "--gen3d", "--partialcharge", "gasteiger"])
        os.remove(sdf_file)

        print(f"Saved {pdbqt_file}")
    except Exception as e:
        print(f"Failed to process ligand {i+1}: {e}")
