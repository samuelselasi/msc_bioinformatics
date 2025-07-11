from Bio.PDB import PDBParser

parser = PDBParser()
structure = parser.get_structure("HBB", "1A3N.pdb")

for model in structure:
    print(f"Model: {model.id}")
    for chain in model:
        print(f"  Chain: {chain.id}")
        for residue in chain:
            print(f"    Residue: {residue.resname} - {residue.id}")
        break
    break
