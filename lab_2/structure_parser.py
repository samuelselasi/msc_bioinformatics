from Bio.PDB import PDBParser

# Initialize the parser and structure
parser = PDBParser()
structure = parser.get_structure("HBB", "1A3N.pdb")

# Prepare output list
output_lines = []

for model in structure:
    output_lines.append(f"Model: {model.id}")
    for chain in model:
        output_lines.append(f"  Chain: {chain.id}")
        for residue in chain:
            output_lines.append(f"    Residue: {residue.resname} - {residue.id}")
        break  # Only the first chain per model
    break  # Only the first model

# Write to file
with open("structure_summary.txt", "w") as f:
    f.write("\n".join(output_lines))

# Optional: also print to console
print("Structure summary saved to structure_summary.txt")
print("\n".join(output_lines))
