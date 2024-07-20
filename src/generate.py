# this is to generate smiles of this format: id, molecule

from rdkit import Chem


def generate_molecules(n):
    mols = []
    for i in range(n):
        mol = f"{str(i)},aspirin,CC(=O)Oc1ccccc1C(=O)O"
        mols.append(mol)
    return mols

def create_file(n, filename):
    mols = generate_molecules(n)
    with open(filename, "w") as f:
        for mol in mols:
            f.write(f"{mol}\n")

create_file(10, "molecules.csv")
