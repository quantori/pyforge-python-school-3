from rdkit import Chem
from rdkit.Chem import Draw

def substructure_search(mols, mol):
    molecule = Chem.MolFromSmiles(mol)
    return [m for m in mols if Chem.MolFromSmiles(m).HasSubstructMatch(molecule)]



# test
mols = ["CCO", "c1ccccc1", "CC(=O)O", "CC(=O)Oc1ccccc1C(=O)O"]
mol = "c1ccccc1"
matches = substructure_search(mols, mol)
print(matches)
img = Draw.MolsToGridImage([Chem.MolFromSmiles(x) for x in matches])
img.show()
