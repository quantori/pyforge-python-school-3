from rdkit import Chem
from rdkit.Chem import Draw


def substructure_search(mols, mol):
    molecule = Chem.MolFromSmiles(mol)
    return list(filter(lambda x: x.HasSubstructMatch(molecule), map(lambda x: Chem.MolFromSmiles(x), mols)))


mols = ["CCO", "c1ccccc1", "CC(=O)O", "CC(=O)Oc1ccccc1C(=O)O"]
mol = "c1ccccc1"
matches = substructure_search(mols, mol)
img = Draw.MolsToGridImage(matches)
img.show()
