import rdkit
from rdkit import Chem
# from rdkit.Chem import Draw


def substructure_search(mols: list, mol: str) -> list:
    substructure = Chem.MolFromSmiles(mol)
    matches = []
    for smiles in mols:
        chem_mol = Chem.MolFromSmiles(smiles)
        match = chem_mol.HasSubstructMatch(substructure)
        if match:
            matches.append(smiles)
    print(matches)

substructure_search(["CCO", "c1ccccc1", "CC(=O)O", "CC(=O)Oc1ccccc1C(=O)O"], "c1ccccc1")
