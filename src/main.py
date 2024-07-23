from rdkit import Chem

def substructure_search(mols: list, mol: str) -> list:
    substructure = Chem.MolFromSmiles(mol)
    matches = []
    for smiles in mols:
        chem_mol = Chem.MolFromSmiles(smiles)
        match = chem_mol.HasSubstructMatch(substructure)
        if match:
            matches.append(smiles)
    return matches


expected_matches = ['c1ccccc1', 'CC(=O)Oc1ccccc1C(=O)O']
assert substructure_search(["CCO", "c1ccccc1", "CC(=O)O", "CC(=O)Oc1ccccc1C(=O)O"], "c1ccccc1") == expected_matches, f"Expected {expected_matches}, but got {substructure_search(["CCO", "c1ccccc1", "CC(=O)O", "CC(=O)Oc1ccccc1C(=O)O"], "c1ccccc1")}"
print("All tests passed!")



