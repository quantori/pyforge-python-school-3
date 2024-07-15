from rdkit import Chem


def substructure_search(structures_smiles: list[str], substructure_smiles: str) -> list[str]: 
    search_result = []
    substructure_mol = Chem.MolFromSmiles(substructure_smiles)

    for structure_smiles in structures_smiles:
        structure_mol = Chem.MolFromSmiles(structure_smiles)

        if structure_mol.HasSubstructMatch(substructure_mol):
            search_result.append(structure_smiles)

    return search_result
