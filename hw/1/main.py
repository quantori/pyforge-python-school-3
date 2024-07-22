from rdkit import Chem


def substructure_search(molecules, substructure):
    substructure_mol = Chem.MolFromSmiles(substructure)
    
    matching_molecules = []

    for smiles in molecules:
        mols = Chem.MolFromSmiles(smiles)
        
        if mols.HasSubstructMatch(substructure_mol):
            matching_molecules.append(smiles)

    return matching_molecules


if __name__ == "__main__":
    molecules = ["CCO", "c1ccccc1", "CC(=O)O", "CC(=O)Oc1ccccc1C(=O)O"]
    substructure = "c1ccccc1"
    print(substructure_search(molecules, substructure))
