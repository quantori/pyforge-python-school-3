from rdkit import Chem


def substructure_search(mols, mol):
    matching_molecules = []
    try:
        substructure_mol = Chem.MolFromSmiles(mol)

        for smiles in mols:
            molecule = Chem.MolFromSmiles(smiles)

            if molecule.HasSubstructMatch(substructure_mol):
                matching_molecules.append(smiles)
    except Exception as e:
        print(e)

    return matching_molecules
