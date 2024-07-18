def substructure_search(mols, mol):
    target = Chem.MolFromSmiles(mol)
    result = []
    
    for elem in mols:
        molecule = Chem.MolFromSmiles(elem)
        if molecule.HasSubstructMatch(target):
            result.append(str(elem))
            Draw.MolToImage(molecule)

    return result