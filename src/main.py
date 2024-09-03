from rdkit import Chem

def substructure_search(molecules, substructure):
    """
    SMILES (Simplified Molecular Input Line Entry System) is a textual representation
    of the structure of a molecule, convenient for storing and transmitting information.
    """
    # Define the molecule and the substructure to search for
    desired_substructure = Chem.MolFromSmiles(substructure)

    # Collect molecules containing desired_substructure
    result = []

    for smile in molecules:

        # Create a molecule from a SMILES string
        molecule = Chem.MolFromSmiles(smile)

        # Perform the substructure search and append to the result
        if molecule.HasSubstructMatch(desired_substructure):
            result.append(smile)

    return result

molecules = ["CCO", "c1ccccc1", "CC(=O)O", "CC(=O)Oc1ccccc1C(=O)O"]
substructure = "c1ccccc1"
print(substructure_search(molecules, substructure))