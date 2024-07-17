from rdkit import Chem
from rdkit.Chem import Draw


def substructure_search(mols, mol):
    # Convert the substructure SMILES to an RDKit molecule
    substructure = Chem.MolFromSmiles(mol)

    # List to store molecules that contain the substructure
    matching_molecules = []
    matching_molecule_structures = []

    # Iterate through the list of molecule SMILES
    for smiles in mols:
        # Convert the SMILES to an RDKit molecule
        molecule = Chem.MolFromSmiles(smiles)

        # Check if the molecule contains the substructure
        if molecule.HasSubstructMatch(substructure):
            # If it matches, add the SMILES to the matching molecules list
            matching_molecules.append(smiles)
            matching_molecule_structures.append(molecule)

    # Visualize the matching molecules
    if matching_molecule_structures:
        img = Draw.MolsToGridImage(matching_molecule_structures, molsPerRow=2, subImgSize=(500, 500), returnPNG=False)
        img.show()

    return matching_molecules


"""
#Example:
substructure_search(["CCO", "c1ccccc1", "CC(=O)O", "CC(=O)Oc1ccccc1C(=O)O"], "c1ccccc1")
"""
