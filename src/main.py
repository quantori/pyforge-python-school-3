
from rdkit import Chem
from rdkit.Chem import Draw

def substructure_search(mols, mol):
   
    
    substructure = Chem.MolFromSmiles(mol)
  
    matched_molecules = []

    for smiles in mols:
        molecule = Chem.MolFromSmiles(smiles)
        if molecule.HasSubstructMatch(substructure):
            matched_molecules.append(smiles)
    
    return matched_molecules


print(substructure_search(["CCO", "c1ccccc1", "CC(=O)O", "CC(=O)Oc1ccccc1C(=O)O"], "c1ccccc1"))