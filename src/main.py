from rdkit import Chem

def substructure_search(mols, mol):
    molecule = Chem.MolFromSmiles(mol)
    molecules_container = [Chem.MolFromSmiles(smiles) for smiles in mols]
    substructure_result = [molec for molec in molecules_container if molec.HasSubstructMatch(molecule)]
    return substructure_result
