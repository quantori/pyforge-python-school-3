from rdkit import Chem

def substructure_search(mols, mol):
    chem_mol = Chem.MolFromSmiles(mol)
    return [m for m in mols if Chem.MolFromSmiles(m).HasSubstructMatch(chem_mol)]