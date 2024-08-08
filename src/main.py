from fastapi import FastAPI, HTTPException, status, UploadFile, File
from models import Molecule
from rdkit import Chem


def substructure_search(mols, mol):
    substructure_mol = Chem.MolFromSmiles(mol)
    match = []
    
    for i in mols:
        each_mol = Chem.MolFromSmiles(i)
        if each_mol.HasSubstructMatch(substructure_mol):
            match.append(i)
    
    return match

molecules = ["CCO", "c1ccccc1", "CC(=O)O", "CC(=O)Oc1ccccc1C(=O)O"]
mol = "c1ccccc1"
print(substructure_search(molecules, mol))


# additional problem I try to solve using rdkit
# here I count the descriptors of the given smiles molecules
''' 
from rdkit import Chem
from rdkit.Chem import Descriptors


def description(smiles):
    mol = Chem.MolFromSmiles(smiles)
    descr = {}
    if mol is None:
        return "invalid smiles"
    
    descr['MolecularWeight'] = Descriptors.MolWt(mol)
    descr['NumHDonors'] = Descriptors.NumHDonors(mol)
    descr['NumHAcceptors'] = Descriptors.NumHAcceptors(mol)
    descr['NumRotatableBonds'] = Descriptors.NumRotatableBonds(mol)
    descr['TPSA'] = Descriptors.TPSA(mol)
    
    return descr

print(compute_description('CCO'))
'''