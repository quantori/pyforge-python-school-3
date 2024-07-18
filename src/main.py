from rdkit import Chem

def substructure_search(smiles_list, sub_smiles):
    sub_mol = Chem.MolFromSmiles(sub_smiles) #for conversion substructure SMILES to RDKIT molecujle
    matches = []        #empty list for storage of molecules
    for smiles in smiles_list:  #converts each smilel to rdkit mol and if it is created it checks if it is valid structure
        mol = Chem.MolFromSmiles(smiles)
        if mol and mol.HasSubstructMatch(sub_mol):      
            #adds molecules that are matched to list
            matches.append(smiles)
    return matches

molecules = ["CCO", "c1ccccc1", "CC(=O)O", "CC(=O)Oc1ccccc1C(=O)O"]
substructure = "c1ccccc1"
result = substructure_search(molecules, substructure)
print(result) 
