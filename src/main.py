# from rdkit import Chem
# from rdkit.Chem import Draw

# def substructure_search(mols, mol):
    
#     substructure = Chem.MolFromSmiles(mol)  # Creating an RDKit molecule object for the substructure
#     matching_mols = []  # Creating an empty list to store molecules that match the substructure

#     # Iterate over the list of molecule SMILES strings (mols)
#     for smiles in mols:
#         molecule = Chem.MolFromSmiles(smiles)  # Creating a RDKit molecule object for the current molecule

#         if molecule.HasSubstructMatch(substructure):  # Checking if the molecule contains the substructure
#             matching_mols.append(smiles)  # If it does, adding the SMILES string to the list of matching molecules
    
    
#     return matching_mols # The list matching molecules will be returned



# if __name__ == "__main__":
    
#     # 2. Basic Molecular Representations
#     # Create a molecule from a SMILES string
#     smiles = "CCO"  # Ethanol
#     molecule = Chem.MolFromSmiles(smiles)

#     # Draw the molecule
#     Draw.MolToImage(molecule)

#     # 3. Substructure Search Example
#     # Define the molecule and the substructure to search for
#     benzene = Chem.MolFromSmiles("c1ccccc1")
#     ethanol = Chem.MolFromSmiles("CCO")

#     # Perform the substructure search
#     match = ethanol.HasSubstructMatch(benzene)
#     print("Benzene ring found in ethanol:", match)

#     # 4. Molecular Visualization
#     # Create a list of molecules
#     smiles_list = ["CCO", "c1ccccc1", "CC(=O)O", "CC(=O)Oc1ccccc1C(=O)O"]
#     molecules = [Chem.MolFromSmiles(smiles) for smiles in smiles_list]

#     # Draw the molecules in a grid
#     img = Draw.MolsToGridImage(molecules, molsPerRow=2, subImgSize=(200, 200), returnPNG=False)
#     img.show()
    
#     # Example usage of substructure search function
#     substructure_smiles = "c1ccccc1"
#     matching_molecules = substructure_search(smiles_list, substructure_smiles)
#     print(matching_molecules)

from fastapi import FastAPI, HTTPException
from pydantic import BaseModel
from typing import Dict, List
from rdkit import Chem
from rdkit.Chem import Draw

app = FastAPI()

molecules: Dict[str, str] = {}


class Molecule(BaseModel):
    identifier: str
    smiles: str


class MoleculeUpdate(BaseModel):
    smiles: str


# Defining API endpoints
@app.post("/molecules/", status_code=201)
def add_molecule(molecule: Molecule):
    if molecule.identifier in molecules:
        raise HTTPException(status_code=400, detail="Identifier already exists")
    molecules[molecule.identifier] = molecule.smiles
    return molecule


@app.get("/molecules/{identifier}", response_model=Molecule)
def get_molecule(identifier: str):
    if identifier not in molecules:
        raise HTTPException(status_code=404, detail="Molecule not found")
    return Molecule(identifier=identifier, smiles=molecules[identifier])


@app.put("/molecules/{identifier}", response_model=Molecule)
def update_molecule(identifier: str, molecule_update: MoleculeUpdate):
    if identifier not in molecules:
        raise HTTPException(status_code=404, detail="Molecule not found")
    molecules[identifier] = molecule_update.smiles
    return Molecule(identifier=identifier, smiles=molecule_update.smiles)


@app.delete("/molecules/{identifier}", status_code=204)
def delete_molecule(identifier: str):
    if identifier not in molecules:
        raise HTTPException(status_code=404, detail="Molecule not found")
    del molecules[identifier]
    return


@app.get("/molecules/", response_model=List[Molecule])
def list_molecules():
    return [Molecule(identifier=identifier, smiles=smiles) for identifier, smiles in molecules.items()]


def substructure_search(substructure_smiles: str) -> List[str]:
    substructure = Chem.MolFromSmiles(substructure_smiles)
    if not substructure:
        raise HTTPException(status_code=400, detail="Invalid substructure SMILES")
    
    matching_ids = []
    for identifier, smiles in molecules.items():
        molecule = Chem.MolFromSmiles(smiles)
        if molecule and molecule.HasSubstructMatch(substructure):
            matching_ids.append(identifier)
    
    return matching_ids


@app.get("/substructure_search/", response_model=List[Molecule])
def search_substructure(substructure: str):
    matching_ids = substructure_search(substructure)
    result = [Molecule(identifier=identifier, smiles=molecules[identifier]) for identifier in matching_ids]
    return result

@app.get("/")
def read_root():
    return {"message": "Hello World!"}

if __name__ == "__main__":
    import uvicorn
    uvicorn.run(app, host="0.0.0.0", port=8000)
