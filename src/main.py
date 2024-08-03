from os import getenv
from rdkit import Chem
from fastapi import FastAPI, HTTPException, UploadFile
from pydantic import BaseModel

app = FastAPI()

class Molecule (BaseModel):
    id: int
    smiles: str

database = []

@app.get("/")
def get_server():
    return {"server_id": getenv("SERVER_ID", "1")}

#Add molecule (smiles) and its identifier
@app.post("/molecules", status_code=201, tags=["Molecules"], summary="Add new molecules to the DB", response_description="Molecule added successfully")
async def add_molecule(molecule:Molecule):
    '''
    Create a molecule with the following details:

    **id**: Each item must have a unique identifier  
    **smiles**: Molecule in SMILES format  
    '''
    for data in database:
        if data['id']==molecule.id:
            raise HTTPException(status_code=409, detail="Molecule identifier already exists")
    substructure = Chem.MolFromSmiles(molecule.smiles)
    if not substructure:
        raise HTTPException(status_code=400, detail="Invalid SMILES molecule")
    database.append({'id': molecule.id, 'name':molecule.smiles})
    return database[-1]

#Get molecule by identifier
@app.get("/molecules/{molecule_id}",  tags=["Molecules"], summary="Retrieve molecule by ID", response_description="Molecule retrieved successfully")
async def get_molecule(molecule_id: int):
    '''
    Get molecule by identifier:

    **id**: Required
    '''
    for data in database:
        if data['id']==molecule_id:
            return data
    raise HTTPException(status_code=404, detail="Molecule not found")
        
#Updating a molecule by identifier
@app.put("/molecules/{molecule_id}",  tags=["Molecules"], summary="Update molecule by ID", response_description="Molecule updated successfully")
async def update_molecule(molecule_id: int, updated_molecule:Molecule):
    '''
    Updating a molecule by identifier

    **id**: Required
    '''
    substructure = Chem.MolFromSmiles(updated_molecule.smiles)
    if not substructure:
        raise HTTPException(status_code=400, detail="Invalid SMILES molecule")
    for data in database:
        if data['id']==molecule_id:
            data['name'] = updated_molecule.smiles
            return data
    raise HTTPException(status_code=404, detail="Molecule not found")

#Delete a molecule by identifier
@app.delete("/molecules/{molecule_id}",  tags=["Molecules"], summary="Delete molecule by ID", response_description="Molecule Deleted")
async def delete_molecule(molecule_id: int):
    '''
    Delete a molecule by identifier

    **id**: Required
    '''
    for index, data in enumerate(database):
        if data['id']==molecule_id:
            deleted_mol = database.pop(index)
            return deleted_mol
    raise HTTPException(status_code=404, detail="Molecule not found")

#List all molecules
@app.get("/molecules", tags=["Search"], summary="Retrieve all molecules", response_description="Molecules retrieved successfully")
async def get_all():
    '''
    List all molecules
    '''
    return database

#Substructure search for all added molecules
@app.get("/substructures", tags=["Search"], summary="Substructure search molecules", response_description="Substructure match molecules retrieved successfully")
async def substructure_search(molecule: str):
    '''
    Substructure search

    **substructure**: smiles molecule required
    '''
    substructure = Chem.MolFromSmiles(molecule)
    if not substructure:
        raise HTTPException(status_code=400, detail="Invalid SMILES molecule")
    molecules = [data['name'] for data in database]
    matches = [smiles for smiles in molecules if Chem.MolFromSmiles(smiles).HasSubstructMatch(substructure)]
    return matches

#[Optional] Upload file with molecules (the choice of format is yours).
@app.post("/upload-manager", status_code=201, tags=["Upload"], summary="Upload new molecules to the DB", response_description="Molecules uploaded successfully")
async def create_upload_file(file: UploadFile):
    contents = await file.read()
    return {"content": contents}
