from fastapi import FastAPI, HTTPException, status, UploadFile, File, Depends
from sqlalchemy.ext.asyncio import AsyncSession
from sqlalchemy.future import select
from rdkit import Chem
import json
import logging
from os import getenv
import uvicorn
from database import get_db
from molecules.models import Molecule as MoleculeModel
from molecules.schema import MoleculeAdd, MoleculeResponse
from molecules.dao import MoleculeDAO
from molecules.request_body import RBMolecule
from typing import List, Dict
from dotenv import load_dotenv
from typing import Union
from database import async_session_maker


logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

app = FastAPI()  

load_dotenv(".env")

DB_URL = getenv("DB_URL")

if DB_URL is None:
    logger.error("Database URL is not in env")
    raise ValueError("Database URL is not in env")

@app.get("/")
def get_server():
    server_id = getenv("SERVER_ID", "1")
    logger.info(f"Server {server_id} received a request to the root endpoint")
    return {"server_id": server_id}, {"message": "Welcome to Molecules app!"}

@app.post("/molecules")
async def add_molecule(molecule: MoleculeAdd):
    try:
        logger.info(f"Adding molecule: {molecule.name}")
        existing_molecule = await MoleculeDAO.add_mol(**molecule.model_dump())
        if existing_molecule:
            logger.warning(f"Molecule with ID {molecule.id} already exists")
            raise HTTPException(status_code=400, detail="Molecule with this ID already exists")
        if not Chem.MolFromSmiles(molecule.name):
            logger.warning(f"Invalid SMILES: {molecule.name}")
            raise HTTPException(status_code=400, detail=f"Invalid SMILES: {molecule.name}")
    except Exception as e:
        logger.error(f"Error occurred while adding molecule: {str(e)}")
        raise HTTPException(status_code=500, detail="Internal Server Error")

@app.get(
    "/molecules",
    tags=["Molecules"],
    status_code=status.HTTP_200_OK,
    response_description="List of all molecules",
)
async def retrieve_molecules(request_body: RBMolecule = Depends(),) -> list[MoleculeResponse]:
    logger.info("Retrieving molecules")
    return await MoleculeDAO.find_all_molecules(**request_body.to_dict())

@app.get(
    "/molecules/{item_id}",
    tags=["Molecules"],
    status_code=status.HTTP_200_OK,
    response_description="Get molecule by ID",
)
async def get_mol_by_id(item_id: int) -> Union[MoleculeResponse, dict]:
    logger.info(f"Retrieving molecule with ID {item_id}")
    result = await MoleculeDAO.find_full_data(mol_id=item_id)
    if not result:
        logger.warning(f"Molecule with ID {item_id} not found")
        raise HTTPException(status_code=404, detail="Molecule not found")
    return result

@app.put(
    "/molecules/{mol_id}",
    status_code=status.HTTP_200_OK,
    tags=["Molecules"],
    response_description="Update molecule by ID",
)
async def update_molecule(mol_id: int, name: str):
    try:
        logger.info(f"Updating molecule with ID {mol_id} to new name {name}")
        await MoleculeDAO.update_mol(mol_id, name)
        return {"message": "Molecule updated successfully"}
    except ValueError as e:
        logger.warning(f"Update failed: {str(e)}")
        raise HTTPException(status_code=404, detail=str(e))
    except Exception as e:
        logger.error(f"Internal server error during molecule update: {str(e)}")
        raise HTTPException(status_code=500, detail="Internal Server Error")

@app.delete(
    "/molecules/{mol_id}",
    status_code=status.HTTP_200_OK,
    tags=["Molecules"],
    response_description="Delete molecule by ID"
)
async def delete_mol(mol_id: int) -> dict:
    logger.info(f"Deleting molecule with ID {mol_id}")
    molecule = await MoleculeDAO.delete_mol_by_id(mol_id=mol_id)
    if molecule:
        logger.info(f"Molecule with ID {mol_id} deleted")
        return {"message": f"The molecule with id {mol_id} is deleted!"}
    else:
        logger.warning(f"Molecule with ID {mol_id} not found for deletion")
        raise HTTPException(status_code=404, detail="Molecule not found")

@app.get(
    "/substructure_search",
    tags=["Molecules"],
    status_code=status.HTTP_200_OK,
    response_model=None  
)
async def substructure_search(substructure_name: str) -> List[Dict]:
    try:
        logger.info(f"Searching for molecules matching substructure: {substructure_name}")
        matches = await MoleculeDAO.find_by_substructure(substructure_name)
        if not matches:
            logger.warning(f"No molecules found matching the substructure: {substructure_name}")
            raise HTTPException(status_code=404, detail="No molecules found matching the substructure")
        return matches
    except ValueError as e:
        logger.warning(f"Bad request for substructure search: {str(e)}")
        raise HTTPException(status_code=400, detail=str(e))
    except Exception as e:
        logger.error(f"Internal server error during substructure search: {str(e)}")
        raise HTTPException(status_code=500, detail="Internal Server Error")

@app.post(
    "/upload_file/",
    status_code=status.HTTP_201_CREATED,
    tags=["File Upload"],
    response_description="File uploaded and molecules parsed successfully"
)
async def upload_file(file: UploadFile = File(...)):
    try:
        logger.info(f"Uploading and processing file: {file.filename}")
        content = await file.read()
        content = content.decode("utf-8")
        molecules = json.loads(content)
    except (UnicodeDecodeError, json.JSONDecodeError) as e:
        logger.warning(f"Failed to process file {file.filename}: {str(e)}")
        raise HTTPException(status_code=400, detail="Invalid file format")

    for molecule in molecules:
        if await MoleculeDAO.find_molecule_by_id(molecule["mol_id"]):
            logger.warning(f"Molecule with ID {molecule['mol_id']} already exists in the database")
            raise HTTPException(status_code=400, detail=f"Molecule with ID {molecule['mol_id']} already exists")
        
        if not Chem.MolFromSmiles(molecule.get("name")):
            logger.warning(f"Invalid SMILES in file: {molecule['name']}")
            raise HTTPException(status_code=400, detail=f"Invalid SMILES: {molecule['name']}")
        
        await MoleculeDAO.add_molecule(mol_id=molecule["mol_id"], name=molecule["name"])
        logger.info(f"Molecule with ID {molecule['mol_id']} added successfully")

    return {
        "message": "File uploaded and molecules parsed successfully",
        "num_molecules": len(molecules)
    }

if __name__ == "__main__":
    port = int(getenv("PORT", 8000))
    logger.info(f"Starting server on port {port}")
    uvicorn.run(app, host="0.0.0.0", port=port)
