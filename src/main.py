from fastapi import FastAPI, HTTPException, status, UploadFile, File, Depends
from sqlalchemy.ext.asyncio import AsyncSession
from sqlalchemy.future import select
from rdkit import Chem
import json
from os import getenv
import uvicorn
from database import get_db
from molecules.models import Molecule as MoleculeModel
from molecules.schema import MoleculeCreate, MoleculeRead
from molecules.dao import MoleculeDAO
from molecules.request_body import RBMolecule
from typing import List
from dotenv import load_dotenv

app = FastAPI()  

load_dotenv(".env")

DB_URL = ("DB_URL")

if DB_URL is None:
    raise ValueError("Database URL is not in env")

@app.get("/")
def get_server():
    return {"server_id": getenv("SERVER_ID", "1")}, {"message": "Welcome to Molecules app!"}

@app.post(
    "/molecules",
    status_code=status.HTTP_201_CREATED,
    tags=["Molecules"],
    response_description="Molecule added successfully",
    response_model=MoleculeRead
)
async def add_molecule(molecule: MoleculeCreate, db: AsyncSession = Depends(get_db)):
    existing_molecule = await MoleculeDAO.find_by_id(molecule.id, db)
    
    if existing_molecule:
        raise HTTPException(
            status_code=400,
            detail="Molecule with this ID already exists"
        )

    if not Chem.MolFromSmiles(molecule.name):
        raise HTTPException(
            status_code=400,
            detail=f"Invalid SMILES: {molecule.name}"
        )
    
    molecules_db = MoleculeModel(id=molecule.id, name=molecule.name)
    await MoleculeDAO.add(molecules_db, db)
    return molecules_db

@app.get(
    "/molecules",
    tags=["Molecules"],
    status_code=status.HTTP_200_OK,
    response_description="List of all molecules",
    response_model=List[MoleculeRead]
)
async def retrieve_molecules(db: AsyncSession = Depends(get_db)):
    molecules = await MoleculeDAO.find_all(db)
    return molecules

@app.get(
    "/molecules/{item_id}",
    tags=["Molecules"],
    status_code=status.HTTP_200_OK,
    response_description="Get molecule by ID",
    response_model=MoleculeRead
)
async def get_mol_by_id(item_id: int, db: AsyncSession = Depends(get_db)):
    molecule = await MoleculeDAO.find_by_id(item_id, db)
    if not molecule:
        raise HTTPException(status_code=404, detail="Molecule not found")
    return molecule

@app.put(
    "/molecules/{mol_id}",
    status_code=status.HTTP_200_OK,
    tags=["Molecules"],
    response_description="Update molecule by ID",
    response_model=MoleculeRead
)
async def update_mol(mol_id: int, updated_mol: MoleculeCreate, db: AsyncSession = Depends(get_db)):
    molecule = await MoleculeDAO.find_by_id(mol_id, db)
    if not molecule:
        raise HTTPException(status_code=404, detail="Molecule not found")

    if not Chem.MolFromSmiles(updated_mol.name):
        raise HTTPException(
            status_code=400,
            detail=f"Invalid SMILES: {updated_mol.name}"
        )

    molecule.name = updated_mol.name
    await MoleculeDAO.update(molecule, db)
    return molecule

@app.delete(
    "/molecules/{mol_id}",
    status_code=status.HTTP_200_OK,
    tags=["Molecules"],
    response_description="Delete molecule by ID",
    response_model=MoleculeRead
)
async def delete_mol(mol_id: int, db: AsyncSession = Depends(get_db)):
    molecule = await MoleculeDAO.find_by_id(mol_id, db)
    if not molecule:
        raise HTTPException(status_code=404, detail="Molecule not found")
    
    await MoleculeDAO.delete(molecule, db)
    return molecule

@app.get(
    "/substructure_search/",
    tags=["Molecules"],
    status_code=status.HTTP_200_OK,
    response_description="Search for molecules with a specific substructure"
)
async def substructure_search(substructure_name: str, db: AsyncSession = Depends(get_db)):
    if not substructure_name:
        raise HTTPException(
            status_code=400,
            detail="Invalid substructure SMILES"
        )
    substructure_mol = Chem.MolFromSmiles(substructure_name)
    if substructure_mol is None:
        raise HTTPException(
            status_code=400,
            detail="Invalid substructure SMILES"
        )
    
    molecules = await MoleculeDAO.find_all(db)
    
    matches = []
    for molecule in molecules:
        each_mol = Chem.MolFromSmiles(molecule.name)
        if each_mol is None:
            print(f"Invalid SMILES: {molecule.name}")
            continue
        if each_mol.HasSubstructMatch(substructure_mol):
            matches.append(molecule)
    return {"molecules": matches}

@app.post(
    "/upload_file/",
    status_code=status.HTTP_201_CREATED,
    tags=["File Upload"],
    response_description="File uploaded and molecules parsed successfully"
)
async def upload_file(file: UploadFile = File(...), db: AsyncSession = Depends(get_db)):
    try:
        content = await file.read()
        content = content.decode("utf-8")
        molecules = json.loads(content)
    except (UnicodeDecodeError, json.JSONDecodeError) as e:
        raise HTTPException(status_code=400, detail="Invalid file format")

    for molecule in molecules:
        if await MoleculeDAO.find_by_id(molecule["id"], db):
            raise HTTPException(status_code=400, detail=f"Molecule with ID {molecule['id']} already exists")
        
        if not Chem.MolFromSmiles(molecule.get("name")):
            raise HTTPException(status_code=400, detail=f"Invalid SMILES: {molecule['name']}")
        
        db_molecule = MoleculeModel(id=molecule["id"], name=molecule["name"])
        await MoleculeDAO.add(db_molecule, db)

    return {
        "message": "File uploaded and molecules parsed successfully",
        "num_molecules": len(molecules)
    }


if __name__ == "__main__":
    port = int(getenv("PORT", 8000))
    uvicorn.run(app, host="0.0.0.0", port=port)
