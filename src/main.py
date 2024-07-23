from typing import List
from fastapi import FastAPI, HTTPException, status, UploadFile, File
from pydantic import BaseModel
from rdkit import Chem
import csv


tags_metadata = [
    {
        "name": "Root",
        "description": "Root endpoint for the application.",
    },
    {
        "name": "Molecule Management",
        "description": "Operations related to managing molecules, including adding, updating, retrieving, and deleting molecules.",
    },
    {
        "name": "Substructure Search",
        "description": "Endpoints for performing substructure searches on molecules.",
    },
    {
        "name": "Upload a file",
        "description": "Endpoint for uploading a CSV file containing molecule data. The CSV should have a 'SMILES' column with molecule SMILES strings. This endpoint adds the molecules from the file to the existing database."
    }
]

app = FastAPI(openapi_tags=tags_metadata)


class Molecule(BaseModel):
    molecule_id: int
    smiles: str


class SubstructureQuery(BaseModel):
    substructure: str


molecules = [
    {"molecule_id": 1, "smiles": "CCO"},
    {"molecule_id": 2, "smiles": "c1ccccc1"},
    {"molecule_id": 3, "smiles": "CC(=O)O"},
    {"molecule_id": 4, "smiles": "CC(=O)Oc1ccccc1C(=O)O"}
]


@app.get("/", tags=["Root"])
def root():
    """Root endpoint for the application"""
    return {"message": "Molecule Fast API application"}


@app.post("/add", status_code=status.HTTP_201_CREATED, response_model=Molecule, tags=["Molecule Management"])
def add_molecule(molecule: Molecule):
    """Add a new molecule to the database"""

    for mol in molecules:
        if mol["molecule_id"] == molecule.molecule_id:
            raise HTTPException(status_code=status.HTTP_400_BAD_REQUEST, detail="Molecule with the ID: {molecule_id} already exists")
    molecules.append(molecule.model_dump())
    return molecule


@app.get("/molecule/{molecule_id}", response_model=Molecule, tags=["Molecule Management"])
def get_molecule(molecule_id: int):
    """Retrieve a molecule by its ID"""

    for mol in molecules:
        if mol["molecule_id"] == molecule_id:
            return mol
    raise HTTPException(status_code=status.HTTP_404_NOT_FOUND, detail="Molecule not found")


@app.put("/molecule/{molecule_id}", response_model=Molecule, tags=["Molecule Management"])
def update_molecule(molecule_id: int, molecule: Molecule):
    """Update a molecule by its ID"""

    for index, mol in enumerate(molecules):
        if mol["molecule_id"] == molecule_id:
            molecules[index] = molecule.model_dump()
            print(molecule.model_dump())
            return molecule.model_dump()
    raise HTTPException(status_code=status.HTTP_404_NOT_FOUND, detail="Molecule not found")


@app.delete("/molecule/{molecule_id}", status_code=status.HTTP_204_NO_CONTENT, tags=["Molecule Management"])
def delete_molecule(molecule_id: int):
    """Delete a molecule by its ID"""

    for index, mol in enumerate(molecules):
        if mol["molecule_id"] == molecule_id:
            molecules.pop(index)
            return
    raise HTTPException(status_code=status.HTTP_404_NOT_FOUND, detail="Molecule not found")


@app.get("/molecules", tags=["Molecule Management"])
def list_molecules():
    """List all molecules currently stored"""
    return molecules



@app.post("/search", status_code=status.HTTP_200_OK, response_model=List[Molecule], tags=["Substructure Search"])
def substructure_search(query: SubstructureQuery):
    """Perform substructure search using a query molecule"""

    substructure_molecule = Chem.MolFromSmiles(query.substructure)
    if not substructure_molecule:
        raise HTTPException(status_code=400, detail="Invalid substructure SMILES")
    
    matching_molecules = []
    for mol in molecules:
        molecule = Chem.MolFromSmiles(mol["smiles"])
        if molecule and molecule.HasSubstructMatch(substructure_molecule):
            matching_molecules.append(mol)

    return matching_molecules


@app.post("/uploadfile", status_code=status.HTTP_200_OK, tags=["Upload a file"])
async def upload_file(file: UploadFile = File(...)):
    """Upload a CSV file containing molecules (SMILES)"""
    
    if file.content_type != "text/csv":
        raise HTTPException(status_code=status.HTTP_400_BAD_REQUEST, detail="Invalid file format. Only CSV files are supported.")
    
    if file.filename == "":
        raise HTTPException(status_code=status.HTTP_400_BAD_REQUEST, detail="No file selected for upload.")
    
    contents = await file.read()
    if not contents:
        raise HTTPException(status_code=status.HTTP_400_BAD_REQUEST, detail="Uploaded file is empty.")
    
    try:
        reader = csv.DictReader(contents.decode().splitlines())
        if "SMILES" not in reader.fieldnames and "smiles" not in reader.fieldnames:
            raise HTTPException(status_code=status.HTTP_400_BAD_REQUEST, detail="CSV file must contain 'SMILES' or 'smiles' header.")
        
        new_molecules = []
        for index, row in enumerate(reader, start=len(molecules) + 1):
            smiles_value = row.get("SMILES") or row.get("smiles")
            if not smiles_value:
                raise HTTPException(status_code=status.HTTP_400_BAD_REQUEST, detail=f"Missing 'SMILES' or 'smiles' value in row {index}.")
            
            molecule = {
                "molecule_id": index,
                "smiles": smiles_value
            }
            new_molecules.append(molecule)
        
        molecules.extend(new_molecules)
        return {"message": f"Successfully uploaded {len(new_molecules)} molecules"}
    
    except csv.Error as e:
        raise HTTPException(status_code=status.HTTP_400_BAD_REQUEST, detail=f"Error processing CSV file: {str(e)}")
    except Exception as e:
        raise HTTPException(status_code=status.HTTP_500_INTERNAL_SERVER_ERROR, detail=f"Unexpected error: {str(e)}")
