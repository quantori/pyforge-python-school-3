import os
from fastapi import FastAPI, HTTPException, Depends, UploadFile, File
from sqlalchemy.orm import Session
from models import MoleculeBase
from db import Molecule, get_db
from rdkit import Chem
import io
import csv

app = FastAPI()

@app.get("/")
def get_server():
    return {"server_id": os.getenv("SERVER_ID", "1")}

@app.get("/smiles")
def retrieve_molecules(db: Session = Depends(get_db)):
    """
    Endpoint to retrieve all smiles from the database.
    Returns a dictionary with all smiles.
    """
    molecules = db.query(Molecule).all()
    return molecules

@app.get("/smiles/{identifier}")
def retrieve_molecule(identifier: str, db: Session = Depends(get_db)):
    """
    Endpoint to retrieve a specific smile by its identifier.
    Args:
        identifier: The unique identifier of the smile.
    Returns:
        The smile with the given identifier if it exists.
    Raises:
        HTTPException: If the smile with the given identifier does not exist.
    """
    molecule = db.query(Molecule).filter(Molecule.identifier == identifier).first()
    if molecule is None:
        raise HTTPException(status_code=404, detail="Molecule not found")
    return molecule

@app.post("/add")
def add_molecule(molecule: MoleculeBase, db: Session = Depends(get_db)):
    """
    Endpoint to add a new smile to the database.
    Args:
        molecule: The smile object to be added.
    Returns:
        The added smile.
    Raises:
        HTTPException: If a smile with the same identifier already exists.
    """
    db_molecule = db.query(Molecule).filter(Molecule.identifier == molecule.identifier).first()
    if db_molecule:
        raise HTTPException(status_code=400, detail="Molecule with this identifier already exists")
    new_molecule = Molecule(identifier=molecule.identifier, smile=molecule.smile)
    db.add(new_molecule)
    db.commit()
    db.refresh(new_molecule)
    return new_molecule

@app.put("/smiles/{identifier}")
def update_molecule(identifier: str, updated_molecule: MoleculeBase, db: Session = Depends(get_db)):
    """
    Endpoint to update an existing smile's information.
    Args:
        identifier: The unique identifier of the smile to be updated.
        updated_molecule: The updated smile object.
    Returns:
        The updated smile.
    Raises:
        HTTPException: If the smile to be updated does not exist or if the new identifier already exists.
    """
    molecule = db.query(Molecule).filter(Molecule.identifier == identifier).first()
    if molecule is None:
        raise HTTPException(status_code=404, detail="Molecule not found")

    if updated_molecule.identifier != identifier:
        db_molecule_with_new_id = db.query(Molecule).filter(Molecule.identifier == updated_molecule.identifier).first()
        if db_molecule_with_new_id:
            raise HTTPException(status_code=400, detail="New identifier already exists")

    molecule.identifier = updated_molecule.identifier
    molecule.smile = updated_molecule.smile
    db.commit()
    db.refresh(molecule)
    return molecule

@app.delete("/delete/{identifier}")
def delete_molecule(identifier: str, db: Session = Depends(get_db)):
    """
    Endpoint to delete a specific smile by its identifier.
    Args:
        identifier: The unique identifier of the smile to be deleted.
    Returns:
        A confirmation message upon successful deletion.
    Raises:
        HTTPException: If the smile with the given identifier does not exist.
    """
    molecule = db.query(Molecule).filter(Molecule.identifier == identifier).first()
    if molecule is None:
        raise HTTPException(status_code=404, detail="Molecule not found")
    db.delete(molecule)
    db.commit()
    return {"detail": "Molecule deleted"}

@app.get("/substructures")
def substructure_search(db: Session = Depends(get_db)):
    """
    Search for substructures within the molecules database.
    Returns a list of molecules that contain substructures from other molecules.
    """
    molecules = db.query(Molecule).all()
    molecule_objects = {mol.identifier: Chem.MolFromSmiles(mol.smile) for mol in molecules}

    results = []
    substructure_matches = {identifier: [] for identifier in molecule_objects}

    for identifier, mol in molecule_objects.items():
        if mol is None:
            continue

        for sub_id, sub_mol in molecule_objects.items():
            if identifier == sub_id or sub_mol is None:
                continue

            if mol.HasSubstructMatch(sub_mol):
                substructure_matches[identifier].append({
                    "identifier": sub_id,
                    "smile": next(m.smile for m in molecules if m.identifier == sub_id)
                })

    for identifier, matches in substructure_matches.items():
        if matches:
            results.append({
                "identifier": identifier,
                "substructures": matches
            })

    return results

@app.post("/upload")
async def upload_file(file: UploadFile = File(...), db: Session = Depends(get_db)):
    """
    Endpoint to upload a CSV file containing smile data.
    The CSV file must have 'identifier' and 'smile' columns and use a tab delimiter.
    Args:
        file: The uploaded CSV file.
    Returns:
        A list of added smile.
    Raises:
        HTTPException: If the file format is incorrect or if there are issues with the CSV content.
    """
    if not file.filename.endswith('.csv'):
        raise HTTPException(status_code=400, detail="Only CSV files are supported")

    content = await file.read()
    data = io.StringIO(content.decode("utf-8"))

    try:
        csv_reader = csv.DictReader(data, delimiter='\t')
    except Exception as e:
        raise HTTPException(status_code=400, detail=f"Failed to read CSV: {str(e)}")

    added_smiles = []

    for row in csv_reader:
        identifier, smile = row.get('identifier'), row.get('smile')
        if not identifier or not smile:
            raise HTTPException(status_code=400, detail="File must contain 'identifier' and 'smile' columns")

        if db.query(Molecule).filter(Molecule.identifier == identifier).first() is None:
            molecule = Molecule(identifier=identifier, smile=smile)
            db.add(molecule)
            db.commit()
            added_smiles.append(molecule)

    return {"added_molecules": added_smiles}
