from fastapi import FastAPI, HTTPException, UploadFile, File
from models import Molecules
from rdkit import Chem
import io
import csv

app = FastAPI()

molecules_db =[
    Molecules(identifier="m1", smile="CCO"),
    Molecules(identifier="m2", smile="CCN"),
    Molecules(identifier="m3", smile="CCCO"),
    Molecules(identifier="m4", smile="CC(=O)O"),
    Molecules(identifier="m5", smile="C1=CC=CC=C1")
]

@app.get("/")
def read_root():
    return {"Hello": "World"}


@app.get("/smiles")
def retrive_molecules():
    return molecules_db


@app.get("/smiles/{identifier}")
def retrieve_molecule(identifier: str):
    for molecule in molecules_db:
        if molecule.identifier == identifier:
            return molecule
    raise HTTPException(status_code=404, detail="Molecule not found")


@app.post("/add")
def add_molecules(molecule: Molecules):
    # I don't know if smiles all need to be uppercase so I will leave it as it is.
    for mol in molecules_db:
        if mol.identifier == molecule.identifier:
            raise HTTPException(status_code=400, detail="Molecule with this identifier already exists")
    molecules_db.append(molecule)
    return molecules_db


@app.put("/smiles/{identifier}")
def update_molecule(identifier:str, updated_smile: Molecules):
    for index, molecule in enumerate(molecules_db):
        if molecule.identifier == identifier:
            molecules_db[index] = updated_smile
            return updated_smile
    raise HTTPException(status_code=404, detail="Molecule not found")


@app.delete("/delete/{identifier}")
def delete_molecule(identifier: str):
    for index, molecule in enumerate(molecules_db):
        if molecule.identifier == identifier:
            molecules_db.pop(index)
            return {"detail": "Molecule deleted"}
    raise HTTPException(status_code=404, detail="Molecule not found")


@app.get("/substructures")
def substructure_search():
    # Checks if any of the Slimes are substructure of other Slimes. I hope this is what you meant.
    results = []

    for molecule in molecules_db:
        identifier = molecule.identifier
        smile = molecule.smile
        mol = Chem.MolFromSmiles(smile)
        if mol is None:
            continue

        substructures = []
        for substructure in molecules_db:
            if substructure.identifier == identifier:
                continue

            sub_mol = Chem.MolFromSmiles(substructure.smile)
            if sub_mol and mol.HasSubstructMatch(sub_mol):
                substructures.append({"identifier": substructure.identifier, "smile": substructure.smile})

        if substructures:
            results.append({
                "identifier": identifier,
                "substructures": substructures
            })

    return results


@app.post("/upload")
async def upload_file(file: UploadFile = File(...)):
    # Make sure to upload CSV with tab as delimiter
    if not file.filename.endswith('.csv'):
        raise HTTPException(status_code=400, detail="Only CSV files are supported")
    
    content = await file.read()
    # Decodes the file
    data = io.StringIO(content.decode("utf-8"))
    
    try:
        csv_reader = csv.DictReader(data, delimiter='\t')
    except Exception as e:
        raise HTTPException(status_code=400, detail=f"Failed to read CSV: {str(e)}")

    added_smiles = []
    
    for row in csv_reader:
        if 'identifier' not in row or 'smile' not in row:
            raise HTTPException(status_code=400, detail="File must contain 'identifier' and 'smile' columns")
        
        identifier = row['identifier']
        smile = row['smile']
        molecule = Molecules(identifier=identifier, smile=smile)
        
        if any(mol.identifier == molecule.identifier for mol in molecules_db):
            continue
        
        molecules_db.append(molecule)
        added_smiles.append(molecule)
    
    return {"added_molecules": added_smiles}