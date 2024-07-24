from fastapi import FastAPI, HTTPException
from models import Molecules

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

