from os import getenv

from fastapi import FastAPI, HTTPException
from pydantic import BaseModel
from typing import List, Dict
from rdkit import Chem

app = FastAPI()


class Molecule(BaseModel):
    identifier: int
    smiles: str


molecules: Dict[int, str] = {}


@app.post("/add")
def add_molecule(molecule: Molecule):
    if molecule.identifier in molecules:
        raise HTTPException(status_code=400, detail="Molecule identifier already exists")
    molecules[molecule.identifier] = molecule.smiles
    return molecule


@app.get("/molecule/{identifier}")
def get_molecule(identifier: int):
    smiles = molecules.get(identifier)
    if not smiles:
        raise HTTPException(status_code=404, detail="Molecule not found")
    return {"identifier": identifier, "smiles": smiles}


@app.put("/molecule/{identifier}")
def update_molecule(identifier: int, updated_molecule: Molecule):
    if identifier not in molecules:
        raise HTTPException(status_code=404, detail="Molecule not found")
    if identifier != updated_molecule.identifier:
        raise HTTPException(status_code=400, detail="Identifier mismatch")
    molecules[identifier] = updated_molecule.smiles
    return updated_molecule


@app.delete("/molecule/{identifier}")
def delete_molecule(identifier: int):
    if identifier not in molecules:
        raise HTTPException(status_code=404, detail="Molecule not found")
    del molecules[identifier]
    return {"detail": "Molecule deleted"}


@app.get("/molecules", response_model=List[Molecule])
def list_molecules():
    return [{"identifier": identifier, "smiles": smiles} for identifier, smiles in molecules.items()]


@app.get("/substructure_search", response_model=Dict[str, List[str]])
def substructure_search():
    results = {}
    for smiles in molecules.values():
        if smiles in results:
            continue
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            match_list = []
            for other_smiles in molecules.values():
                other_mol = Chem.MolFromSmiles(other_smiles)
                if other_mol and other_mol.HasSubstructMatch(mol):
                    match_list.append(other_smiles)
            results[smiles] = match_list

    return results
@app.get("/")
def get_server():
    return {"server_id": getenv("SERVER_ID", "1")}