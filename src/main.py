from fastapi import FastAPI, HTTPException
from pydantic import BaseModel, ValidationInfo, field_validator
from chem import substructure_search, valid_smile
import random

app = FastAPI()

db = []


class Molecule(BaseModel):
    smile: str

    @field_validator("smile")
    def check_valid_smile(cls, v: str, info: ValidationInfo):
        if not valid_smile(v):
            raise ValueError(f"{v} is not a valid SMILES string.")
        return v


def find_by_id(id: int):
    for m in db:
        if m["id"] == id:
            return m
    return None


def find_index_by_id(id: int):
    for i, m in enumerate(db):
        if m["id"] == id:
            return i
    return -1


def delete_by_id(id: int):
    molecule = find_by_id(id)
    if molecule:
        db.remove(molecule)
        return True
    else:
        return False


def get_all():
    return db


def get_filtered(substructre: str):
    return substructure_search(db, substructre)


def update_by_id(id: int, new_molecule: Molecule):
    index = find_index_by_id(id)

    if index == -1:
        return None

    molecule = new_molecule.model_dump()
    molecule["id"] = id
    db[index] = molecule
    return molecule


@app.get("/molecules", tags=["molecules"])
def get_all_molecules(substructre: str | None = None):
    if substructre:
        return get_filtered(substructre)
    return get_all()


@app.get("/molecules/{molecule_id}", tags=["molecules"])
def get_molecule_by_id(molecule_id: int):
    molecule = find_by_id(molecule_id)
    if not molecule:
        raise HTTPException(
            status_code=404, detail=f"Molecule not found by id: {molecule_id}"
        )
    return molecule


@app.post("/molecules", tags=["molecules"], status_code=201)
def create_molecule(request: Molecule):
    molecule = request.model_dump()
    molecule["id"] = random.randint(1, 10000)
    db.append(molecule)
    return molecule


@app.put("/molecules/{molecule_id}", tags=["molecules"])
def update_molecule_by_id(molecule_id: int, request: Molecule):
    updated_molecule = update_by_id(molecule_id, request)
    if not updated_molecule:
        raise HTTPException(
            status_code=404, detail=f"Molecule not found by id: {molecule_id}"
        )
    return updated_molecule


@app.delete("/molecules/{molecule_id}", tags=["molecules"], status_code=204)
def delete_molecule_by_id(molecule_id: int):
    deleted = delete_by_id(molecule_id)
    if not deleted:
        raise HTTPException(
            status_code=404, detail=f"Molecule not found by id: {molecule_id}"
        )
