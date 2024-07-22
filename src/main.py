from fastapi import FastAPI, HTTPException
from pydantic import BaseModel
from typing import List
import random

app = FastAPI()

db = []


class Molecule(BaseModel):
    smile: str


def find_by_id(id: int):
    for m in db:
        if m["id"] == id:
            return m
    return None


def delete_by_id(id: int):
    molecule = find_by_id(id)
    if molecule:
        db.remove(molecule)
        return True
    else:
        return False


def get_all():
    return db


def update_by_id(id: int, new_molecule: Molecule):
    pass


@app.get("/molecules", tags=["molecules"])
def get_all_molecules():
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
    pass


@app.delete("/molecules/{molecule_id}", tags=["molecules"], status_code=204)
def delete_molecule_by_id(molecule_id: int):
    deleted = delete_by_id(molecule_id)
    if not deleted:
        raise HTTPException(
            status_code=404, detail=f"Molecule not found by id: {molecule_id}"
        )
