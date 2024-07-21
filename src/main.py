from fastapi import FastAPI
from pydantic import BaseModel
from typing import List

app = FastAPI()

db = {}


class Molecule(BaseModel):
    smile: str


@app.get("/molecules", tags=["molecules"])
def get_all()->List[str]:
    pass


@app.get("/molecules/{molecule_id}", tags=["molecules"])
def get_by_id(molecule_id: int):
    pass


@app.post("/molecules", tags=["molecules"])
def create(request: Molecule):
    pass


@app.put("/molecules/{molecule_id}", tags=["molecules"])
def update_by_id(request: Molecule, molecule_id: int):
    pass


@app.delete("/molecules/{molecule_id}", tags=["molecules"])
def delete_by_id(molecule_id: int):
    pass
