from typing import List
from fastapi import APIRouter, HTTPException
from rdkit import Chem
import logging

from src.utils import substructure_search
from src.models import Molecule


logger = logging.getLogger(__name__)

router = APIRouter()

molecules: dict[str, str] = {}

@router.post("/add_molecule")
def add_molecule(molecule: Molecule):
    if molecule.identifier in molecules:
        raise HTTPException(status_code=400, detail="Identifier already exists")
    try:
        mol = Chem.MolFromSmiles(molecule.smiles)
        if mol is None:
            raise ValueError("Invalid SMILES string")
    except Exception as e:
        logger.exception(f"Error adding molecule with SMILES: {molecule.smiles}")
        raise HTTPException(status_code=400, detail="Invalid SMILES string")
    molecules[molecule.identifier] = molecule.smiles
    return {"message": "Molecule added successfully"}


@router.get("/get_molecule/{identifier}")
def get_molecule(identifier: str):
    if identifier not in molecules:
        raise HTTPException(status_code=404, detail="Molecule not found")
    return {"identifier": identifier, "smiles": molecules[identifier]}

@router.put("/update_molecule/{identifier}")
def update_molecule(identifier: str, molecule: Molecule):
    if identifier not in molecules:
        raise HTTPException(status_code=404, detail="Molecule not found")
    try:
        mol = Chem.MolFromSmiles(molecule.smiles)
        if mol is None:
            raise ValueError("Invalid SMILES string")
    except Exception as e:
        logger.exception(f"Error updating molecule with SMILES: {molecule.smiles}")
        raise HTTPException(status_code=400, detail="Invalid SMILES string")
    molecules[identifier] = molecule.smiles
    return {"message": "Molecule updated successfully"}

@router.delete("/delete_molecule/{identifier}")
def delete_molecule(identifier: str):
    if identifier not in molecules:
        raise HTTPException(status_code=404, detail="Molecule not found")
    del molecules[identifier]
    return {"message": "Molecule deleted successfully"}

@router.get("/list_molecules")
def list_molecules():
    return molecules

@router.get("/substructure_search")
def substructure_search_api(smiles: str):
    result = substructure_search(list(molecules.values()), smiles)
    return {"matches": result}
