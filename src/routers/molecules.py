from fastapi import APIRouter, HTTPException

from src.db.molecule import (
    get_filtered,
    get_all,
    find_by_id,
    create,
    update_by_id,
    delete_by_id,
)
from src.models.molecule import Molecule

router = APIRouter()


@router.get("")
async def get_all_molecules(substructre: str | None = None):
    if substructre:
        return get_filtered(substructre)
    return get_all()


@router.get("/{molecule_id}")
async def get_molecule_by_id(molecule_id: int):
    molecule = find_by_id(molecule_id)
    if not molecule:
        raise HTTPException(
            status_code=404, detail=f"Molecule not found by id: {molecule_id}"
        )
    return molecule


@router.post("", status_code=201)
async def create_molecule(request: Molecule):
    return create(request)


@router.put("/{molecule_id}")
async def update_molecule_by_id(molecule_id: int, request: Molecule):
    updated_molecule = update_by_id(molecule_id, request)
    if not updated_molecule:
        raise HTTPException(
            status_code=404, detail=f"Molecule not found by id: {molecule_id}"
        )
    return updated_molecule


@router.delete("/{molecule_id}", status_code=204)
async def delete_molecule_by_id(molecule_id: int):
    deleted = delete_by_id(molecule_id)
    if not deleted:
        raise HTTPException(
            status_code=404, detail=f"Molecule not found by id: {molecule_id}"
        )
