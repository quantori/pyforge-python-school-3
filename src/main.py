from fastapi import FastAPI, Depends
from src.schemas import MoleculeRequest, MoleculeResponse
from src.dependencies import get_molecule_service
from src.service import MoleculeService


app = FastAPI()


@app.post("/molecules/", status_code=201)
def add_molecule(
    molecule_request: MoleculeRequest,
    service: MoleculeService = Depends(get_molecule_service),
) -> MoleculeResponse:
    return service.save(molecule_request)


@app.get("/molecules/{molecule_id}", status_code=200)
def get_molecule(
    molecule_id: int, service: MoleculeService = Depends(get_molecule_service)
) -> MoleculeResponse:
    return service.find_by_id(molecule_id)


@app.get("/molecules/", status_code=200)
def get_molecules(
    service: MoleculeService = Depends(get_molecule_service),
) -> list[MoleculeResponse]:
    return service.find_all()
