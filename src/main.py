from fastapi import FastAPI

from src.molecules.molecule_exceptions import DuplicateSmilesException
from src.molecules.router import router as molecule_router
from src.drugs.router import router as drug_router
from src.handler import register_exception_handlers
from src.molecules.schemas import MoleculeRequest
from src.molecules.service import get_molecule_service

app = FastAPI()

# register the routers
app.include_router(molecule_router, prefix="/molecules")
app.include_router(drug_router, prefix="/drugs")
register_exception_handlers(app)
service = get_molecule_service()


@app.get("/")
def get_server_id():
    from os import getenv
    return "Hello from  server" + getenv("SERVER_ID", "1")


@app.on_event("startup")
def add_3_molecules():
    # add caffeine, surcose and water molecules
    try:
        service.save(
            MoleculeRequest.model_validate(
                {"smiles": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C", "name": "Caffeine"}
            )
        )
    except DuplicateSmilesException:
        pass

    try:
        service.save(
            MoleculeRequest.model_validate(
                {"smiles": "C(C1C(C(C(C(O1)O)O)O)O)O", "name": "Sucrose"}
            )
        )
    except DuplicateSmilesException:
        pass

    try:
        service.save(MoleculeRequest.model_validate({"smiles": "O", "name": "Water"}))
    except DuplicateSmilesException:
        pass
