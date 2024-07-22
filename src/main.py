from fastapi import FastAPI, status, HTTPException, UploadFile
from rdkit import Chem
from src.models import Molecule, AddMoleculeRequest, MoleculeResponse
from src.repository import Repository, InMemoryMoleculesRepository
from src.utils import chem_from_smiles_http_error_if_invalid

app = FastAPI()
molecules_repository: Repository[int, Molecule] = InMemoryMoleculesRepository()


@app.post("/molecules",
          status_code=status.HTTP_201_CREATED,
          responses={status.HTTP_201_CREATED: {"description": "Molecule added successfully"},
                     status.HTTP_400_BAD_REQUEST: {"description": "Smiles string is invalid or missing"}}
          )
def add_molecule(add_molecule_request: AddMoleculeRequest) -> MoleculeResponse:
    """
    Add a new molecule to the repository.

    :param add_molecule_request: The request body.

    :return: The response containing the details of the added molecule.

    :raises HTTPException: If the provided SMILES string is invalid.
    """

    chem_from_smiles_http_error_if_invalid(add_molecule_request.smiles)
    molecule = Molecule.from_add_molecule_request(add_molecule_request)
    molecules_repository.add(molecule)
    return MoleculeResponse.from_molecule(molecule)


# @app.get("/molecules{molecule_id}",
#          status_code=status.HTTP_200_OK,
#          responses={status.HTTP_200_OK: {"description": "Molecule found successfully"},
#                     status.HTTP_404_NOT_FOUND: {"description": "Molecule with the provided ID is not found"}}
#          )
# def get_molecule(molecule_id: int) -> MoleculeResponse:
#     if not molecules_repository.exists_by_id(molecule_id):
#         raise HTTPException(status_code=status.HTTP_404_NOT_FOUND,
#                             detail=f"Molecule with the id {molecule_id} is not found")
#
#     molecule = molecules_repository.find_by_id(molecule_id)
#     return MoleculeResponse.from_molecule(molecule)

# @app.put("/molecules{molecule_id}", status_code=status.HTTP_200_OK)
# async def update_molecule(molecule_id: str, update_molecule_request: UpdateMoleculeRequest) -> Molecule:
#     if not molecules_repository.exists_by_id(molecule_id):
#         raise HTTPException(status_code=status.HTTP_404_NOT_FOUND,
#                             detail=f"Molecule with the id {molecule_id} is not found")
#
#     mol = molecules_repository.find_by_id(molecule_id)
#     mol.molecule_name = update_molecule_request.molecule_name
#     mol.smiles = update_molecule_request.smiles
#
#     return mol

#
# @app.delete("/molecules{molecule_id}", status_code=status.HTTP_200_OK)
# async def delete_molecule(molecule_id: str) -> None:
#     if not molecules_repository.exists_by_id(molecule_id):
#         raise HTTPException(status_code=status.HTTP_404_NOT_FOUND,
#                             detail=f"Molecule with the id {molecule_id} is not found")
#
#     molecules_repository.delete_by_id(molecule_id)
#
#     return
#
#
# @app.get("/molecules", status_code=status.HTTP_200_OK)
# async def list_molecules(skip: int = 0, limit: int = -1) -> list[Molecule]:
#     find_all = molecules_repository.find_all()
#
#     if limit == -1:
#         return find_all[skip:]
#
#     return find_all[skip:limit]
#
#
# @app.get("/molecules/substructure_search")
# async def get_substructure_search(smiles: str) -> list[Molecule]:
#     molecule = Chem.MolFromSmiles(smiles)
#     if molecule is None:
#         raise HTTPException(status_code=status.HTTP_400_BAD_REQUEST, detail=f"Invalid smiles format: {smiles}")
#
#     find_all = molecules_repository.find_all()
#
#     return [mol for mol in find_all if molecule.HasSubstructMatch(Chem.MolFromSmiles(mol.smiles))]


# @app.post("/upload_molecules", status_code=status.HTTP_201_CREATED)
# async def upload_molecules(file: UploadFile, background_tasks: BackgroundTasks):
#     contents = await file.read()
#     lines = contents.decode().splitlines()
#
#     background_tasks.add_task(process_molecules, lines)
#     return {"message": "Molecules are being processed in the background"}
#
#
# def process_molecules(lines: list[str]):
#     for line in lines:
#         parts = line.split(',')
#         if len(parts) != 3:
#             raise HTTPException(status_code=status.HTTP_400_BAD_REQUEST, detail="Invalid file format")
#
#         molecule_id, molecule_name, smiles = parts
#
#         if not is_valid_smiles(smiles):
#             raise HTTPException(status_code=status.HTTP_400_BAD_REQUEST, detail=f"Invalid SMILES string: {smiles}")
#
#         molecule = Molecule(molecule_id=molecule_id.strip(), molecule_name=molecule_name.strip(), smiles=smiles.strip())
#
#         if molecules_repository.exists_by_id(molecule.molecule_id):
#             raise HTTPException(status_code=status.HTTP_409_CONFLICT,
#                                 detail=f"Molecule with ID {molecule.molecule_id} already exists")
#
#         molecules_repository.add(molecule)
#
#
# def is_valid_smiles(smiles: str) -> bool:
#     mol = Chem.MolFromSmiles(smiles)
#     return mol is not None
#
