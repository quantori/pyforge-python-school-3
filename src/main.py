from fastapi import FastAPI, status, HTTPException, UploadFile
from rdkit import Chem
from src.repository import MoleculesRepository, InMemoryMoleculesRepository
from src.models import Molecule, UpdateMoleculeRequest

app = FastAPI()

molecules_repository: MoleculesRepository = InMemoryMoleculesRepository()


@app.post("/molecules", status_code=status.HTTP_201_CREATED)
def add_molecule(molecule: Molecule) -> Molecule:
    if molecules_repository.exists_by_id(molecule.molecule_id):
        raise HTTPException(status_code=status.HTTP_400_BAD_REQUEST,
                            detail=f"molecule with the same molecule_id {molecule.molecule_id} already exists")

    chem_molecule = Chem.MolFromSmiles(molecule.smiles)

    if chem_molecule is None:
        raise HTTPException(status_code=status.HTTP_400_BAD_REQUEST, detail=f"Invalid smiles format: {molecule.smiles}")

    molecules_repository.add(molecule)
    return molecule


@app.get("/molecules{molecule_id}", status_code=status.HTTP_200_OK)
def get_molecule(molecule_id: str) -> Molecule:
    if not molecules_repository.exists_by_id(molecule_id):
        raise HTTPException(status_code=status.HTTP_404_NOT_FOUND,
                            detail=f"Molecule with the id {molecule_id} is not found")

    return molecules_repository.find_by_id(molecule_id)


@app.put("/molecules{molecule_id}", status_code=status.HTTP_200_OK)
def update_molecule(molecule_id: str, update_molecule_request: UpdateMoleculeRequest) -> Molecule:
    if not molecules_repository.exists_by_id(molecule_id):
        raise HTTPException(status_code=status.HTTP_404_NOT_FOUND,
                            detail=f"Molecule with the id {molecule_id} is not found")

    mol = molecules_repository.find_by_id(molecule_id)
    mol.molecule_name = update_molecule_request.molecule_name
    mol.smiles = update_molecule_request.smiles

    return mol


@app.delete("/molecules{molecule_id}", status_code=status.HTTP_200_OK)
def delete_molecule(molecule_id: str) -> None:
    if not molecules_repository.exists_by_id(molecule_id):
        raise HTTPException(status_code=status.HTTP_404_NOT_FOUND,
                            detail=f"Molecule with the id {molecule_id} is not found")

    molecules_repository.delete_by_id(molecule_id)

    return


@app.get("/molecules", status_code=status.HTTP_200_OK)
async def list_molecules(skip: int = 0, limit: int = -1) -> list[Molecule]:
    find_all = molecules_repository.find_all()

    if limit == -1:
        return find_all[skip:]

    return find_all[skip:limit]


@app.get("/molecules/substructure_search")
def get_substructure_search(smiles: str) -> list[Molecule]:
    molecule = Chem.MolFromSmiles(smiles)
    if molecule is None:
        raise HTTPException(status_code=status.HTTP_400_BAD_REQUEST, detail=f"Invalid smiles format: {smiles}")

    find_all = molecules_repository.find_all()

    return [mol for mol in find_all if molecule.HasSubstructMatch(Chem.MolFromSmiles(mol.smiles))]


@app.post("/upload_molecules", status_code=status.HTTP_201_CREATED)
def upload_molecules(file: UploadFile):
    contents = file.file.read()
    lines = contents.decode().splitlines()

    for line in lines:
        parts = line.split(',')
        if len(parts) != 3:
            raise HTTPException(status_code=status.HTTP_400_BAD_REQUEST, detail="Invalid file format")

        molecule_id, molecule_name, smiles = parts

        if not is_valid_smiles(smiles):
            raise HTTPException(status_code=status.HTTP_400_BAD_REQUEST, detail=f"Invalid SMILES string: {smiles}")

        molecule = Molecule(molecule_id=molecule_id.strip(), molecule_name=molecule_name.strip(), smiles=smiles.strip())

        if molecules_repository.exists_by_id(molecule.molecule_id):
            raise HTTPException(status_code=status.HTTP_409_CONFLICT,
                                detail=f"Molecule with ID {molecule.molecule_id} already exists")

        molecules_repository.add(molecule)

    return {"message": "Molecules uploaded successfully"}


def is_valid_smiles(smiles: str) -> bool:
    mol = Chem.MolFromSmiles(smiles)
    return mol is not None


