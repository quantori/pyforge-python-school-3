from pydantic import BaseModel


class Molecule(BaseModel):
    molecule_id: str
    molecule_name: str | None
    smiles: str


class UpdateMoleculeRequest(BaseModel):
    molecule_name: str | None
    smiles: str
