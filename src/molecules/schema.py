from pydantic import BaseModel, Field

class MoleculeCreate(BaseModel):
    id: int = Field(..., description="Molecule ID")
    name: str = Field(..., min_length=1, max_length=100, description="Molecule name")

class MoleculeRead(BaseModel):
    id: int
    name: str = Field(..., min_length=1, max_length=100, description="Molecule name")
