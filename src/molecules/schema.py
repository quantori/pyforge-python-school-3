from pydantic import BaseModel, Field

name_field = Field(..., min_length=1, max_length=100, description="Molecule name")

class MoleculeResponse(BaseModel):
    id: int
    name: str = name_field

class MoleculeAdd(BaseModel):
    id: int
    name: str = name_field
