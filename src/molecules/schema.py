import re
from pydantic import BaseModel, Field, field_validator


class MoleculeResponse(BaseModel):
    id: int
    name: str = Field(..., min_length=1, max_length=100, description="Molecule name")


class MoleculeAdd(BaseModel):
    id: int
    name: str = Field(..., min_length=1, max_length=100, description="Molecule name")

