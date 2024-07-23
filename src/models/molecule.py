from pydantic import BaseModel, ValidationInfo, field_validator

from src.utils.chem import valid_smile


class Molecule(BaseModel):
    smile: str

    @field_validator("smile")
    def check_valid_smile(cls, v: str, info: ValidationInfo):
        if not valid_smile(v):
            raise ValueError(f"{v} is not a valid SMILES string.")
        return v
