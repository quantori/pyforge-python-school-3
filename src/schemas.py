from typing import Annotated

from black.linegen import Optional
from pydantic import BaseModel, Field, field_validator

from src.exceptions import InvalidSmilesException
from src.utils import is_valid_smiles


class MoleculeRequest(BaseModel):
    smiles: Annotated[
        str,
        Field(
            ...,
            min_length=1,
            description="SMILES string of the molecule, should be unique",
        ),
    ]
    name: Annotated[Optional[str], Field(description="Name of the molecule")]

    @classmethod
    @field_validator("smiles")
    def validate_smiles(cls, smiles: str):
        if not is_valid_smiles(smiles):
            raise InvalidSmilesException(smiles)
        return smiles

    model_config = {
        "examples": {
            "aspirin": {
                "smiles": "CC(=O)Oc1ccccc1C(=O)O",
                "name": "Aspirin",
            },
            "Methane": {
                "smiles": "C",
            },
        }
    }


class MoleculeResponse(MoleculeRequest):
    molecule_id: Annotated[
        int, Field(..., description="Unique identifier for the molecule")
    ]

    model_config = {
        "examples": {
            "aspirin": {
                "molecule_id": 1,
                "smiles": "CC(=O)Oc1ccccc1C(=O)O",
                "name": "Aspirin",
            },
            "Methane": {
                "molecule_id": 2,
                "smiles": "C",
            },
        }
    }
