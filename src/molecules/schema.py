import datetime
from typing import Annotated
from black.linegen import Optional
from pydantic import BaseModel, Field, field_validator
from src.molecules.exception import InvalidSmilesException
from src.molecules.utils import is_valid_smiles
from src.schema import BaseResponse


class MoleculeRequest(BaseModel):
    smiles: Annotated[
        str,
        Field(
            min_length=1,
            description="SMILES string of the molecule, should be unique",
        ),
    ]
    name: Annotated[Optional[str], Field(description="Name of the molecule")]

    @field_validator("smiles")
    @classmethod
    def validate_smiles(cls, smiles: str):
        if not is_valid_smiles(smiles):
            raise InvalidSmilesException(smiles)
        return smiles

    model_config = {
        "json_schema_extra": {
            "examples": [
                {"smiles": "CC(=O)Oc1ccccc1C(=O)O", "name": "Aspirin"},
                {"smiles": "C", "name": "methane"},
            ]
        }
    }


class MoleculeResponse(MoleculeRequest, BaseResponse):
    molecule_id: Annotated[int, Field(description="Unique identifier for the molecule")]
    created_at: Annotated[
        datetime.datetime, Field(description="Timestamp when the molecule was created")
    ] = None
    updated_at: Annotated[
        datetime.datetime, Field(description="Timestamp when the molecule was updated")
    ] = None

    model_config = {
        "json_schema_extra": {
            "examples": [
                {
                    "molecule_id": 1,
                    "smiles": "CC(=O)Oc1ccccc1C(=O)O",
                    "name": "Aspirin",
                },
                {
                    "molecule_id": 2,
                    "smiles": "C",
                },
            ]
        }
    }


