import datetime
from typing import Annotated

from black.linegen import Optional
from fastapi import Query
from pydantic import BaseModel, Field, field_validator

from src.molecules.molecule_exceptions import InvalidSmilesException
from src.molecules.utils import is_valid_smiles
from src.schemas import BaseResponse


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
                {"smiles": "C"},
            ]
        }
    }


class MoleculeResponse(MoleculeRequest, BaseResponse):
    molecule_id: Annotated[
        int, Field(description="Unique identifier for the molecule")
    ]
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


class PaginationQueryParams(BaseModel):
    """Query parameters for paginated responses. Page is 0-indexed."""

    page: Annotated[int, Query(0, description="Page number", ge=0)] = 0
    page_size: Annotated[
        int, Query(1000, description="Number of items per page", ge=1)
    ] = 1000

    model_config = {"json_schema_extra": {"examples": [{"page": 0, "limit": 1000}]}}
