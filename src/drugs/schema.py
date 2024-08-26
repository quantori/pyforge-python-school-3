import enum
from typing import Annotated, Optional

from pydantic import BaseModel, Field

from src.schemas import BaseResponse


class MolecularConcentrationUnit(enum.Enum):
    MOLAR = "MOLAR"
    MASS = "MASS"
    VOLUME = "VOLUME"


class DrugRequest(BaseModel):
    name: Annotated[str, Field(..., min_length=1)]
    description: Annotated[Optional[str], Field()]

    model_config = {
        "json_schema_extra": {
            "examples": [
                {
                    "name": "Aspirin",
                    "description": "An analgesic medication",
                    "molecule_ids": [69],
                }
            ]
        }
    }


class DrugResponse(DrugRequest, BaseResponse):
    id: Annotated[
        int,
        Field(
            ...,
        ),
    ]


class DrugMoleculeRequest:
    molecule_id: Annotated[
        int,
        Field(
            ..., description="Reference to the molecule, should exist in the database"
        ),
    ]
    molecule_concentration: Annotated[
        float, Field(..., description="Concentration of the molecule in the drug")
    ]


class DrugMoleculeRequests:
    molecule_concentration_unit: Annotated[
        MolecularConcentrationUnit, Field(..., description="Unit of the concentration")
    ]
    molecules: Annotated[
        list[DrugMoleculeRequest],
        Field(..., description="List of molecules in the drug"),
    ]





