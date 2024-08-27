import enum
from typing import Annotated, Optional

from pydantic import BaseModel, Field

from src.schemas import BaseResponse
from src.drugs.model import QuantityUnit

class DrugRequest(BaseModel):
    name: Annotated[str, Field(..., min_length=1)]
    description: Annotated[Optional[str], Field()]

    model_config = {
        "json_schema_extra": {
            "examples": [
                {
                    "name": "Aspirin",
                    "description": "An analgesic medication",
                    "molecule_ids": [69, 420],
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
    quantity: Annotated[
        float, Field(..., description="Quantity of the molecule in the single dose")
    ]


class DrugMoleculeRequests:
    quantity_unit: Annotated[
        QuantityUnit, Field(..., description="Unit ")
    ]
    molecules: Annotated[
        list[DrugMoleculeRequest],
        Field(..., description="List of molecules in the drug"),
    ]


