from typing import Annotated, Optional
from pydantic import BaseModel, Field
from src.schemas import BaseResponse
from src.drugs.model import QuantityUnit


class DrugMoleculeRequest(BaseModel):
    """
    Single molecule in the drug.

    Molecule can be referenced by molecule_id

    """

    molecule_id: Annotated[
        Optional[int],
        Field(description="Reference to the molecule, should exist in the database"),
    ]
    quantity: Annotated[
        float, Field(..., description="Quantity of the molecule in the single dose")
    ]
    quantity_unit: Annotated[
        QuantityUnit,
        Field(
            ..., description="Unit of the quantity of the molecule in the single dose"
        ),
    ]


class DrugRequest(BaseModel):
    name: Annotated[str, Field(..., min_length=1)]
    description: Annotated[Optional[str], Field()]
    molecules: Annotated[list[DrugMoleculeRequest], Field(..., min_length=0)]

    model_config = {
        "json_schema_extra": {
            "examples": [
                {
                    "name": "Coffe",
                    "description": "Best drug ever",
                    "molecules": [
                        {
                            "smiles": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
                            "quantity": 0.2,
                            "quantity_unit": "GRAM",
                        },
                        {
                            "smiles": "OC[C@H]1O[C@@H](O[C@H]2[C@@H](O)[C@@H](CO)O[C@H](O)[C@H]2O)C(O)[C@@H](O)[C@H]1O",
                            "quantity": 5,
                            "quantity_unit": "GRAM",
                        },
                        {"smiles": "O", "quantity": 150, "quantity_unit": "ML"},
                    ],
                }
            ]
        }
    }


class DrugMoleculeResponse(BaseResponse):
    molecule_id: Annotated[int, Field()]
    quantity: Annotated[float, Field()]
    quantity_unit: Annotated[QuantityUnit, Field()]


class DrugResponse(BaseResponse):
    """ """

    drug_id: Annotated[
        int,
        Field(
            ...,
        ),
    ]
    name: Annotated[str, Field()]
    description: Annotated[Optional[str], Field()]
    molecules: Annotated[list[DrugMoleculeResponse], Field()]
