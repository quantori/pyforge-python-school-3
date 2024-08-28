from src.drugs.model import Drug
from src.drugs.schema import DrugResponse


def drug_to_response(drug: Drug) -> DrugResponse:
    return DrugResponse(
        drug_id=drug.drug_id,
        name=drug.name,
        description=drug.description,
    )
