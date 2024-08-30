from src.drugs.model import Drug
from src.drugs.schema import DrugResponse, DrugMoleculeResponse
from src.molecules import mapper as molecule_mapper
from src.schema import Link


def generate_links(drug_id: int):
    return {
        "self": Link.model_validate(
            {
                "href": f"/drugs/{drug_id}",
                "rel": "self",
                "type": "GET",
            }
        ),
        "molecules": Link.model_validate(
            {
                "href": f"/drugs/{drug_id}/molecules",
                "rel": "molecules",
                "type": "GET",
            }
        ),
    }


def drug_to_response(drug: Drug) -> DrugResponse:
    molecule_responses = []

    for molecule in drug.molecules:
        molecule_responses.append(
            DrugMoleculeResponse(
                molecule_id=molecule.molecule_id,
                quantity=molecule.quantity,
                quantity_unit=molecule.quantity_unit,
                links=molecule_mapper.generate_links_from_id(molecule.molecule_id),
            )
        )

    return DrugResponse(
        drug_id=drug.drug_id,
        name=drug.name,
        description=drug.description,
        molecules=molecule_responses,
        links=generate_links(drug.drug_id),
    )
