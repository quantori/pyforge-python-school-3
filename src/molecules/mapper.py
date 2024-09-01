from src.molecules.schema import MoleculeResponse
from src.schema import Link


def generate_links_from_id(molecule_id: int):
    return {
        "self": Link.model_validate(
            {
                "href": f"/molecules/{molecule_id}",
                "rel": "self",
                "type": "GET",
            }
        ),
        "molecules/search/substructures": Link.model_validate(
            {
                "href": f"/molecules/search/substructures?molecule_id={molecule_id}",
                "rel": "substructures",
                "type": "GET",
            }
        ),
        "molecules/search/is_substructure_of": Link.model_validate(
            {
                "href": f"/molecules/search/is_substructure_of?molecule_id={molecule_id}",
                "rel": "substructures",
                "type": "GET",
            }
        ),
    }


def generate_links_from_id_and_smiles(molecule_id: int, smiles: str):
    return {
        "self": Link.model_validate(
            {
                "href": f"/molecules/{molecule_id}",
                "rel": "self",
                "type": "GET",
            }
        ),
        "molecules/search/substructures": Link.model_validate(
            {
                "href": f"/molecules/search/substructures?smiles={smiles}",
                "rel": "substructures",
                "type": "GET",
            }
        ),
        "molecules/search/is_substructure_of": Link.model_validate(
            {
                "href": f"/molecules/search/superstructures?smiles={smiles}",
                "rel": "substructures",
                "type": "GET",
            }
        ),
    }


def model_to_response(molecule):
    return MoleculeResponse(
        molecule_id=molecule.molecule_id,
        smiles=molecule.smiles,
        name=molecule.name,
        mass=molecule.mass,
        created_at=molecule.created_at,
        updated_at=molecule.updated_at,
        links=generate_links_from_id_and_smiles(molecule.molecule_id, molecule.smiles)
    )
