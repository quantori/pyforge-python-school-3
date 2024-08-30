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
