

def generate_links_from_id(molecule_id: int):
    return {
        "self": {
            "href": f"/molecules/{molecule_id}",
            "rel": "self",
            "type": "GET",
        },
        "substructures": {
            "href": f"/substructure_search?molecule_id={molecule_id}",
            "rel": "substructures",
            "type": "GET",
        },
    }
