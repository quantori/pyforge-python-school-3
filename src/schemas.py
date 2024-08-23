from pydantic import BaseModel


class Link(BaseModel):
    href: str
    rel: str
    type: str

    model_config = {
        "json_schema_extra": {
            "examples": [
                {"href": "/molecules/1", "rel": "self", "type": "GET"},
                {
                    "href": "/substructure_search?smiles=C",
                    "rel": "substructures",
                    "type": "GET",
                },
            ]
        }
    }


class BaseResponse(BaseModel):
    links: dict[str, Link] = {}

    model_config = {
        "json_schema_extra": {
            "examples": {
                "self": {"href": "/molecules/1", "rel": "self", "type": "GET"},
                "substructures": {
                    "href": "/substructure_search?smiles=C",
                    "rel": "substructures",
                    "type": "GET",
                },
            }
        }
    }
