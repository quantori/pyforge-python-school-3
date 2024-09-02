import collections

from src.molecules.model import Molecule
import logging

logger = logging.getLogger(__name__)

# simplest organic compounds
alkanes = collections.OrderedDict(
    {
        "methane": {"smiles": "C", "name": "Methane", "molecule_id": 100},
        "ethane": {"smiles": "CC", "name": "Ethane", "molecule_id": 101},
        "propane": {"smiles": "CCC", "name": "Propane", "molecule_id": 102},
        "butane": {"smiles": "CCCC", "name": "Butane", "molecule_id": 103},
        "pentane": {"smiles": "CCCCC", "name": "Pentane", "molecule_id": 104},
        "hexane": {"smiles": "CCCCCC", "name": "Hexane", "molecule_id": 105},
        "heptane": {"smiles": "CCCCCCC", "name": "Heptane", "molecule_id": 106},
        "octane": {"smiles": "CCCCCCCC", "name": "Octane", "molecule_id": 107},
        "nonane": {"smiles": "CCCCCCCCC", "name": "Nonane", "molecule_id": 108},
        "decane": {"smiles": "CCCCCCCCCC", "name": "Decane", "molecule_id": 109},
    }
)


def is_equal(molecule: Molecule, molecule_dict: dict):
    """
    :param molecule: Molecule model
    :param molecule_dict: dictionary representing essential args
    :return: True if the molecule and molecule_dict are equal, False otherwise
    """

    return (
            molecule.molecule_id == molecule_dict["molecule_id"]
            and molecule.smiles == molecule_dict["smiles"]
            and molecule.name == molecule_dict["name"]
    )


def is_equal_dict_without_id(molecule_dict1: dict, molecule_dict2: dict):
    """
    :param molecule_dict1: dictionary representing essential args
    :param molecule_dict2: dictionary representing essential args
    :return: True if the molecule_dict1 and molecule_dict2 are equal, False otherwise
    """

    return (
            molecule_dict1["smiles"] == molecule_dict2["smiles"]
            and molecule_dict1["name"] == molecule_dict2["name"]
    )


def to_molecule_request_dict(molecule_dict: dict):
    """
    Just get rid of the molecule_id key

    :param molecule_dict: dictionary representing essential args
    :return: MoleculeRequest instance
    """

    copyo = molecule_dict.copy()
    del copyo["molecule_id"]
    return copyo


alkane_request_jsons = {
    i: (lambda: {"name": f"Alkane {i}", "smiles": "C" * i})() for i in range(1, 100)
}


def validate_response_dict_for_alkane(response_dict, alkane):

    if response_dict["name"] != alkane["name"]:
        logger.error(f"response_dict['name'] != alkane['name']")
        return False

    if response_dict["smiles"] != alkane["smiles"]:
        logger.error(f"response_dict['smiles'] != alkane['smiles']")
        return False

    if response_dict["molecule_id"] is None:
        logger.error(f"response_dict['molecule_id'] is None")
        return False

    if response_dict["created_at"] is None:
        logger.error(f"response_dict['created_at'] is None")
        return False

    if response_dict["updated_at"] is None:
        logger.error(f"response_dict['updated_at'] is None")
        return False

    expected_links = {
        "self": {
            "href": f"/molecules/{response_dict['molecule_id']}",
            "rel": "self",
            "type": "GET",
        },
        "substructures": {
            "href": f"/molecules/search/substructures?smiles={response_dict['smiles']}",
            "rel": "substructures",
            "type": "GET",
        },
        "superstructures": {
            "href": f"/molecules/search/superstructures?smiles={response_dict['smiles']}",
            "rel": "superstructures",
            "type": "GET",
        },
    }

    if response_dict["links"] != expected_links:
        # print first difference
        for key in response_dict["links"]:
            if response_dict["links"][key] != expected_links[key]:
                logger.error(f"response_dict['links'][{key}] != expected_links[{key}]")
                return False
        return False

    return True


def validate_response_dict_for_ith_alkane(response_dict, i):
    alkane = alkane_request_jsons[i]
    return validate_response_dict_for_alkane(response_dict, alkane)