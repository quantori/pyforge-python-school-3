import collections

from src.models import Molecule

# simples organic compounds
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
