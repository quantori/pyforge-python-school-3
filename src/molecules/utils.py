from rdkit import Chem

from src.molecules.exception import InvalidSmilesException


def is_valid_smiles(smiles: str) -> bool:
    """
    Check if a SMILES string represents a valid molecule.

    :param smiles: SMILES string
    :return: True if the SMILES string not empty and represents a valid molecule , False otherwise
    """
    if not smiles:
        return False
    return Chem.MolFromSmiles(smiles) is not None


def get_chem_molecule_from_smiles_or_raise_exception(smiles: str):
    """
    Get a RDKit molecule object from a SMILES string or raise an InvalidSmilesException.

    :param smiles: SMILES string
    :return: RDKit molecule object
    :raises InvalidSmilesException: if the SMILES string
    """

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise InvalidSmilesException(smiles)
    return mol
