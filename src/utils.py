from fastapi import HTTPException, status
from rdkit import Chem

from src.exception import InvalidSmilesException


def chem_from_smiles_error_if_invalid(smiles: str) -> Chem.Mol:
    """
    Convert a SMILES string to a RDKit molecule object. If the SMILES string is invalid, raise an HTTPException
    with status code 400 and a message indicating the invalid format.

    :param smiles: The SMILES string.

    :return: The RDKit molecule object.

    :raises HTTPException: If the SMILES string is invalid.
    """

    chem_molecule = Chem.MolFromSmiles(smiles)

    if chem_molecule is None:
        raise InvalidSmilesException(smiles=smiles)

    return chem_molecule


def validate_smiles(smiles: str) -> str:
    """
    Validate the SMILES string. If the SMILES string is invalid, raise an InvalidSmilesException.

    :param smiles: The SMILES string.

    :return: The validated SMILES string.

    :raises InvalidSmilesException: If the SMILES string is invalid or empty.
    """
    if not smiles or Chem.MolFromSmiles(smiles) is None:
        raise InvalidSmilesException(smiles=smiles)
    return smiles