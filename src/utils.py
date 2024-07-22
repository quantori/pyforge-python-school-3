from fastapi import HTTPException, status
from rdkit import Chem


def chem_from_smiles_http_error_if_invalid(smiles: str) -> Chem.Mol:
    chem_molecule = Chem.MolFromSmiles(smiles)

    if chem_molecule is None:
        raise HTTPException(status_code=status.HTTP_400_BAD_REQUEST, detail=f"Invalid smiles format: {smiles}")

    return chem_molecule


