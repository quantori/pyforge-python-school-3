from rdkit import Chem


def is_valid_smiles(smiles: str) -> bool:
    """Check if a SMILES string represents a valid molecule."""
    return Chem.MolFromSmiles(smiles) is not None
