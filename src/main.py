from typing import List
from rdkit import Chem
from rdkit.Chem import Draw


def substructure_search(mols: List[str], mol: str) -> List[str]:
    """
    Find molecules containing a specified substructure.

    Args:
    - mols (list): List of molecules as SMILES strings.
    - mol (str): Substructure to search for as a SMILES string.

    Returns:
    - list: List of SMILES strings of molecules containing the substructure.
    """
    substructure = Chem.MolFromSmiles(mol)

    matches = [smile for smile in mols if Chem.MolFromSmiles(smile).HasSubstructMatch(substructure)]

    if not matches:
        print(f"No molecules found containing the substructure '{mol}'.")

    return matches


# Find molecules that contain the substructure "c1ccccc1" (Cyclohexane)
result = substructure_search(["CCO", "c1ccccc1", "CC(=O)O", "CC(=O)Oc1ccccc1C(=O)O"], "c1ccccc1")
print("Molecules containing the substructure 'c1ccccc1':", result)


# Display an image of molecules found
if result:
    img = Draw.MolsToGridImage([Chem.MolFromSmiles(smile) for smile in result])
    img.show()