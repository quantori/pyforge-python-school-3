from typing import List
from rdkit import Chem
import logging

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

def substructure_search(mols: List[str], mol: str) -> List[str]:
    try:
        substructure = Chem.MolFromSmiles(mol)
        if substructure is None:
            logger.error(f"Invalid SMILES substructure: {mol}")
            return []
    except Exception as e:
        logger.exception(f"Error creating substructure molecule from SMILES: {mol}")
        return []
    
    substructure_matches = []

    for smiles in mols:
        try:
            molecule = Chem.MolFromSmiles(smiles)
            if molecule is None:
                logger.error(f"Invalid SMILEs molecule: {smiles}")
                continue
            if molecule.HasSubstructMatch(substructure):
                substructure_matches.append(smiles)
        except Exception as e:
            logger.exception(f"Error processing molecule SMILES: {smiles}")

    return substructure_matches

