from celery_worker import celery
from rdkit import Chem
from molecules.dao import MoleculeDAO


@celery.task
def substructure_search_task(substructure_name: str, limit: int):
    substructure_mol = Chem.MolFromSmiles(substructure_name)
    if substructure_mol is None:
        raise ValueError("Invalid SMILES string")

    iterator = MoleculeDAO.find_by_substructure_iterator(
        substructure_name,
        limit
    )
    return [match for match in iterator]
