from src.celery_worker import celery
from src.database import SessionLocal
from src.crud import get_all_molecules
from rdkit import Chem
import logging

FORMAT = '%(asctime)s | %(name)s | %(levelname)s | %(message)s'
logging.basicConfig(
    level=logging.INFO,
    format=FORMAT,
    handlers=[
        logging.FileHandler('app.log'),
        logging.StreamHandler(),
    ]
)

@celery.task
def async_substructure_search_task(substructure: str):
    logging.info(f"Starting async substructure search task for substructure: {substructure}")
    session = SessionLocal()

    try:
        desired_substructure = Chem.MolFromSmiles(substructure)
        if not desired_substructure:
            logging.warning(f'Substructure {substructure} is not valid.')
            raise ValueError("Invalid substructure SMILES")

        molecules = get_all_molecules(session)
        result = []

        for molecule in molecules:
            try:
                smile = molecule.smiles
                rdkit_molecule = Chem.MolFromSmiles(smile)

                if rdkit_molecule and rdkit_molecule.HasSubstructMatch(desired_substructure):
                    result.append({
                        "id": molecule.id,
                        "smiles": smile,
                        "name": molecule.name,
                        "weight": molecule.weight,
                        "formula": molecule.formula
                    })
            except Exception as e:
                logging.error(f'Error processing molecule with ID {molecule.id}: {e}')

        logging.info(f'Finished async substructure search for {substructure} with {len(result)} molecules.')
        return result

    except Exception as e:
        logging.error(f'Error in async substructure search for {substructure}: {e}')
        raise
    finally:
        session.close()
