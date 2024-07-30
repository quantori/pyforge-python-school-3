import csv
from io import StringIO
import os
from typing import Any
from fastapi import FastAPI, HTTPException, UploadFile
from rdkit import Chem


def substructure_search(structures_smiles: list[str], substructure_smiles: str) -> list[str]: 
    search_result = []
    substructure_mol = Chem.MolFromSmiles(substructure_smiles)

    for structure_smiles in structures_smiles:
        structure_mol = Chem.MolFromSmiles(structure_smiles)

        if structure_mol.HasSubstructMatch(substructure_mol):
            search_result.append(structure_smiles)

    return search_result

def load_molecules_db_from_file(file):
    buffer = StringIO(file.read().decode('utf-8'))  # read file as bytes and decode bytes into text stream
    reader = csv.DictReader(buffer)

    molecules_db = {}

    for row in reader:
        molecule = {
            'smiles': row['SMILES'],
            'molecule_formula': row['FORMULA'],
            'molecule_weight': int(row['WEIGHT']),
        }
        molecules_db[int(row['INDEX_ID'])] = molecule    

    return molecules_db


app = FastAPI()

molecules_db: dict[int, dict[str, Any]] = {
    1: {'smiles': 'CCO', 'molecule_formula': 'C2H5OH', 'molecule_weight': 25},
    2: {'smiles': 'c1ccccc1', 'molecule_formula': 'C6H6', 'molecule_weight': 78},
    3: {'smiles': 'CC(=O)Oc1ccccc1C(=O)O', 'molecule_formula': 'C9H8O4', 'molecule_weight': 180},
}

# Add molecule (smiles) and its identifier.
@app.post('/add', status_code=201)
def add_molecule(molecule: dict):
    molecules_db.update(molecule)
    return molecule


# Get molecule by identifier.
@app.get('/molecules/{molecule_id}')
def retrieve_molecule(molecule_id: int):
    if molecule_id in molecules_db:
        return molecules_db[molecule_id]
    else:
        raise HTTPException(status_code=404, detail='Molecule is not found.')


# Updating a molecule by identifier.
@app.put('/molecules/{molecule_id}')
def update_molecule(molecule_id: int, updated_molecule: dict[str, Any]):
    if molecule_id in molecules_db:
        molecules_db[molecule_id].update(updated_molecule)
    else:
        raise HTTPException(status_code=404, detail='Molecule is not found.')


# Delete a molecule by identifier.
@app.delete('/molecules/{molecule_id}')
def delete_molecule(molecule_id: int):
    if molecule_id in molecules_db:
        del molecules_db[molecule_id]
    else:
        raise HTTPException(status_code=404, detail='Molecule is not found.')


# List all molecules.
@app.get('/molecules/')
def retrieve_all_molecules():
    return molecules_db


# Substructure search for all added molecules.
@app.get('/substructure_search/{substructure_smiles}')
def substructure_search_molecules(substructure_smiles: str):
    smiles_from_db = []
    smiles_x_db_index = {}
    
    for index, molecule in molecules_db.items():
        smiles_from_db.append(molecule['smiles'])
        smiles_x_db_index[molecule['smiles']] = index
    
    smiles_search = substructure_search(smiles_from_db, substructure_smiles)
    
    smiles_indexes = [smiles_x_db_index[smiles] for smiles in smiles_search]
    
    result = {smiles_index: molecules_db[smiles_index] for smiles_index in smiles_indexes}

    return result


# [Optional] Upload file with molecules (the choice of format is yours).
@app.post('/create_db/')
async def upload_molecules_from_file(file: UploadFile):
    molecules_db.update(load_molecules_db_from_file(file.file))
    return molecules_db


# [Optional] Load balancer. Add method to check balancing.
@app.get('/check_balancing/')
def get_server():
    return {'server_id': os.getenv('SERVER_ID', 'UNKNOWN')}
