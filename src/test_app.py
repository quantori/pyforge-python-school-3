from fastapi.testclient import TestClient
import pytest
from main import app, molecules_db


@pytest.fixture
def client():
    return TestClient(app)


def test_add_molecule(client):
    molecule = {
        '4': {'smiles': '[Na+].[Cl-]', 'molecule_formula': 'NaCl', 'molecule_weight': 58},
    }
    response = client.post('/add', json=molecule)
    assert response.status_code == 201
    assert response.json() == molecule
    assert '4' in molecules_db
    assert molecule['4'] == molecules_db['4']


expected_response_molecule = {'smiles': 'CCO', 'molecule_formula': 'C2H5OH', 'molecule_weight': 25}
expected_response_error = {'detail': 'Molecule is not found.'}

@pytest.mark.parametrize(
        'molecule_id, status_code, expected_response',
        [('1', 200, expected_response_molecule), ('10', 404, expected_response_error)],
) 
def test_retrieve_molecule(client, molecule_id, status_code, expected_response):    
    response = client.get(f'/molecules/{molecule_id}')
    assert response.status_code == status_code
    assert response.json() == expected_response


@pytest.mark.parametrize(
        'molecule_id, status_code, expected_response',
        [('1', 200, None), ('10', 404, expected_response_error)],
) 
def test_update_molecule(client, molecule_id, status_code, expected_response):
    molecule = {'smiles': 'CCO', 'molecule_formula': 'C2H5OH', 'molecule_weight': 27}
    response = client.put(f'/molecules/{molecule_id}', json=molecule)
    assert response.status_code == status_code
    assert response.json() == expected_response


@pytest.mark.parametrize('molecule_id, status_code', [('1', 200), ('10', 404)]) 
def test_delete_molecule(client, molecule_id, status_code):    
    response = client.delete(f'/molecules/{molecule_id}')
    assert response.status_code == status_code


def test_retrieve_all_molecules(client):
    response = client.get('/molecules/')
    assert response.status_code == 200
    assert response.json() == molecules_db

expected_molecules_response = {    
    '2': {'smiles': 'c1ccccc1', 'molecule_formula': 'C6H6', 'molecule_weight': 78},
    '3': {'smiles': 'CC(=O)Oc1ccccc1C(=O)O', 'molecule_formula': 'C9H8O4', 'molecule_weight': 180},
}

@pytest.mark.parametrize(
        'status_code, substructure_smiles, expected_response',
        [(200, 'c1ccccc1', expected_molecules_response), (200, 'CBr', {})],
) 
def test_substructure_search_molecules(client, status_code, substructure_smiles, expected_response):    
    response = client.get(f'/substructure_search/{substructure_smiles}')
    assert response.status_code == status_code
    assert response.json() == expected_response


