import requests
import json

ENDPOINT = "http://localhost:8000"  # Update this to match your actual server URL

def upload_molecules_json(filename='molecules.json'):
    """Helper function to upload a molecules JSON file."""
    with open(filename, 'rb') as file:
        files = {'file': (filename, file, 'application/json')}
        response = requests.post(ENDPOINT + "/upload_file/", files=files)
    return response

def test_upload_file_success():
    response = upload_molecules_json('test_molecules.json')
    assert response.status_code == 201
    assert response.json() == {
        "message": "File uploaded and molecules parsed successfully",
        "num_molecules": 10
    }

def test_get_molecule_by_id():
    response = requests.get(ENDPOINT + "/molecules/4")
    assert response.status_code == 200
    assert response.json() == {"id": 4, "name": "CNC"}

def test_update_molecule():
    response = requests.put(
        ENDPOINT + "/molecules/2",
        json={"name": "CNO"}
    )
    assert response.status_code == 200
    assert response.json() == {"message": "Molecule updated successfully"}

def test_delete_molecule():
    response = requests.delete(ENDPOINT + "/molecules/5")
    assert response.status_code == 200
    assert response.json() == {"message": "The molecule with id 5 is deleted!"}

def test_substructure_search_valid_smiles():
    upload_molecules_json('test_molecules.json')
    response = requests.get(
        ENDPOINT + "/substructure_search/",
        params={"substructure_name": "C"}
    )
    assert response.status_code == 200
    data = response.json()
    assert "molecules" in data
    assert len(data["molecules"]) > 0

def test_upload_file_invalid_json():
    files = {
        'file': (
            'invalid_molecules.json',
            '{"mol_id": 5, "name": "C1=CC=CC=C1"',  # Missing closing bracket
            'application/json'
        )
    }
    response = requests.post(ENDPOINT + "/upload_file/", files=files)
    assert response.status_code == 400
    assert response.json() == {"detail": "Invalid file format"}
