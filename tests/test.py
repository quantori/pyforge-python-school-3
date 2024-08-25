import requests
import json

ENDPOINT = "http://localhost:8011"


def test_upload_file_success():
    """Test successful file upload."""
    response = upload_molecules_json()
    assert response.status_code == 201
    assert response.json() == {
        "message": "File uploaded and molecules parsed successfully",
        "num_molecules": 10
    }


def upload_molecules_json(filename='src/molecules.json'):
    """Helper function to upload a molecules JSON file."""
    with open(filename, 'rb') as file:
        files = {'file': (filename, file, 'application/json')}
        response = requests.post(ENDPOINT + "/upload_file/", files=files)
    return response


def test_get_server():
    """Test to check if the server is up and running."""
    response = requests.get(ENDPOINT)
    assert response.status_code == 200
    assert "server_id" in response.json()


def test_get_molecule_by_id():
    response = requests.get(ENDPOINT + "/molecules/4")
    assert response.status_code == 200
    assert response.json() == {"id": 4, "name": "CNC"}


def test_update_molecule():
    response = requests.put(
        ENDPOINT + "/molecules/2",
        params={"name": "CNO"}  
    )
    assert response.status_code == 200
    assert response.json() == {"message": "Molecule updated successfully"}


def test_delete_molecule():
    response = requests.delete(ENDPOINT + "/molecules/7")
    assert response.status_code == 200
    assert response.json() == {"message": "The molecule with id 7 is deleted!"}
    

def test_substructure_search_valid_smiles():
    upload_molecules_json()
    response = requests.get(
        ENDPOINT + "/substructure_search/",
        params={"substructure_name": "C"}
    )
    
    assert response.status_code == 200
    data = response.json()
    assert "molecules" in data
    assert isinstance(data["molecules"], list)
    assert len(data["molecules"]) > 0


def test_substructure_search_invalid_smiles():
    """Test substructure search with invalid SMILES."""
    upload_molecules_json()
    response = requests.get(
        ENDPOINT + "/substructure_search/",
        params={"substructure_name": "INVALID"}
    )
    assert response.status_code == 400
    assert response.json() == {"detail": "Invalid SMILES string"}


def test_substructure_search_empty_smiles():
    """Test substructure search with empty SMILES."""
    upload_molecules_json()
    response = requests.get(
        ENDPOINT + "/substructure_search/",
        params={"substructure_name": ""}
    )
    assert response.status_code == 400
    assert response.json() == {"detail": "Invalid SMILES string"}


def test_substructure_search_special_characters():
    """Test substructure search with special characters."""
    upload_molecules_json()
    response = requests.get(
        ENDPOINT + "/substructure_search/",
        params={"substructure_name": "CCO$"}
    )
    assert response.status_code == 400
    assert response.json() == {"detail": "Invalid SMILES string"}


def test_upload_file_invalid_json():
    """Test uploading an invalid JSON file."""
    files = {
        'file': (
            'molecules.json',
            '{"mol_id": 5, "name": "C1=CC=CC=C1"',
            'application/json'
        )
    }
    response = requests.post(ENDPOINT + "/upload_file/", files=files)
    assert response.status_code == 400
    assert response.json() == {"detail": "Invalid JSON file"}
