from fastapi.testclient import TestClient

import src.main
from src.main import app

client = TestClient(app)


# 3 test cases for the add_molecule endpoint

def test_add_molecule_success():
    request_data = {
        "smiles": "CC(=O)Oc1ccccc1C(=O)O",
        "molecule_name": "Aspirin",
        "description": "Aspirin is in a group of medications called salicylates..."
    }
    response = client.post("/molecules", json=request_data)
    assert response.status_code == 201

    response_data = response.json()
    assert response_data["smiles"] == request_data["smiles"]
    assert response_data["molecule_name"] == request_data["molecule_name"]
    assert response_data["description"] == request_data["description"]
    assert "molecule_id" in response_data
    assert "links" in response_data


def test_add_molecule_invalid_smiles():
    request_data = {
        "smiles": "invalid_smiles_string",
        "molecule_name": "Invalid",
        "description": "This is an invalid SMILES string."
    }
    response = client.post("/molecules", json=request_data)
    assert response.status_code == 400


def test_add_molecule_missing_smiles():
    request_data = {
        "molecule_name": "No SMILES",
        "description": "This request is missing the SMILES string."
    }
    response = client.post("/molecules", json=request_data)
    assert response.status_code == 422


# 2 test cases for the get_molecule endpoint

def test_get_molecule_success():
    request_data = {
        "smiles": "CC(=O)Oc1ccccc1C(=O)O",
        "molecule_name": "Aspirin",
        "description": "Aspirin is in a group of medications called salicylates..."
    }
    add_response = client.post("/molecules", json=request_data)
    assert add_response.status_code == 201
    print(add_response.json())

    add_response_data = add_response.json()
    molecule_id = add_response_data["molecule_id"]

    get_response = client.get(add_response_data["links"]["self"]["href"])
    assert get_response.status_code == 200

    get_response_data = get_response.json()
    assert get_response_data["smiles"] == request_data["smiles"]
    assert get_response_data["molecule_name"] == request_data["molecule_name"]
    assert get_response_data["description"] == request_data["description"]
    assert get_response_data["molecule_id"] == molecule_id
    assert "links" in get_response_data


def test_get_molecule_not_found():
    response = client.get("/molecules/69")
    assert response.status_code == 404


# 2 test cases for the update_molecule endpoint


def test_update_molecule_success():
    request_data = {
        "smiles": "CC(=O)Oc1ccccc1C(=O)O",
        "molecule_name": "Aspirin",
        "description": "Aspirin is in a group of medications called salicylates..."
    }
    add_response = client.post("/molecules", json=request_data)
    assert add_response.status_code == 201

    add_response_data = add_response.json()
    molecule_id = add_response_data["molecule_id"]

    update_request_data = {
        "smiles": "C",
        "molecule_name": "Carbon",
        "description": "Base of all life on Earth"
    }
    update_response = client.put(add_response_data["links"]["self"]["href"], json=update_request_data)
    assert update_response.status_code == 200

    update_response_data = update_response.json()
    assert update_response_data["smiles"] == update_request_data["smiles"]
    assert update_response_data["molecule_name"] == update_request_data["molecule_name"]
    assert update_response_data["description"] == update_request_data["description"]
    assert update_response_data["molecule_id"] == molecule_id
    assert "links" in update_response_data


def test_update_molecule_not_found():
    update_request_data = {
        "smiles": "C",
        "molecule_name": "Carbon",
        "description": "Base of all life on Earth"
    }
    response = client.put("/molecules/69", json=update_request_data)
    assert response.status_code == 404


# 1 test cases for the delete_molecule endpoint

def test_delete_molecule():
    request_data = {
        "smiles": "CC(=O)Oc1ccccc1C(=O)O",
        "molecule_name": "Aspirin",
        "description": "Aspirin is in a group of medications called salicylates..."
    }
    add_response = client.post("/molecules", json=request_data)
    assert add_response.status_code == 201

    add_response_data = add_response.json()
    molecule_id = add_response_data["molecule_id"]

    delete_response = client.delete(add_response_data["links"]["self"]["href"])
    assert delete_response.status_code == 200

    get_response = client.get(add_response_data["links"]["self"]["href"])
    assert get_response.status_code == 404


# 2 test cases for the list_molecules endpoint

def test_list_molecules():
    src.main.molecules_repository.clear()

    request_data = {
        "smiles": "CC(=O)Oc1ccccc1C(=O)O",
        "molecule_name": "Aspirin",
        "description": "Aspirin is in a group of medications called salicylates..."
    }
    add_response = client.post("/molecules", json=request_data)
    assert add_response.status_code == 201

    request_data = {
        "smiles": "C",
        "molecule_name": "Carbon",
        "description": "Base of all life on Earth"
    }
    add_response = client.post("/molecules", json=request_data)
    assert add_response.status_code == 201

    get_response = client.get("/molecules")
    assert get_response.status_code == 200
    # TODO assert the response data


    # TODO test for different skip and limit query parameters
