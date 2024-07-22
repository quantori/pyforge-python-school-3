from fastapi.testclient import TestClient

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
