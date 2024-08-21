from fastapi.testclient import TestClient
from unittest.mock import patch, AsyncMock
from src.molecules.router import router
from src.molecules.dao import MoleculeDAO
from src.molecules.schema import MoleculeAdd, MoleculeUpdate

client = TestClient(router)

# Test for adding a molecule
@patch.object(MoleculeDAO, 'add_molecule', new_callable=AsyncMock)
def test_add_molecule(mock_add_molecule):
    mock_add_molecule.return_value = 1
    response = client.post("/molecules", json={"smiles": "C1CCCCC1"})
    
    assert response.status_code == 201
    assert response.json() == {"message": "The molecule is added!", "molecule": {"smiles": "C1CCCCC1"}}
    mock_add_molecule.assert_called_once()

# Test for getting a molecule by ID
@patch.object(MoleculeDAO, 'find_full_data', new_callable=AsyncMock)
def test_get_molecule_by_id(mock_find_full_data):
    mock_find_full_data.return_value = {"id": 1, "smiles": "C1CCCCC1"}
    response = client.get("/molecules/1")
    
    assert response.status_code == 200
    assert response.json() == {"id": 1, "smiles": "C1CCCCC1"}
    mock_find_full_data.assert_called_once_with(molecule_id=1)

# Test for updating a molecule
@patch.object(MoleculeDAO, 'update', new_callable=AsyncMock)
def test_update_molecule(mock_update):
    mock_update.return_value = {"id": 1, "smiles": "C1CCCCC1"}
    response = client.put("/molecules/1", json={"smiles": "C1CCCCC1"})
    
    assert response.status_code == 200
    assert response.json() == {"id": 1, "smiles": "C1CCCCC1"}
    mock_update.assert_called_once_with(1, MoleculeUpdate(smiles="C1CCCCC1"))

# Test for deleting a molecule
@patch.object(MoleculeDAO, 'delete_molecule_by_id', new_callable=AsyncMock)
def test_delete_molecule(mock_delete_molecule_by_id):
    mock_delete_molecule_by_id.return_value = 1
    response = client.delete("/molecules/1")
    
    assert response.status_code == 200
    assert response.json() == {"message": "The molecule with id 1 is deleted!"}
    mock_delete_molecule_by_id.assert_called_once_with(molecule_id=1)

# Test for substructure search
@patch.object(MoleculeDAO, 'substructure_search', new_callable=AsyncMock)
def test_substructure_search(mock_substructure_search):
    mock_substructure_search.return_value = ["C1CCCCC1"]
    response = client.get("/substructures", params={"smiles": "C1CCCCC1"})
    
    assert response.status_code == 200
    assert response.json() == ["C1CCCCC1"]
    mock_substructure_search.assert_called_once_with("C1CCCCC1")
