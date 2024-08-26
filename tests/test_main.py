import pytest
from fastapi.testclient import TestClient
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
from src.main import app
from database.database import Base
from database.models import Molecule

# Test database URL
DATABASE_URL = "postgresql://postgres:password@localhost:5432/test_database"

# Create an engine and session local
test_engine = create_engine(DATABASE_URL)

# Initialize the database
Molecule.metadata.create_all(bind=test_engine)
TestingSessionLocal = sessionmaker(autocommit=False, autoflush=False, bind=test_engine)

# Create a TestClient instance
client = TestClient(app)


@pytest.fixture(scope="module")
def setup_teardown():

    # Populate the database with initial data
    with TestingSessionLocal() as session:

        molecules = [
            Molecule(identifier="water", smiles="O"),
            Molecule(identifier="methane", smiles="C"),
            Molecule(identifier="ethanol", smiles="CCO"),
            Molecule(identifier="benzene", smiles="c1ccccc1"),
            Molecule(identifier="acetic_acid", smiles="CC(O)=O"),
            Molecule(identifier="aspirin", smiles="CC(=O)Oc1ccccc1C(=O)O")
        ]
        session.add_all(molecules)
        session.commit()

    yield

    # Teardown: Drop the database tables
    Base.metadata.drop_all(bind=test_engine)


def test_get_server():
    response = client.get("/")
    assert response.status_code == 200
    assert "server_id" in response.json()


def test_add_molecule(setup_teardown):
    response = client.post(
        "/molecule",
        json={"identifier": "new_test", "smiles": "CCO"}
    )
    assert response.status_code == 200
    assert response.json() == {"identifier": "new_test", "smiles": "CCO"}

    # Test adding a molecule with the same identifier (should fail)
    response = client.post(
        "/molecule",
        json={"identifier": "new_test", "smiles": "C"}
    )
    assert response.status_code == 400

    # Cleanup
    client.delete("/molecule/new_test")


def test_get_molecule(setup_teardown):
    response = client.get("/molecule/water")
    assert response.status_code == 200
    assert response.json() == {"identifier": "water", "smiles": "O"}

    # Test getting a non-existing molecule
    response = client.get("/molecule/non_existing")
    assert response.status_code == 400


def test_update_molecule(setup_teardown):
    response = client.put(
        "/molecule/acetic_acid",
        json={"identifier": "acetic_acid", "smiles": "CC(O)=OO"}
    )
    assert response.status_code == 200
    assert response.json() == {"identifier": "acetic_acid", "smiles": "CC(O)=OO"}

    # Test updating a non-existing molecule
    response = client.put(
        "/molecule/non_existing",
        json={"identifier": "non_existing", "smiles": "C"}
    )
    assert response.status_code == 400


def test_delete_molecule(setup_teardown):
    response = client.delete("/molecule/ethanol")
    assert response.status_code == 200
    assert response.json() == {"message": "Molecule deleted successfully."}

    # Check the deletion
    response = client.get("/molecule/ethanol")
    assert response.status_code == 400


def test_list_molecules(setup_teardown):
    response = client.get("/molecules/")
    assert response.status_code == 200
    result = response.json()
    assert {"identifier": "ethanol", "smiles": "CCO"} in result


def test_search_substructure(setup_teardown):
    response = client.post("/search/", json={"substructure": "CCO"})
    assert response.status_code == 200
    result = response.json()
    assert len(result) > 0
    assert {"identifier": "ethanol", "smiles": "CCO"} in result
    assert {"identifier": "aspirin", "smiles": "CC(=O)Oc1ccccc1C(=O)O"} in result

    # Test invalid substructure
    response = client.post("/search/", json={"substructure": "invalid"})
    assert response.status_code == 400
