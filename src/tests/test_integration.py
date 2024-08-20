"""
Let's test the integration of the whole system
"""

import pytest
from sqlalchemy.orm import sessionmaker
from sqlalchemy import create_engine, StaticPool
from src.database import Base
from src.dependencies import get_molecule_service
from src.repositories import MoleculeRepository
from src.models import Molecule
from src.service import MoleculeService
from src.tests.sample_data import (
    alkanes,
    to_molecule_request_dict,
    is_equal_dict_without_id,
)
from fastapi.testclient import TestClient
from src.main import app

# engine = create_engine("postgresql://user:password@localhost:5432/db_test")
engine = create_engine(
    "sqlite:///:memory:",
    connect_args={"check_same_thread": False},
    poolclass=StaticPool,
)
session_factory = sessionmaker(bind=engine)
molecule_repository = MoleculeRepository(session_factory)
molecule_service = MoleculeService(molecule_repository)

client = TestClient(app)
app.dependency_overrides[get_molecule_service] = lambda: molecule_service


@pytest.fixture
def init_db_3_alkanes():
    """
    Create the database schema and add first 3 alkanes to the database.
    """
    Base.metadata.drop_all(engine)
    Base.metadata.create_all(engine)
    # add data to the database, DO NOT FORGET TO LOOK AT BULK INSERT LATER
    with session_factory() as session:
        for molecule in list(alkanes.values())[0:3]:
            session.add(Molecule(**molecule))
            session.commit()
    yield
    Base.metadata.drop_all(engine)


def test_add_molecule_find_all(init_db_3_alkanes):
    decane_request = to_molecule_request_dict(alkanes["decane"])
    response = client.post("/molecules/", json=decane_request)
    response_json = response.json()

    assert response.status_code == 201
    assert is_equal_dict_without_id(response_json, alkanes["decane"])

    response = client.get("/molecules/")
    response_json = response.json()
    assert response.status_code == 200
    assert len(response_json) == 4


def test_add_find_all_pagination(init_db_3_alkanes):
    decane_request = to_molecule_request_dict(alkanes["decane"])
    response = client.post("/molecules/", json=decane_request)
    response_json = response.json()

    assert response.status_code == 201
    assert is_equal_dict_without_id(response_json, alkanes["decane"])

    response = client.get("/molecules/")
    response_json = response.json()
    assert response.status_code == 200
    assert len(response_json) == 4

    response = client.get("/molecules/?page=1&page_size=2")
    response_json = response.json()
    assert response.status_code == 200
    assert len(response_json) == 2
    assert is_equal_dict_without_id(response_json[0], alkanes["propane"])
    assert is_equal_dict_without_id(response_json[1], alkanes["decane"])


def test_add_get_molecule_by_id(init_db_3_alkanes):
    decane_request = to_molecule_request_dict(alkanes["decane"])
    response = client.post("/molecules/", json=decane_request)
    response_json = response.json()
    decane_id = response_json["molecule_id"]

    response = client.get(f"/molecules/{decane_id}")
    response_json = response.json()
    assert response.status_code == 200
    assert is_equal_dict_without_id(response_json, alkanes["decane"])


# @pytest.mark.xfail
def test_add_duplicate_smiles(init_db_3_alkanes):
    methane_request = to_molecule_request_dict(alkanes["methane"])
    response = client.post("/molecules/", json=methane_request)
    assert response.status_code == 400


# @pytest.mark.xfail
def test_add_molecule_invalid_smiles(init_db_3_alkanes):
    invalid_smiles = {"smiles": "i hate my life and i want to die"}
    response = client.post("/molecules/", json=invalid_smiles)
    assert response.status_code == 400


# @pytest.mark.xfail
def test_get_molecule_by_id_not_found(init_db_3_alkanes):
    response = client.get("/molecules/1000")
    assert response.status_code == 404


def test_update_molecule(init_db_3_alkanes):
    decane_request = to_molecule_request_dict(alkanes["decane"])
    response = client.post("/molecules/", json=decane_request)
    response_json = response.json()
    decane_id = response_json["molecule_id"]

    updated_decane = alkanes["decane"].copy()
    updated_decane["name"] = "Decane Updated"
    response = client.put(f"/molecules/{decane_id}", json=updated_decane)
    response_json = response.json()

    assert response.status_code == 200
    assert is_equal_dict_without_id(response_json, updated_decane)

    response = client.get(f"/molecules/{decane_id}")
    response_json = response.json()
    assert response.status_code == 200
    assert is_equal_dict_without_id(response_json, updated_decane)


# @pytest.mark.xfail
def test_update_molecule_not_found(init_db_3_alkanes):
    updated_decane = alkanes["decane"].copy()
    updated_decane["name"] = "Decane Updated"
    response = client.put("/molecules/1000", json=updated_decane)
    assert response.status_code == 404


def test_find_all_delete_by_id(init_db_3_alkanes):
    response = client.get("/molecules/")
    response_json = response.json()
    assert response.status_code == 200
    assert len(response_json) == 3

    for molecule in list(alkanes.values())[0:3]:
        response = client.delete(f"/molecules/{molecule['molecule_id']}")
        assert response.status_code == 200

    response = client.get("/molecules/")
    response_json = response.json()
    assert response.status_code == 200
    assert len(response_json) == 0


def test_substructure_search(init_db_3_alkanes):
    response = client.get("/substructure_search?smiles=CC")
    response_json = response.json()
    assert response.status_code == 200
    assert len(response_json) == 2
    assert is_equal_dict_without_id(response_json[0], alkanes["methane"])
    assert is_equal_dict_without_id(response_json[1], alkanes["ethane"])

