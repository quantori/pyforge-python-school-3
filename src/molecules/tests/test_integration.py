import pytest
from sqlalchemy.orm import sessionmaker
from sqlalchemy import create_engine
from src.database import Base
from src.molecules.service import get_molecule_service, MoleculeService
from src.molecules.molecule_repository import MoleculeRepository
from src.molecules.models import Molecule
from src.molecules.tests.sample_data import (
    alkanes,
    to_molecule_request_dict,
    is_equal_dict_without_id,
)
from fastapi.testclient import TestClient
from src.main import app
from src.configs import get_settings
from src.molecules.tests.generate_csv_file import generate_testing_files

# engine = create_engine("postgresql://user:password@localhost:5432/db_test")
engine = create_engine(
    get_settings().TEST_DB_URL,
)
session_factory = sessionmaker(bind=engine)
molecule_repository = MoleculeRepository()
molecule_service = MoleculeService(molecule_repository, session_factory)

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


@pytest.fixture
def create_testing_files():
    generate_testing_files()
    yield
    import os

    os.remove("alkanes.csv")
    os.remove("invalid_header.csv")
    os.remove("decane_nonane_invalid_smiles.csv")


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


def test_substructures(init_db_3_alkanes):
    response = client.get("molecules/search/substructures?smiles=CC")
    response_json = response.json()
    assert response.status_code == 200
    assert len(response_json) == 2
    assert is_equal_dict_without_id(response_json[0], alkanes["methane"])
    assert is_equal_dict_without_id(response_json[1], alkanes["ethane"])


def test_substructures_of(init_db_3_alkanes):
    response = client.get("molecules/search/substructure_of?smiles=CC")
    response_json = response.json()
    assert response.status_code == 200
    assert len(response_json) == 2
    assert is_equal_dict_without_id(response_json[1], alkanes["propane"])
    assert is_equal_dict_without_id(response_json[0], alkanes["ethane"])


def test_file_upload(init_db_3_alkanes, create_testing_files):
    """
    This tests alkanes.csv file which contains 10 alkanes, starting from methane to decane.

    Since the first 3 alkanes are already in the database, the number of molecules added should be 7.


    """
    with open("alkanes.csv", "rb") as file:
        response = client.post(
            "/molecules/upload/upload_molecules_csv", files={"file": file}
        )
        response_json = response.json()
        assert response.status_code == 201
        assert response_json["number_of_molecules_added"] == 7


def test_file_upload_invalid_header(init_db_3_alkanes, create_testing_files):
    """
    This tests invalid_header.csv file which has an invalid header.

    The response should be 400.
    """
    with open("invalid_header.csv", "rb") as file:
        response = client.post(
            "/molecules/upload/upload_molecules_csv", files={"file": file}
        )
        assert response.status_code == 400


def test_decane_nonane_invalid_smiles(init_db_3_alkanes, create_testing_files):
    """
    This tests decane_nonane_invalid_smiles.csv file which has invalid smiles for decane and nonane.

    Only decane and nonane should not be added to the database, so the number of molecules added should be 5
    """
    with open("decane_nonane_invalid_smiles.csv", "rb") as file:
        response = client.post(
            "/molecules/upload/upload_molecules_csv", files={"file": file}
        )
        assert response.status_code == 201
        response_json = response.json()
        assert response_json["number_of_molecules_added"] == 5
