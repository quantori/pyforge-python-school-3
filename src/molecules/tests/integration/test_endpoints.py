import os
import random

import pytest
from fastapi.testclient import TestClient
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker

from src.config import get_settings
from src.database import Base
from src.main import app
from src.molecules.repository import MoleculeRepository
from src.molecules.service import get_molecule_service, MoleculeService
from src.molecules.tests.generate_csv_file import generate_testing_files
from src.molecules.tests.testing_utils import (
    alkane_request_jsons,
    validate_response_dict_for_ith_alkane,
    validate_response_dict_for_alkane,
)

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
def init_db():
    """
    Create the database schema and add first 100 alkanes to the database.
    """
    Base.metadata.drop_all(engine)
    Base.metadata.create_all(engine)
    yield
    Base.metadata.drop_all(engine)


def post_consecutive_alkanes(start_number, amount):
    """
    This is helper method to set up the rest of the tests. It posts alkanes from start_number to end_number
    and returns the responses.

    Every response is validated to ensure that the response is correct, so can be used in the rest of the tests
    directly.
    :param start_number: the number of the first alkane to post
    :param amount: the number of alkanes to post, auto-capped so that start_number + amount < 100,
    """
    amount = min(amount, 100 - start_number + 1)
    responses = []
    for i in range(start_number, start_number + amount):
        alkane_request = alkane_request_jsons[i]
        response = client.post("/molecules/", json=alkane_request)
        assert response.status_code == 201
        response_json = response.json()
        assert validate_response_dict_for_ith_alkane(response_json, i)
        responses.append(response_json)
    return responses


@pytest.mark.parametrize("i", [random.randint(1, 99) for _ in range(5)])
def test_save_molecule(i, init_db):
    post_consecutive_alkanes(i, 1)


@pytest.mark.parametrize("i", [random.randint(1, 99) for _ in range(5)])
def test_save_duplicate_smiles(i, init_db):
    res = post_consecutive_alkanes(i, 1)[0]
    duplicate_request = {"name": "Gaiozi", "smiles": res["smiles"]}
    response = client.post("/molecules/", json=duplicate_request)
    assert response.status_code == 400


@pytest.mark.parametrize("i", [random.randint(1, 99) for _ in range(5)])
def test_find_by_id(i, init_db):
    response = post_consecutive_alkanes(i, 1)[0]
    response_id = response["molecule_id"]
    response = client.get(f"/molecules/{response_id}")
    assert response.status_code == 200
    assert validate_response_dict_for_ith_alkane(response.json(), i)


@pytest.mark.parametrize("i", [random.randint(1, 99) for _ in range(5)])
def test_find_by_id_not_found(i, init_db):
    response = post_consecutive_alkanes(i, 1)[0]
    response_id = response["molecule_id"]
    response = client.get(f"/molecules/{response_id + 1}")
    assert response.status_code == 404


@pytest.mark.parametrize("i", [random.randint(1, 10) for _ in range(5)])
def test_update_molecule(i, init_db):
    """
    Update the molecule 2 times, and check that the molecule is updated correctly.
    """
    responses = post_consecutive_alkanes(1, 11)
    update_request1 = {"name": "UpdatedName"}
    update_request2 = {"name": "UpdatedName2"}

    for req in [update_request1, update_request2]:
        response = client.patch(f"/molecules/{responses[i]['molecule_id']}", json=req)
        assert response.status_code == 200
        assert response.json()["name"] == req["name"]
        js = response.json()
        js["name"] = responses[i]["name"]
        assert validate_response_dict_for_alkane(js, responses[i])


@pytest.mark.parametrize("i", [random.randint(1, 10) for _ in range(5)])
def test_update_molecule_not_found(i, init_db):
    post_consecutive_alkanes(1, 11)
    update_request = {"name": "UpdatedName"}
    response = client.patch(f"/molecules/{100000}", json=update_request)
    assert response.status_code == 404


@pytest.mark.parametrize("i", [random.randint(1, 10) for _ in range(5)])
def test_delete_molecule(i, init_db):
    responses = post_consecutive_alkanes(1, 11)
    response = client.delete(f"/molecules/{responses[i]['molecule_id']}")
    assert response.status_code == 200
    response = client.get(f"/molecules/{responses[i]['molecule_id']}")
    assert response.status_code == 404


@pytest.mark.parametrize("i", [random.randint(1, 19) for _ in range(10)])
def test_substructures(i, init_db):
    responses = post_consecutive_alkanes(1, 20)
    response = client.get(
        f"/molecules/search/substructures?smiles={responses[i]['smiles']}"
    )
    assert response.status_code == 200
    response_json = response.json()
    # every alkane up to i should be in the response
    for j in range(1, i + 1):
        assert validate_response_dict_for_ith_alkane(response_json[j], j + 1)


@pytest.mark.parametrize("i", [random.randint(1, 20) for _ in range(10)])
def test_superstructures(i, init_db):
    responses = post_consecutive_alkanes(1, 20)
    response = client.get(
        f"/molecules/search/superstructures?smiles={responses[i - 1]['smiles']}"
    )
    assert response.status_code == 200
    response_json = response.json()
    # every alkane from i should be in the response
    for j in range(i, 20):
        assert validate_response_dict_for_ith_alkane(response_json[j - i], j)


@pytest.fixture
def create_testing_files():
    generate_testing_files()
    yield

    os.remove("alkanes.csv")
    os.remove("invalid_header.csv")
    os.remove("decane_nonane_invalid_smiles.csv")


def test_file_upload(init_db, create_testing_files):
    """
    This tests alkanes.csv file which contains 10 alkanes, starting from methane to decane.

    Since the first 3 alkanes are already in the database, the number of molecules added should be 7.


    """
    post_consecutive_alkanes(1, 3)
    with open("alkanes.csv", "rb") as file:
        response = client.post(
            "/molecules/upload", files={"file": file}
        )
        response_json = response.json()
        assert response.status_code == 201
        assert response_json["number_of_molecules_added"] == 7


def test_file_upload_invalid_header(init_db, create_testing_files):
    """
    This tests invalid_header.csv file which has an invalid header.

    The response should be 400.
    """
    with open("invalid_header.csv", "rb") as file:
        response = client.post(
            "/molecules/upload", files={"file": file}
        )
        assert response.status_code == 400


def test_decane_nonane_invalid_smiles(init_db, create_testing_files):
    """
    This tests decane_nonane_invalid_smiles.csv file which has invalid smiles for decane and nonane.

    Only decane and nonane should not be added to the database, so the number of molecules added should be 5
    """
    post_consecutive_alkanes(1, 3)
    with open("decane_nonane_invalid_smiles.csv", "rb") as file:
        response = client.post(
            "/molecules/upload/", files={"file": file}
        )
        assert response.status_code == 201
        response_json = response.json()
        assert response_json["number_of_molecules_added"] == 5
