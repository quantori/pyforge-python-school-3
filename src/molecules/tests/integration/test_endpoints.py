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
from src.molecules.tests.sample_data import (
    alkane_request_jsons,
    validate_response_dict_for_ith_alkane,
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

import logging


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

    Every response is validated to ensure that the response is correct, so can be used in the rest of the tests directly.
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


@pytest.mark.parametrize("i", [random.randint(1, 100) for _ in range(5)])
def test_save_molecule(i, init_db):
    post_consecutive_alkanes(i, 1)


@pytest.mark.parametrize("i", [random.randint(1, 99) for _ in range(5)])
def test_save_duplicate_smiles(i, init_db):
    res = post_consecutive_alkanes(i, 1)[0]
    duplicate_request = {"name": "Gaiozi", "smiles": res["smiles"]}
    response = client.post("/molecules/", json=duplicate_request)
    assert response.status_code == 400


@pytest.mark.parametrigze("i", [random.randint(1, 19) for _ in range(5)])
def test_find_by_id(i, init_db):
    responses = post_consecutive_alkanes(1, 20)
    response = client.get(f"/molecules/{responses[i]['molecule_id']}")
    assert response.status_code == 200
    assert validate_response_dict_for_ith_alkane(response.json(), i)
