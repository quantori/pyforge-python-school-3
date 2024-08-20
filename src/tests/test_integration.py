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
    is_equal,
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


def test_add_molecule(init_db_3_alkanes):
    decane_request = to_molecule_request_dict(alkanes["decane"])
    response = client.post("/molecules/", json=decane_request)
    response_json = response.json()

    assert response.status_code == 201
    assert is_equal_dict_without_id(response_json, alkanes["decane"])




