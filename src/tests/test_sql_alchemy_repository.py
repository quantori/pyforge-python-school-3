"""
Unit test Molecule Repository. Create sqlite for testing. Override the
get_session_maker function to return a session maker that uses the in-memory SQLite database.

Not fully Unit test, as it uses the database and does not mock, but that is sqlite in memory. I think it is fine, and
the best way to test the repository. Service will be mocked though.
"""

import pytest
from sqlalchemy.orm import sessionmaker
from sqlalchemy import create_engine
from src.database import Base
from src.repositories import MoleculeRepository
from src.models import Molecule
from src.tests.sample_data import alkanes, is_equal

engine = create_engine("sqlite:///:memory:", echo=True)
session_factory = sessionmaker(bind=engine)
molecule_repository = MoleculeRepository(session_factory)


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


def test_find_by_id(init_db_3_alkanes):
    molecule = molecule_repository.find_by_id(100)
    assert is_equal(molecule, alkanes["methane"])


def test_find_all(init_db_3_alkanes):
    """
    Warning: This test assumes that the order of the retrieved molecules is the same as the insertion order(Which it
    should be). Also, that the initial test data has less tan 1000 rows(which is the default page size currently).
    """

    molecules = molecule_repository.find_all()
    assert len(molecules) == 3
    assert is_equal(molecules[0], alkanes["methane"])
    assert is_equal(molecules[1], alkanes["ethane"])
    assert is_equal(molecules[2], alkanes["propane"])


@pytest.mark.parametrize(
    "page, page_size, expected_alkanes",
    [(0, 2, ("methane", "ethane")), (1, 2, ("propane",)), (1, 1, ("ethane",))],
)
def test_find_all_pagination(page, page_size, expected_alkanes, init_db_3_alkanes):
    molecules = molecule_repository.find_all(page=page, page_size=page_size)
    assert len(molecules) == len(expected_alkanes)
    for m in expected_alkanes:
        assert is_equal(molecules.pop(0), alkanes[m])


def test_save(init_db_3_alkanes):
    butane = alkanes["butane"]
    molecule = molecule_repository.save(butane)
    assert is_equal(molecule, butane)
    assert len(molecule_repository.find_all()) == 4
    assert is_equal(molecule_repository.find_by_id(103), butane)


def test_delete(init_db_3_alkanes):
    assert molecule_repository.delete(100)
    assert len(molecule_repository.find_all()) == 2
    assert molecule_repository.find_by_id(100) is None


def test_filter_by_name(init_db_3_alkanes):
    #     add another isomer of methane https://www.daylight.com/meetings/summerschool98/course/dave/smiles-isomers.html
    methane_isomer_dict = {"smiles": "[13CH4]", "name": "Methane", "molecule_id": 110}
    methane_isomer = Molecule(**methane_isomer_dict)

    with session_factory() as session:
        session.add(methane_isomer)
        session.commit()

    molecules = molecule_repository.filter(name="Methane")
    assert len(molecules) == 2
    assert is_equal(molecules[0], alkanes["methane"])
    assert is_equal(molecules[1], methane_isomer_dict)


def test_update(init_db_3_alkanes):
    methane_updated = alkanes["methane"].copy()
    methane_updated["name"] = "Methane2"
    molecule = molecule_repository.update(100, methane_updated)

    assert is_equal(molecule, methane_updated)
    assert len(molecule_repository.find_all()) == 3
    assert is_equal(molecule_repository.find_by_id(100), methane_updated)
