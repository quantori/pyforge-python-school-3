"""
Unit test Molecule Service. Create sqlite for testing. Override the
get_session_maker function to return a session maker that uses the in-memory SQLite database.

Using mocks to test the service without the database.
"""

import pytest
from src.molecules.service import MoleculeService
from src.exceptions import UnknownIdentifierException
from src.molecules.models import Molecule
from src.molecules.tests.sample_data import alkanes
from src.molecules.tests.sample_data import is_equal
from src.molecules.schemas import MoleculeRequest


# Fixture to create a mock repository
@pytest.fixture
def mock_repository(mocker):
    return mocker.Mock()


# Fixture to create the service with the mock repository
@pytest.fixture
def molecule_service(mock_repository):
    return MoleculeService(repository=mock_repository)


def test_exists_by_id_when_molecule_exists(molecule_service, mock_repository):
    mock_repository.find_by_id.return_value = Molecule(**alkanes["methane"])
    result = molecule_service.exists_by_id(1)
    assert result is True
    mock_repository.find_by_id.assert_called_once_with(1)


def test_exists_by_id_when_molecule_does_not_exist(molecule_service, mock_repository):
    mock_repository.find_by_id.return_value = None
    result = molecule_service.exists_by_id(1)
    assert result is False
    mock_repository.find_by_id.assert_called_once_with(1)


def test_find_by_id_when_molecule_exists(molecule_service, mock_repository):
    mock_repository.find_by_id.return_value = Molecule(**alkanes["methane"])
    molecule = molecule_service.find_by_id(1)
    print(molecule, "tabata")
    assert is_equal(molecule, alkanes["methane"])
    mock_repository.find_by_id.assert_called_with(1)


def test_find_by_id_when_molecule_does_not_exist(molecule_service, mock_repository):
    mock_repository.find_by_id.return_value = None
    with pytest.raises(UnknownIdentifierException):
        molecule_service.find_by_id(1)


def test_save_molecule(molecule_service, mock_repository, mocker):
    # implementation has to filter by smiles to check for duplicates, so we need to mock the filter method
    mock_repository.filter = mocker.Mock(return_value=[])
    ethane_request = alkanes["ethane"].copy()
    del ethane_request["molecule_id"]
    molecule_request = MoleculeRequest(**ethane_request)
    molecule_service.save(molecule_request)
    mock_repository.save.assert_called_once_with(ethane_request)


def test_update_molecule(molecule_service, mock_repository):
    ethane_request = alkanes["ethane"].copy()
    del ethane_request["molecule_id"]
    molecule_request = MoleculeRequest(**ethane_request)
    mock_repository.find_by_id.return_value = Molecule(**alkanes["ethane"])
    molecule_service.update(101, molecule_request)
    mock_repository.update.assert_called_once_with(101, ethane_request)


def test_delete_molecule(molecule_service, mock_repository):
    mock_repository.delete.return_value = True
    result = molecule_service.delete(1)
    assert result is True
    mock_repository.delete.assert_called_once_with(1)


def test_delete_molecule_that_does_not_exist(molecule_service, mock_repository):
    mock_repository.find_by_id.return_value = None
    with pytest.raises(UnknownIdentifierException):
        molecule_service.delete(1)


# def test_find_all_molecules(molecule_service, mock_repository):
#     mock_repository.find_all.return_value = [Molecule(**molecule) for molecule in alkanes.values()]
#     molecules = molecule_service.find_all()
#     assert len(molecules) == 3
#     for molecule in molecules:
#         assert is_equal(molecule, alkanes[molecule["name"]])
#
