import pytest
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker

from src.configs import get_settings
from src.drugs import mapper
from src.drugs.repository import DrugRepository
from src.drugs.service import DrugService
from src.database import Base
from src.molecules.molecule_repository import MoleculeRepository
from src.molecules.service import MoleculeService
import src.drugs.tests.samle_data as sample_data

engine = create_engine(get_settings().TEST_DB_URL)
session_factory = sessionmaker(bind=engine)
repository = DrugRepository(session_factory)
service = DrugService(repository)

molecule_repository = MoleculeRepository(session_factory)
molecule_service = MoleculeService(molecule_repository)


@pytest.fixture
def init_db():
    Base.metadata.drop_all(engine)
    Base.metadata.create_all(engine)
    # add caffeine molecule
    molecule_service.save(sample_data.caffeine_request)
    # add water and sugar molecules
    molecule_service.save(sample_data.water_request)
    molecule_service.save(sample_data.sugar_request)
    yield
    Base.metadata.drop_all(engine)


def test_save(init_db):
    response = service.save(sample_data.coffe_request)
    assert response.name == sample_data.coffe_request.name
    assert response.description == sample_data.coffe_request.description




