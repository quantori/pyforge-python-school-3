import pytest
from sqlalchemy import create_engine, text
from sqlalchemy.orm import sessionmaker
from starlette.testclient import TestClient
from src.config import get_settings
from src.drugs.repository import DrugRepository
from src.drugs.service import DrugService
from src.drugs.service import get_drug_service
from src.main import app
from src.molecules.repository import MoleculeRepository
from src.molecules.service import MoleculeService, get_molecule_service
from src.drugs.tests import sample_data
from src.database import Base

engine = create_engine(get_settings().TEST_DB_URL)

session_factory = sessionmaker(bind=engine)
drug_repository = DrugRepository()
drug_service = DrugService(drug_repository, session_factory)

mol_repository = MoleculeRepository()
molecule_service = MoleculeService(mol_repository, session_factory)

test_client = TestClient(app)

app.dependency_overrides[get_drug_service] = lambda: drug_service
app.dependency_overrides[get_molecule_service] = lambda: molecule_service


def sth():
    Base.metadata.drop_all(engine)
    Base.metadata.create_all(engine)
    molecule_service.save(sample_data.caffeine_request)
    molecule_service.save(sample_data.water_request)
    molecule_service.save(sample_data.sugar_request)
    drug_service.save(sample_data.coffe_request)


def select():
    with session_factory() as session:
        print("molecules \n")
        result = session.execute(text("SELECT * FROM molecules"))
        for row in result:
            print(row)

        print("drugs \n")
        result = session.execute(text("SELECT * FROM drugs"))
        for row in result:
            print(row)
        result = session.execute(text("SELECT * FROM drug_molecule"))

        print("drug_molecules \n")
        for row in result:
            print(row)

        print("\n")


drug = drug_service.find_by_id(1)
for mol in drug.molecules:
    print(mol.molecule_id)
    print(mol.quantity)
    print(mol.quantity_unit)
    print("\n")

# sth()
# select()
# with session_factory() as session:
#     session.execute(text("DELETE FROM molecules where name = 'Caffeine'"))


