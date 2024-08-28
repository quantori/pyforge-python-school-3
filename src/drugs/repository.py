from functools import lru_cache

from sqlalchemy.orm import sessionmaker

from src.database import get_session_factory
from src.drugs.model import Drug, DrugMolecule
from src.repositories import SQLAlchemyRepository


class DrugRepository(SQLAlchemyRepository):
    def __init__(self, session_factory: sessionmaker):
        super().__init__(Drug, session_factory)

    def save(self, data: dict):
        """
        First create a drug, then drug_molecules. they will automatically fail if molecule_id is not found

        :param data:
        :return:
        """
        session = self._get_session()
        drug = Drug(name=data["name"], description=data.get("description"))
        session.add(drug)
        session.flush()
        for molecule in data["molecules"]:
            drug_molecule = DrugMolecule(
                drug_id=drug.drug_id,
                molecule_id=molecule["molecule_id"],
                quantity=molecule["quantity"],
                quantity_unit=molecule["quantity_unit"],
            )
            session.add(drug_molecule)

        session.commit()
        session.refresh(drug)
        session.close()
        return drug


@lru_cache
def get_drug_repository():
    return DrugRepository(get_session_factory())
