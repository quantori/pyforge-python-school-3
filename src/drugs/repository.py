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

        I changed the usual behavour of the save method to return the drug object and the session object.
        This session is used in the service to close the session after the transaction is done.
        This became necessary because if I close the session in the repository, the service will not be able to use
        drug. molecules to get the molecules of the drug, as they are lazy loaded.

        I realize now why it's better to pass the session object to every repository methods, so that the service can control the
        session lifecycle.

        :param data:
        :return:
        """
        session = super()._get_session()
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
        return drug, session


@lru_cache
def get_drug_repository():
    return DrugRepository(get_session_factory())
