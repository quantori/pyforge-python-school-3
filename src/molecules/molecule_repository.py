from functools import lru_cache
from sqlalchemy.orm import sessionmaker
from src.database import get_session_factory
from src.molecules.models import Molecule
from src.repositories import SQLAlchemyRepository


class MoleculeRepository(SQLAlchemyRepository):
    def __init__(self, session_factory: sessionmaker):
        super().__init__(Molecule, session_factory)

    # def bulk_save(self, molecules: list[MoleculeRequest]):
    #     session = self._get_session()
    #     ans = session.scalars(
    #         insert(Molecule).returning(Molecule.molecule_id), molecules
    #     )
    #     session.commit()
    #     session.close()
    #     return list(ans)


@lru_cache
def get_molecule_repository():
    return MoleculeRepository(get_session_factory())
