from functools import lru_cache

from sqlalchemy.orm import sessionmaker

from src.database import get_session_factory
from src.molecules.models import Molecule
from src.repositories import SQLAlchemyRepository


class MoleculeRepository(SQLAlchemyRepository):
    def __init__(self, session_factory: sessionmaker):
        super().__init__(Molecule, session_factory)


@lru_cache
def get_molecule_repository():
    return MoleculeRepository(get_session_factory())
