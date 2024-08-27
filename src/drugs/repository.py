from functools import lru_cache

from sqlalchemy.orm import sessionmaker

from src.database import get_session_factory
from src.drugs.model import Drug
from src.repositories import SQLAlchemyRepository


class DrugRepository(SQLAlchemyRepository):
    def __init__(self, session_factory: sessionmaker):
        super().__init__(Drug, session_factory)


@lru_cache
def get_drug_repository():
    return DrugRepository(get_session_factory())
