from typing import Type

from sqlalchemy.orm import sessionmaker

from src.database import Base
from src.models import Molecule


class SQLAlchemyRepository:
    def __init__(self, model_type: Type[Base], session_factory: sessionmaker):
        self._model_type = model_type
        self._session_factory = session_factory

    def find_by_id(self, obj_id):
        with self.__get_session() as session:
            return session.query(self._model_type).get(obj_id)

    def save(self, data: dict):
        instance = self._model_type(**data)
        with self.__get_session() as session:
            session.add(instance)
            session.commit()
            session.refresh(instance)
        return instance

    def find_all(self):
        with self.__get_session() as session:
            return session.query(self._model_type).all()

    def __get_session(self):
        return self._session_factory()


class MoleculeRepository(SQLAlchemyRepository):
    def __init__(self, session_factory: sessionmaker):
        super().__init__(Molecule, session_factory)
