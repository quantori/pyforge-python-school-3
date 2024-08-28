from typing import Type
from sqlalchemy import select
from sqlalchemy.orm import sessionmaker
from src.database import Base


class SQLAlchemyRepository:
    """
    Base class for all repositories that use SQLAlchemy.
    Updates should be done at the service level.
    https://docs.sqlalchemy.org/en/20/tutorial/orm_data_manipulation.html#updating-orm-objects-using-the-unit-of-work-pattern
    """

    def __init__(self, model_type: Type[Base], session_factory: sessionmaker):
        self._model_type = model_type
        self._session_factory = session_factory

    def find_by_id(self, obj_id):
        with self._get_session() as session:
            return session.get(self._model_type, obj_id)

    def find_all(self, page=0, page_size=1000):
        """

        Find all instances of the model with pagination support

        :param page: Zero indexed page number
        :param page_size: Items per page
        :return: List of instances
        """

        session = self._get_session()
        stmt = select(self._model_type).limit(page_size).offset(page * page_size)

        result = session.execute(stmt).scalars().all()
        session.flush()
        session.close()
        return result

    def filter(self, **kwargs):
        with self._get_session() as session:
            stmt = select(self._model_type).filter_by(**kwargs)
            ans = session.execute(stmt).scalars().all()
            return ans

    def save(self, data: dict):
        session = self._get_session()
        instance = self._model_type(**data)
        session.add(instance)
        session.commit()
        session.refresh(instance)
        session.close()
        return instance

    def update(self, obj_id, data: dict):
        session = self._get_session()
        instance = session.get(self._model_type, obj_id)
        for key, value in data.items():
            setattr(instance, key, value)
        session.commit()
        session.refresh(instance)
        session.close()
        return instance

    def delete(self, obj_id):
        """
        Delete an instance
        :param obj_id: id of the instance to be deleted
        :return:  True if the instance is deleted, False otherwise
        """

        try:
            session = self._get_session()
            instance = session.get(self._model_type, obj_id)
            session.delete(instance)
            session.commit()
            session.close()
            return True
        except Exception:
            return False

    def _get_session(self):
        return self._session_factory()
