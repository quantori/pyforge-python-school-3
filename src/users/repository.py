from functools import lru_cache
from typing import Annotated
from fastapi import Depends
from sqlalchemy.orm import sessionmaker
from src.database import get_session_factory
from src.repositories import SQLAlchemyRepository
from src.users.models import User


class UserRepository(SQLAlchemyRepository):

    def __init__(self, session_factory):
        super().__init__(User, session_factory)


@lru_cache
def get_user_repository(
    session_factory: Annotated[sessionmaker, Depends(get_session_factory)]
):
    return UserRepository(session_factory)
