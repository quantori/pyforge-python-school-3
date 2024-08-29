from datetime import datetime
from functools import lru_cache
from typing import Annotated
from sqlalchemy import func, create_engine
from sqlalchemy.orm import (
    DeclarativeBase,
    declared_attr,
    mapped_column,
    Mapped,
    sessionmaker,
)

from src.configs import Settings

created_at = Annotated[datetime, mapped_column(server_default=func.now())]
updated_at = Annotated[
    datetime, mapped_column(server_default=func.now(), onupdate=datetime.now)
]


class Base(DeclarativeBase):
    __abstract__ = True

    @declared_attr.directive
    def __tablename__(cls) -> str:
        return f"{cls.__name__.lower()}s"

    created_at: Mapped[created_at]
    updated_at: Mapped[updated_at]


@lru_cache
def get_database_url():
    return Settings().database_url


@lru_cache
def get_session_factory():
    return sessionmaker(bind=create_engine(get_database_url()))

