from functools import lru_cache
from src.configs import Settings
from src.repositories import MoleculeRepository
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
from src.service import MoleculeService
from src.schemas import PaginationQueryParams


@lru_cache
def get_database_url():
    return Settings().database_url


@lru_cache
def get_session_factory():
    return sessionmaker(bind=create_engine(get_database_url()))


@lru_cache()
def get_molecule_repository():
    return MoleculeRepository(get_session_factory())


@lru_cache
def get_molecule_service():
    return MoleculeService(get_molecule_repository())


@lru_cache
def get_pagination_query_params(page: int = 0, page_size: int = 1000):
    return PaginationQueryParams(page=page, page_size=page_size)
