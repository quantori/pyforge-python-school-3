from functools import lru_cache
from src.configs import Settings
from src.repositories import MoleculeRepository
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
from src.service import MoleculeService


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
