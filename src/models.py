from sqlalchemy import Column, Integer, String, Float
from .database import Base


class Molecule(Base):
    __tablename__ = 'molecules'

    id = Column(Integer, primary_key=True, index=True)
    name = Column(String, index=True)
    smiles = Column(String, index=True)
    weight = Column(Float)
    formula = Column(String)
