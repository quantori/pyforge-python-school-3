from __future__ import annotations
import enum
from typing import List, Optional
from sqlalchemy import ForeignKey, Float
from sqlalchemy.orm import Mapped, mapped_column, relationship
from src.database import Base


class QuantityUnit(enum.Enum):
    MOLAR = "MOLAR"
    MASS = "GRAM"
    VOLUME = "ML"


class DrugMolecule(Base):
    __tablename__ = "drug_molecule"

    drug_id: Mapped[int] = mapped_column(ForeignKey("drugs.drug_id"), primary_key=True)
    molecule_id: Mapped[int] = mapped_column(
        ForeignKey("molecules.molecule_id"), primary_key=True
    )
    quantity: Mapped[float] = mapped_column(Float, nullable=False)
    quantity_unit: Mapped[QuantityUnit] = mapped_column(nullable=False)
    drug = relationship("Drug", back_populates="molecules")
    # Define relationships if needed
    # Example: molecule = relationship("Molecule", back_populates="drug_molecules")


class Drug(Base):
    __tablename__ = "drugs"

    drug_id: Mapped[int] = mapped_column(primary_key=True, autoincrement=True)
    name: Mapped[str] = mapped_column(nullable=False)
    description: Mapped[Optional[str]] = mapped_column(nullable=True)

    # Define the relationship with DrugMolecule
    molecules: Mapped[List[DrugMolecule]] = relationship(
        "DrugMolecule", back_populates="drug", cascade="all, " "delete-orphan"
    )

    # Define relationships if needed
    # Example: drug_molecules = relationship("DrugMolecule", back_populates="drug")
