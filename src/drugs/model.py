from __future__ import annotations

from dataclasses import dataclass
from typing import List, Annotated, Optional
from sqlalchemy import ForeignKey
from sqlalchemy.orm import Mapped
from sqlalchemy.orm import mapped_column
from sqlalchemy.orm import relationship
from src.database import Base
from src.drugs.schema import DrugResponse, MolecularConcentrationUnit


class DrugMolecule(Base):
    __tablename__ = "drug_molecule"

    drug_id: Mapped[
        Annotated[int, mapped_column(ForeignKey("drugs.drug_id"), primary_key=True)]
    ]
    molecule_id: Mapped[
        Annotated[
            int, mapped_column(ForeignKey("molecules.molecule_id"), primary_key=True)
        ]
    ]
    molecule_concentration: Mapped[Annotated[float, mapped_column(nullable=False)]]
    molecule_concentration_unit: Mapped[
        Annotated[MolecularConcentrationUnit, mapped_column(nullable=False)]
    ]


@dataclass
class Drug(Base):
    drug_id: Mapped[Annotated[int, mapped_column(primary_key=True, autoincrement=True)]]
    name: Mapped[Annotated[str, mapped_column(nullable=False)]]
    description: Mapped[Annotated[Optional[str], mapped_column(nullable=True)]]

    molecules: Mapped[List[DrugMolecule]] = relationship()

    def to_response(self):
        return DrugResponse(
            id=self.drug_id,
            name=self.name,
            description=self.description,
            molecule_ids=[molecule.molecule_id for molecule in self.molecules],
            links=[],
        )
