from typing import Annotated, Optional
from sqlalchemy.orm import Mapped, mapped_column
from src.database import Base


class Molecule(Base):
    molecule_id: Mapped[
        Annotated[int, mapped_column(primary_key=True, autoincrement=True)]
    ]
    smiles: Mapped[Annotated[str, mapped_column(unique=True)]]
    name: Mapped[Annotated[Optional[str], mapped_column()]]

    def __repr__(self):
        return f"Molecule(molecule_id={self.molecule_id}, smiles={self.smiles})"
