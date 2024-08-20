from typing import Annotated, Optional
from sqlalchemy.orm import Mapped, mapped_column
from src.database import Base
from src.schemas import MoleculeResponse
from rdkit import Chem


class Molecule(Base):
    molecule_id: Mapped[
        Annotated[int, mapped_column(primary_key=True, autoincrement=True)]
    ]
    smiles: Mapped[Annotated[str, mapped_column(unique=True)]]
    name: Mapped[Annotated[Optional[str], mapped_column()]]

    def __repr__(self):
        return f"Molecule(molecule_id={self.molecule_id}, smiles={self.smiles}, name={self.name})"

    def to_response(self) -> MoleculeResponse:
        return MoleculeResponse(**self.__dict__)

    def to_chem(self):
        return Chem.MolFromSmiles(self.smiles)

