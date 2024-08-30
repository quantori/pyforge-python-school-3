from typing import Annotated, Optional
from sqlalchemy.orm import Mapped, mapped_column
from src.database import Base
from src.molecules.schema import MoleculeResponse
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
        links = {
            "self": {
                "href": f"/molecules/{self.molecule_id}",
                "rel": "self",
                "type": "GET",
            },
            "substructures": {
                "href": f"/substructure_search?smiles={self.smiles}",
                "rel": "substructures",
                "type": "GET",
            },
        }

        return MoleculeResponse(links=links, **self.__dict__)

    def to_chem(self):
        return Chem.MolFromSmiles(self.smiles)
