from functools import lru_cache
from src.molecules.model import Molecule
from src.repository import SQLAlchemyRepository


class MoleculeRepository(SQLAlchemyRepository):
    def __init__(self):
        super().__init__(Molecule)

    # def bulk_save(self, molecules: list[MoleculeRequest]):
    #     session = self._get_session()
    #     ans = session.scalars(
    #         insert(Molecule).returning(Molecule.molecule_id), molecules
    #     )
    #     session.commit()
    #     session.close()
    #     return list(ans)


@lru_cache
def get_molecule_repository():
    return MoleculeRepository()
