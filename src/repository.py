from abc import ABC, abstractmethod
from src.models import Molecule


class MoleculesRepository(ABC):

    @abstractmethod
    def find_by_id(self, molecule_id: str) -> Molecule | None:
        pass

    @abstractmethod
    def find_all(self) -> list[Molecule]:
        pass

    @abstractmethod
    def exists_by_id(self, molecule_id: str) -> bool:
        pass

    @abstractmethod
    def add(self, molecule: Molecule) -> Molecule:
        pass

    @abstractmethod
    def delete_by_id(self, molecule_id: str) -> None:
        pass


class InMemoryMoleculesRepository(MoleculesRepository):

    def __init__(self):
        self._molecules: dict[str, Molecule] = {}

    def find_by_id(self, molecule_id: str) -> Molecule | None:
        return self._molecules.get(molecule_id, None)

    def find_all(self) -> list[Molecule]:
        return list(self._molecules.values())

    def exists_by_id(self, molecule_id: str) -> bool:
        return molecule_id in self._molecules

    def add(self, molecule: Molecule) -> Molecule:
        self._molecules[molecule.molecule_id] = molecule
        return molecule

    def delete_by_id(self, molecule_id: str) -> None:
        if molecule_id in self._molecules:
            del self._molecules[molecule_id]

