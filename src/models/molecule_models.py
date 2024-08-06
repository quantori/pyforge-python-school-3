from src.repository.abstract_repository import RepositoryItem, KEY


class Molecule(RepositoryItem[int]):

    def __init__(self, smiles: str, molecule_id: int = 0, molecule_name: str | None = None, description: str | None = None):
        self.smiles: str = smiles
        self.molecule_id: int = molecule_id
        self.molecule_name: str | None = molecule_name
        self.description: str | None = description

    def get_id(self) -> int:
        return self.molecule_id

    def set_id(self, new_id: int) -> None:
        self.molecule_id = new_id
