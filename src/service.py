from src.repositories import MoleculeRepository
from src.schemas import MoleculeRequest


class MoleculeService:
    def __init__(self, repository: MoleculeRepository):
        self._repository = repository

    def find_by_id(self, obj_id: int):
        return self._repository.find_by_id(obj_id)

    def save(self, molecule_request: MoleculeRequest):
        return self._repository.save(molecule_request.dict())

    def find_all(self):
        return self._repository.find_all()
