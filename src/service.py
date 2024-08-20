from src.exceptions import UnknownIdentifierException, DuplicateSmilesException
from src.repositories import MoleculeRepository
from src.schemas import MoleculeRequest


class MoleculeService:

    def __init__(self, repository: MoleculeRepository):
        self._repository = repository

    def exists_by_id(self, obj_id: int):
        """
        Check if a molecule with the given id exists.
        :param obj_id: molecule id
        :return: True if the molecule exists, False otherwise
        """
        return self._repository.find_by_id(obj_id) is not None

    def find_by_id(self, obj_id: int):
        """
        Find a molecule by its id. Calls exists_by_id to check if the molecule exists, resulting in two database calls.
        Not vert impressive, but I am trying to keep it simple.

        :param obj_id:  molecule id
        :return: found molecule
        :raises UnknownIdentifierException: if the molecule with the given id does not exist
        """

        if not self.exists_by_id(obj_id):
            raise UnknownIdentifierException(obj_id)
        return self._repository.find_by_id(obj_id)

    def save(self, molecule_request: MoleculeRequest):
        """
        Simply save a new molecule to the database. If the smiles is not unique, the database will raise an exception.

        :param molecule_request: Molecule data
        :return: Saved molecule
        """
        same_smiles = self._repository.filter(smiles=molecule_request.smiles)
        if same_smiles:
            raise DuplicateSmilesException(molecule_request.smiles)
        return self._repository.save(molecule_request.dict())

    def update(self, obj_id: int, molecule_request: MoleculeRequest):
        """
        Update a molecule with the given id.
        This is suitable for put request

        :param obj_id: Identifier of the molecule to be updated
        :param molecule_request: New data for the molecule
        :return: Updated molecule
        :raises UnknownIdentifierException: if the molecule with the given id does not exist
        """
        if not self.exists_by_id(obj_id):
            raise UnknownIdentifierException(obj_id)
        return self._repository.update(obj_id, molecule_request.dict())

    def find_all(self, page: int = 0, page_size: int = 1000):
        """
        Find all molecules in the database. Can be paginated. Default page size is 1000.

        :param page: Zero indexed page number, default is 0
        :param page_size: Items per page, default is 1000
        :return: List of all molecules
        """
        return self._repository.find_all(page, page_size)

    def delete(self, obj_id: int):
        """
        Delete a molecule with the given id. If the molecule does not exist, raise an exception.

        :param obj_id: Identifier of the molecule to be deleted
        :return: True
        :raises UnknownIdentifierException: if the molecule with the given id does not exist
        """
        if not self.exists_by_id(obj_id):
            raise UnknownIdentifierException(obj_id)
        return self._repository.delete(obj_id)
