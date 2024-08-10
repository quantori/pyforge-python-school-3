import json

from fastapi import Depends
from src import exceptions
from HTTP_database_client.src.database_client import HTTPDatabaseClient

import logging
from abc import ABC, abstractmethod
from src.models.molecule_models import MoleculeInDB


class IDGenerator(ABC):

    @abstractmethod
    def get_next_id(self) -> int:
        pass


class AbstractMoleculeRepository(ABC):

    def __init__(self, id_generator: IDGenerator):
        self._id_generator = id_generator

    @abstractmethod
    def find_by_id(self, obj_id: int) -> MoleculeInDB:
        pass

    @abstractmethod
    def find_all(self) -> list[MoleculeInDB]:
        pass

    @abstractmethod
    def exists_by_id(self, obj_id: int) -> bool:
        pass

    @abstractmethod
    def add(self, obj: MoleculeInDB) -> MoleculeInDB:
        pass

    @abstractmethod
    def delete_by_id(self, obj_id: int) -> None:
        pass

    @abstractmethod
    def size(self) -> int:
        pass

    @abstractmethod
    def clear(self) -> None:
        pass


class DefaultMoleculeIDGenerator(IDGenerator):
    def __init__(self):
        self.__next_id = 0

    def get_next_id(self) -> int:
        next_id = self.__next_id
        self.__next_id += 1
        return next_id


class InMemoryMoleculesRepository(AbstractMoleculeRepository):
    def __init__(self, id_generator: IDGenerator = DefaultMoleculeIDGenerator()):
        super().__init__(id_generator)
        self.molecules = {}

    def find_by_id(self, obj_id: int) -> MoleculeInDB:
        """
        :return:  Model object if found, None otherwise.
        """
        if not self.exists_by_id(obj_id):
            raise exception.RepositoryItemNotFountException(obj_id)
        return self.molecules.get(obj_id)

    def find_all(self) -> list[MoleculeInDB]:
        """
            Its important that the insertion order is preserved.
            Inner implementation of dictionary is ordered in that way since Python 3.7+
        """
        return list(self.molecules.values())

    def exists_by_id(self, obj_id: int) -> bool:
        return obj_id in self.molecules

    def add(self, obj: MoleculeInDB) -> MoleculeInDB:
        """
            If the provided molecule has non-None molecule_id, it will be replaced with a new one.
        """
        obj.set_id(self._id_generator.get_next_id())
        self.molecules[obj.get_id()] = obj
        return obj

    def delete_by_id(self, obj_id: int) -> None:
        """
        :raises RepositoryItemNotFountException: If the object with the provided id does not exist.
        """
        if not self.exists_by_id(obj_id):
            raise exception.RepositoryItemNotFountException(obj_id)
        del self.molecules[obj_id]

    def size(self):
        return len(self.molecules)

    def clear(self) -> None:
        self.molecules.clear()


class HTTPMoleculeRepository(AbstractMoleculeRepository):
    def __init__(self, base_url: str = Depends, id_generator: IDGenerator= DefaultMoleculeIDGenerator()):
        super().__init__(id_generator)
        self.__client = HTTPDatabaseClient(base_url)
        self.__base_url = base_url

        # initialize the repository by creating the molecules collection in the database.
        # If the collection already exists it will not be overwritten.
        response = self.__client.create_collection("molecules")

        if response.status not in [201, 400]:
            raise exception.HTTPClientException(response)

    def find_by_id(self, obj_id: int) -> MoleculeInDB:
        pass

    def find_all(self) -> list[MoleculeInDB]:
        response = self.__client.get_documents("molecules")
        if response.status != 200:
            raise exception.HTTPClientException(response)
        documents = json.loads(response.read())
        return [MoleculeInDB(**document) for document in documents]

    def exists_by_id(self, obj_id: int) -> bool:
        pass

    def add(self, obj: MoleculeInDB) -> MoleculeInDB:
        pass

    def delete_by_id(self, obj_id: int) -> None:
        pass

    def size(self) -> int:
        pass

    def clear(self) -> None:
        pass
