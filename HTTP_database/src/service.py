import os
import json
import HTTP_database.src.exception as exception


class CollectionService:
    """
    TODO: This implementation is not efficient. We don't use memory to store the collections at all.
    TODO: Currently every deletion removes the entry from the file and shifts the rest of the entries. This is bad. Fix it.
    """
    def __init__(self, base_path: str):
        self.base_path = base_path
        self.collections_dir_path = os.path.join(self.base_path, "collections")
        self.collections_file_path = os.path.join(self.collections_dir_path, "collections.json")

        # Ensure the given base path exists
        if not os.path.exists(self.base_path):
            raise FileNotFoundError(f"The directory {self.base_path} does not exist.")

        # Create the collections directory if it doesn't exist
        if not os.path.exists(self.collections_dir_path):
            os.makedirs(self.collections_dir_path)

        # Create the collections.json file if it doesn't exist
        if not os.path.exists(self.collections_file_path):
            with open(self.collections_file_path, "w") as file:
                file.write("")

    def exists_by_name(self, collection_name: str):
        # at this point we are sure that the collections.json file exists
        with open(self.collections_file_path, "r") as file:
            for line in file:
                if line.strip():  # Ensure the line is not empty
                    collection = json.loads(line)
                    print(collection["name"])
                    if collection["name"] == collection_name:
                        return True
        return False

    def create_collection(self, collection_name: str):
        if self.exists_by_name(collection_name):
            raise exception.CollectionAlreadyExistsException(collection_name)

        # This should also create the file with the same name as the collection
        with open(os.path.join(self.collections_dir_path, collection_name+".json"), "w") as file:
            file.write("")

        with open(self.collections_file_path, "a") as file:
            file.write(json.dumps({"name": collection_name, "size": 0}) + "\n")

    def get_collection(self, collection_name: str):
        with open(self.collections_file_path, "r") as file:
            for line in file:
                if line.strip():
                    collection = json.loads(line)
                    if collection["name"] == collection_name:
                        return collection
        raise exception.NoSuchCollectionException(collection_name)

    def get_collections(self):
        collections = []
        with open(self.collections_file_path, "r") as file:
            for line in file:
                if line.strip():
                    collections.append(json.loads(line))
        return collections

    def delete_collection(self, collection_name: str):
        if not self.exists_by_name(collection_name):
            raise exception.NoSuchCollectionException(collection_name)

        collections = self.get_collections()
        with open(self.collections_file_path, "w") as file:
            for collection in collections:
                if collection["name"] != collection_name:
                    file.write(json.dumps(collection) + "\n")

        # Delete the file associated with the collection
        os.remove(os.path.join(self.collections_dir_path, collection_name+".json"))
        return None
