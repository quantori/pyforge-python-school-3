
class CollectionAlreadyExistsException(Exception):

    def __init__(self, name: str):
        self.message = f"Collection with name '{name}' already exists."
        super().__init__(self.message)


class NoSuchCollectionException(Exception):

    def __init__(self, name: str):
        self.message = f"Collection with name '{name}' does not exist."
        super().__init__(self.message)


