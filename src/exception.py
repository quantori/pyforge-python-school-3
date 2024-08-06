class UnknownIdentifierException(Exception):
    def __init__(self, identifier):
        super().__init__(f"Unknown identifier: {identifier}")
        self.identifier = identifier


class RepositoryItemNotFountException(Exception):
    def __init__(self, key):
        super().__init__(f"Repository item with the key: {key} is not found")
        self.key = key


class InvalidSmilesException(Exception):
    def __init__(self, smiles):
        super().__init__(f"Invalid SMILES string: {smiles}")
        self.smiles = smiles
