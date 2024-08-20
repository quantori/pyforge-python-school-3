class InvalidSmilesException(Exception):
    def __init__(self, smiles):
        self.smiles = smiles
        self.message = (
            f"Smiles string {self.smiles} does not represent a valid molecule"
        )
        super().__init__(self.message)


class UnknownIdentifierException(Exception):
    def __init__(self, identifier):
        self.identifier = identifier
        self.message = f"Unknown identifier {self.identifier}"
        super().__init__(self.message)
