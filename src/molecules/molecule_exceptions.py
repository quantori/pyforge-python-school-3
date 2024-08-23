from src.exceptions import BadRequestException


class InvalidSmilesException(BadRequestException):
    def __init__(self, smiles):
        self.smiles = smiles
        self.message = (
            f"Smiles string {self.smiles} does not represent a valid molecule"
        )
        super().__init__(self.message)


class DuplicateSmilesException(BadRequestException):
    def __init__(self, smiles):
        self.smiles = smiles
        self.message = f"Smiles string {self.smiles} is not unique"
        super().__init__(self.message)
