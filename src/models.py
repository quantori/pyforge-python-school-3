from pydantic import BaseModel, Field


class AddMoleculeRequest(BaseModel):
    smiles: str
    molecule_name: str | None = None
    description: str | None = None

    model_config = {
        "json_schema_extra": {
            "examples": [
                {
                    "smiles": "CC(=O)Oc1ccccc1C(=O)O",
                    "molecule_name": "Aspirin",
                    "description": "Aspirin is in a group of medications called salicylates... so on and so forth."
                }
            ]
        }
    }


class MoleculeResponse(AddMoleculeRequest):
    molecule_id: int
    smiles: str
    molecule_name: str | None
    description: str | None = None

    @staticmethod
    def from_molecule(molecule: 'Molecule') -> 'MoleculeResponse':
        return MoleculeResponse(molecule_id=molecule.molecule_id, smiles=molecule.smiles,
                                molecule_name=molecule.molecule_name, description=molecule.description)

    model_config = {
        "json_schema_extra": {
            "examples": [
                {
                    "molecule_id": 1,
                    "smiles": "CC(=O)Oc1ccccc1C(=O)O",
                    "molecule_name": "Aspirin",
                    "description": "Aspirin is in a group of medications called salicylates... so on and so forth."
                }
            ]
        }
    }


class Molecule:
    __next_id = 0

    def __init__(self, smiles: str, molecule_name: str | None = None, description: str | None = None):
        self._molecule_id: int = Molecule.__next_id
        Molecule.__next_id += 1
        self.molecule_name: str = molecule_name
        self.smiles: str = smiles
        self.description: str | None = description

    @property
    def molecule_id(self) -> int:
        return self._molecule_id

    @staticmethod
    def from_add_molecule_request(add_molecule_request: AddMoleculeRequest) -> 'Molecule':
        return Molecule(smiles=add_molecule_request.smiles, molecule_name=add_molecule_request.molecule_name,
                        description=add_molecule_request.description)
