from typing import Union


class RBMolecule:
    def __init__(
        self,
        mol_id: Union[int, None] = None,
        name: Union[int, None] = None,
    ):
        self.id = mol_id
        self.name = name

    def to_dict(self) -> dict:
        data = {
            "id": self.id,
            "name": self.name,
        }
        filtered_data = {key: value for key, value in data.items() if value is not None}
        return filtered_data