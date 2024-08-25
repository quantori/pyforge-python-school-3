from typing import Union


class RBMolecule:
    def __init__(
        self,
        id: Union[int, None] = None,
        name: Union[str, None] = None,
    ):
        self.id = id
        self.name = name

    def to_dict(self) -> dict:
        data = {
            "id": self.id,
            "name": self.name,
        }
        filtered_data = {key: value for key, value in data.items() if value is not None}
        return filtered_data