from pydantic import BaseModel

class MoleculeBase(BaseModel):
    identifier: str
    smile: str

    class Config:
        orm_mode = True
