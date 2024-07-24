from pydantic import BaseModel

class Molecules(BaseModel):
    identifier: str
    smile: str