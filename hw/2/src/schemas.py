from typing import List, Optional, Generic, TypeVar
from pydantic import BaseModel, Field


T = TypeVar('T')
class MoleculeSchema(BaseModel):
    identifier: str
    smile: str
    
    class Config:
        orm_mode = True


class RequestMolecule(BaseModel):
    parameters: MoleculeSchema = Field(...)
    
class Response(BaseModel, Generic[T]):
    code: str
    status: str
    message: str
    result: Optional[T]
    
