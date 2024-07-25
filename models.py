from pydentic import BaseModel

class User(BaseModel):
    user_id: int
    name: str
    location: str