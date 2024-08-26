from pydantic import BaseModel

from src.users.security import Role


class Token(BaseModel):
    access_token: str
    token_type: str


class TokenData(BaseModel):
    username: str | None = None


class RegisterRequest(BaseModel):
    username: str
    password: str
    full_name: str
    role: Role

    model_config = {
        "json_schema_extra": {
            "examples": [
                {
                    "username": "gagamagaria",
                    "password": "gaioz2003",
                    "full_name": "Gaioz Tabatadze",
                    "role": "super_admin"
                }
            ]
        }
    }


class UserResponse(BaseModel):
    username: str
    full_name: str
