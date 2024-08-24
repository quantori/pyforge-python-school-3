from pydantic import BaseModel


class Token(BaseModel):
    access_token: str
    token_type: str


class TokenData(BaseModel):
    username: str | None = None


class RegisterRequest(BaseModel):
    username: str
    password: str
    full_name: str

    model_config = {
        "json_schema_extra": {
            "examples": [
                {
                    "username": "gagamagaria",
                    "password": "gaioz2003",
                    "full_name": "Gaioz Tabatadze",
                }
            ]
        }
    }


class UserResponse(BaseModel):
    username: str
    full_name: str
