from dataclasses import dataclass
from typing import Annotated

from pydantic import BaseModel, EmailStr, Field

from src.users.security import Role


class Token(BaseModel):
    access_token: str
    token_type: str


class TokenData(BaseModel):
    username: str | None = None


class RegisterRequest(BaseModel):
    email: Annotated[EmailStr, Field(..., )]
    password: Annotated[str, Field(..., min_length=6)]
    full_name: Annotated[str, Field(..., min_length=3)]
    role: Annotated[Role, Field(..., )]

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


class UserRequest(BaseModel):
    """
    Used for updating user data by higher level users
    """

    email: Annotated[EmailStr, Field(..., )]
    is_active: Annotated[bool, Field(..., )]
    role: Annotated[Role, Field(..., )]
    full_name: Annotated[str, Field(..., )]
    password: Annotated[str, Field(..., )]


class UserResponse(BaseModel):
    user_id: Annotated[int, Field(..., )]
    email: Annotated[EmailStr, Field(..., )]
    is_active: Annotated[bool, Field(..., )]
    role: Annotated[Role, Field(..., )]
    full_name: Annotated[str, Field(..., )]

    model_config = {
        "from_attributes": True,
        "extra": "ignore",
    }


@dataclass
class RegistrationResponse(BaseModel):
    email: Annotated[EmailStr, Field(..., )]

    model_config = {
        "from_attributes": True,
        "extra": "allow",
    }
