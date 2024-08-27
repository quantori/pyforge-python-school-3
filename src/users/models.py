from dataclasses import dataclass
from typing import Annotated

from sqlalchemy.orm import mapped_column, Mapped

from src.database import Base
from src.users.schemas import UserResponse
from src.users.security import Role


@dataclass
class User(Base):
    __tablename__ = "users_"

    user_id: Mapped[Annotated[int, mapped_column(primary_key=True, autoincrement=True)]]
    email: Mapped[Annotated[str, mapped_column(unique=True)]]
    full_name: Mapped[Annotated[str, mapped_column(nullable=False)]]
    password: Mapped[Annotated[str, mapped_column(nullable=False)]]
    role: Mapped[Annotated[Role, mapped_column(nullable=False)]]
    is_active: Mapped[Annotated[bool, mapped_column(nullable=False, default=True)]]

    def to_response(self) -> UserResponse:
        return UserResponse.model_validate(
            {
                "user_id": self.user_id,
                "email": self.email,
                "is_active": self.is_active,
                "role": self.role,
                "full_name": self.full_name,
            }
        )
