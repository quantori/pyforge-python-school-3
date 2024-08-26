from typing import Annotated

from sqlalchemy.orm import mapped_column, Mapped

from src.database import Base
from src.users.schemas import UserResponse
from src.users.security import Role


class User(Base):
    __tablename__ = "users_"

    username: Mapped[Annotated[str, mapped_column(primary_key=True)]]
    full_name: Mapped[Annotated[str, mapped_column(nullable=False)]]
    password: Mapped[Annotated[str, mapped_column(nullable=False)]]
    role: Mapped[Annotated[Role, mapped_column(nullable=False)]]

    def to_response(self):
        return UserResponse(username=self.username, full_name=self.full_name)

