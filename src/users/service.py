from functools import lru_cache
from typing import Annotated
from fastapi import Depends
from jwt import InvalidTokenError
from src.users.exceptions import CredentialsException
from src.users.models import User
from src.users.repository import UserRepository, get_user_repository
from src.users.schemas import RegisterRequest
from src.users.security import (
    verify_password,
    decode_access_token,
    create_access_token,
    get_password_hash,
)


class UserService:
    def __init__(self, user_repository: UserRepository):
        self._repository = user_repository

    def register(self, register_request: RegisterRequest):
        # TODO validate username
        data = register_request.model_dump()
        data["password"] = get_password_hash(data["password"])
        user = self._repository.save(data)
        return user.to_response()

    def find_by_username(self, username: str):
        pass

    def find_by_id(self, user_id):
        return self._repository.find_by_id(user_id).to_response()

    def authenticate_user(self, username: str, password: str) -> User | None:
        """
        :param username:
        :param password:
        :return: Return the **database entity**, not response model,  if the user is authenticated, None otherwise
        """
        user = self._repository.find_by_id(username)

        if not user:
            return None

        if not verify_password(password, user.password):
            return None

        return user

    def get_current_user(self, token: str):
        try:
            payload = decode_access_token(token)
            username = payload.get("sub")
            if username is None:
                raise CredentialsException()
        except InvalidTokenError:
            raise CredentialsException()
        user = self.find_by_id(username)
        if user is None:
            raise CredentialsException()
        return user

    def login(self, username: str, password: str):
        user = self.authenticate_user(username, password)
        if user is None:
            raise CredentialsException()
        token = create_access_token({"sub": user.username})
        return {"access_token": token, "token_type": "bearer"}


@lru_cache
def get_user_service(
    user_repository: Annotated[UserRepository, Depends(get_user_repository)]
):
    return UserService(user_repository)
