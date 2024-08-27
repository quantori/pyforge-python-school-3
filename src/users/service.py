from functools import lru_cache
from typing import Annotated
from fastapi import Depends
from fastapi.security import SecurityScopes
from jwt import InvalidTokenError

from src.exceptions import UnknownIdentifierException
from src.users.exceptions import (
    CredentialsException,
    DuplicateEmailException,
    EmailNotFoundException,
    NotEnoughPermissionException,
)
from src.users.models import User
from src.users.repository import UserRepository, get_user_repository
from src.users.schemas import (
    RegisterRequest,
    RegistrationResponse,
    UserResponse,
    UserRequest,
)
from src.users.security import (
    verify_password,
    decode_access_token,
    create_access_token,
    get_password_hash,
    oauth2_scheme,
    role_scopes,
    get_superuser_email_and_password,
)
from src.utils import mark_as_tested


class UserService:
    def __init__(self, user_repository: UserRepository):
        self._repository = user_repository

    @mark_as_tested
    def register(self, register_request: RegisterRequest) -> RegistrationResponse:
        """
        :raises DuplicateEmailException: if the email is already registered
        """
        if self.exists_by_email(register_request.email):
            raise DuplicateEmailException(register_request.email)
        data = register_request.model_dump()
        data["password"] = get_password_hash(data["password"])
        user = self._repository.save(data)
        return RegistrationResponse.model_validate({"email": user.email})

    @mark_as_tested
    def find_by_id(self, user_id) -> UserResponse:
        """
        :raises UnknownIdentifierException: if the user is not found
        """
        if not self.exists_by_id(user_id):
            raise UnknownIdentifierException(user_id)
        user = self._repository.find_by_id(user_id)
        return user.to_response()

    def create_super_admin(self, register_request: RegisterRequest):
        print("Creating super admin")
        if self.exists_by_email(register_request.email):
            print("Super admin already exists")
        self.register(register_request)

    def initiate(self):
        email, password = get_superuser_email_and_password()
        if not self.exists_by_email(email):
            self.create_super_admin(
                RegisterRequest(
                    email=email,
                    password=password,
                    full_name="Super Admin",
                    role="SUPER_ADMIN",
                )
            )

    def update(self, user_id, user_update_request: UserRequest):
        raise NotImplementedError()

    def delete(self, user_id):
        raise NotImplementedError()

    def find_by_email(self, email: str) -> UserResponse:
        """
        :raises: EmailNotFoundException if the email is not found
        """
        found = self._repository.filter(email=email)
        if len(found) == 0:
            raise EmailNotFoundException(email)
        resp = found[0].to_response()
        return resp

    def exists_by_id(self, user_id: int) -> bool:
        return self._repository.find_by_id(user_id) is not None

    def exists_by_email(self, email: str) -> bool:
        return len(self._repository.filter(email=email)) > 0

    def find_all(self, page: int = 0, page_size: int = 1000):
        return self._repository.find_all(page, page_size)

    def authenticate_user(self, username: str, password: str) -> User:
        """
        :param username: email
        :param password: password
        :return: Return the **database entity**, not response model,  if the user is authenticated, None otherwise
        :raises: CredentialsException if the user is not found or the password is incorrect
        """
        print("Authenticating user", username, password)
        found = self._repository.filter(email=username)
        if len(found) == 0:
            raise CredentialsException(message="User not found.")

        user = found[0]

        if not verify_password(password, user.password):
            raise CredentialsException(message="Invalid password.")

        return user

    def get_user_and_scopes_from_token(self, token: str):
        try:
            payload = decode_access_token(token)
            username = payload.get("sub")
            if username is None:
                raise CredentialsException("Invalid token.")
            sop = payload.get("scopes")

        except InvalidTokenError:
            raise CredentialsException("Invalid token.")

        try:
            user = self.find_by_email(username)
        except EmailNotFoundException:
            raise CredentialsException("Invalid token.")

        return user, sop

    def login(self, username: str, password: str, scopes=None):
        if scopes is None:
            scopes = []
        user = self.authenticate_user(username, password)
        if user is None:
            raise CredentialsException("Invalid credentials.")
        user_scopes = role_scopes[user.role]
        for scope in scopes:
            if scope not in user_scopes:
                raise NotEnoughPermissionException(scope)
        token = create_access_token({"sub": user.email, "scopes": scopes})
        return {"access_token": token, "token_type": "bearer"}


@lru_cache
def get_user_service(
    user_repository: Annotated[UserRepository, Depends(get_user_repository)]
):
    return UserService(user_repository)


def get_current_user(
    scopes: SecurityScopes,
    token: Annotated[str, Depends(oauth2_scheme)],
    service: Annotated[UserService, Depends(get_user_service)],
):
    user, sop = service.get_user_and_scopes_from_token(token)
    for scope in scopes.scopes:
        if scope not in sop:
            raise NotEnoughPermissionException(scope)
    return user


def get_current_active_user(
    current_user: User = Depends(get_current_user),
):
    if not current_user.is_active:
        raise CredentialsException("Inactive user")
    return current_user
