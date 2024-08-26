from typing import Annotated
from fastapi import APIRouter, Depends
from fastapi.security import OAuth2PasswordRequestForm
from src.users.schemas import RegisterRequest
from src.users.service import UserService, get_user_service

router = APIRouter()


@router.post("/register")
def register(
    register_request: RegisterRequest,
    service: Annotated[UserService, Depends(get_user_service)],
):
    return service.register(register_request)


@router.post("/token")
def login_for_access_token(
    form_data: Annotated[OAuth2PasswordRequestForm, Depends()],
    service: Annotated[UserService, Depends(get_user_service)],
):
    return service.login(form_data.username, form_data.password, form_data.scopes)


# @router.get("/me")
# def get_current_user(
#     token: Annotated[str, Depends(oauth2_scheme)],
#     service: Annotated[UserService, Depends(get_user_service)],
# ):
#     return service.get_current_user(token)
