from typing import Annotated
from fastapi import APIRouter, Depends, Security
from fastapi.security import OAuth2PasswordRequestForm

from src.users.models import User
from src.users.schemas import RegisterRequest, UserResponse, RegistrationResponse
from src.users.security import Scope
from src.users.service import UserService, get_user_service, get_current_active_user

router = APIRouter()


@router.post(
    "/register",
    status_code=201,
    responses={
        201: {
            "model": RegistrationResponse,
            "description": "User registered successfully",
        },
        400: {"model": str, "description": "probably due to email already registered"},
        403: {
            "model": str,
            "description": "Forbidden. You have to be a super admin to register a user, from this "
            "endpoint. Use patient or doctor registration endpoint instead.",
        },
    },
)
def register(
    _: Annotated[None, Security(get_current_active_user, scopes=[Scope.USERS_WRITE])],
    register_request: RegisterRequest,
    service: Annotated[UserService, Depends(get_user_service)],
) -> RegistrationResponse:
    """
    This endpoint is only for super admin to register a new user. patients and doctors use the respective
    endpoints, that are accessible to all users.

    :param register_request:
    :return RegistrationResponse:
    """
    return service.register(register_request)


@router.post("/token", status_code=201)
def login_for_access_token(
    form_data: Annotated[OAuth2PasswordRequestForm, Depends()],
    service: Annotated[UserService, Depends(get_user_service)],
):
    """

    :param form_data: username means email, and you can ask for various scopes(permissions) such as:

        "molecules:write"

        "drugs:write"

        "users:read"

        "users:write"

        Use openApi documentation for better details

    :param service:
    :return:
    """
    # this line is necessary to initialize the service and create the super admin
    service.initiate()
    return service.login(form_data.username, form_data.password, form_data.scopes)


@router.get("/me")
def get_current_user(
    user: Annotated[User, Security(get_current_active_user, scopes=[Scope.USERS_READ])],
    service: Annotated[UserService, Depends(get_user_service)],
) -> UserResponse:
    """
    Get all the details of the currently logged-in user. This endpoint is accessible to
    higher roles only, like super admin, lab admin, hospital admin.

    For patient and doctor, use the respective endpoints.

    :param user:
    :param token:
    :param service:
    :return:
    """
    return service.find_by_id(user.user_id)
