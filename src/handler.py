from fastapi.encoders import jsonable_encoder
from starlette import status
from starlette.responses import JSONResponse
from fastapi import Request
from src.exceptions import BadRequestException, UnknownIdentifierException


def bad_request_exception_handler(request: Request, exc: BadRequestException):
    return JSONResponse(
        status_code=status.HTTP_400_BAD_REQUEST,
        content=jsonable_encoder({"detail": exc.message}),
    )


def unknown_identifier_exception_handler(
    request: Request, exc: UnknownIdentifierException
):
    return JSONResponse(
        status_code=status.HTTP_404_NOT_FOUND,
        content=jsonable_encoder({"detail": exc.message}),
    )
