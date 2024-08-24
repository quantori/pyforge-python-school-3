from fastapi import FastAPI

from src.exceptions import BadRequestException, UnknownIdentifierException
from src.handler import bad_request_exception_handler, unknown_identifier_exception_handler
from src.users.router import router as user_router
from src.molecules.router import router as molecule_router

app = FastAPI()

# register the routers
app.include_router(user_router, prefix="/users")
app.include_router(molecule_router, prefix="/molecules")


app.add_exception_handler(BadRequestException, bad_request_exception_handler)
app.add_exception_handler(UnknownIdentifierException, unknown_identifier_exception_handler)
