from typing import Annotated

from fastapi import FastAPI, Depends, status, Body, Path, Request
from fastapi.encoders import jsonable_encoder
from starlette.responses import JSONResponse

from src.exceptions import BadRequestException, UnknownIdentifierException
from src.schemas import MoleculeRequest, MoleculeResponse
from src.dependencies import get_molecule_service, get_pagination_query_params
from src.schemas import PaginationQueryParams
from src.service import MoleculeService

app = FastAPI()


@app.exception_handler(BadRequestException)
def validation_exception_handler(request: Request, exc: BadRequestException):
    return JSONResponse(
        status_code=status.HTTP_400_BAD_REQUEST,
        content=jsonable_encoder({"detail": exc.message}),
    )


@app.exception_handler(UnknownIdentifierException)
def unknown_identifier_exception_handler(
    request: Request, exc: UnknownIdentifierException
):
    return JSONResponse(
        status_code=status.HTTP_404_NOT_FOUND,
        content=jsonable_encoder({"detail": exc.message}),
    )


@app.post(
    "/molecules/",
    status_code=201,
    responses={
        status.HTTP_201_CREATED: {"model": MoleculeResponse},
        status.HTTP_400_BAD_REQUEST: {
            "model": str,
            "description": "Probably due to Invalid SMILES string, or smiles uniqueness violation",
        },
    },
)
def add_molecule(
    molecule_request: Annotated[MoleculeRequest, Body(...)],
    service: Annotated[MoleculeService, Depends(get_molecule_service)],
) -> MoleculeResponse:
    return service.save(molecule_request)


@app.get(
    "/molecules/{molecule_id}",
    status_code=200,
    responses={
        status.HTTP_200_OK: {"model": MoleculeResponse},
        status.HTTP_404_NOT_FOUND: {
            "model": str,
            "description": "Molecule with the given ID not found",
        },
    },
)
def get_molecule(
    molecule_id: Annotated[
        int, Path(..., description="Unique identifier for the molecule")
    ],
    service: Annotated[MoleculeService, Depends(get_molecule_service)],
) -> MoleculeResponse:
    return service.find_by_id(molecule_id)


@app.get(
    "/molecules/",
    status_code=200,
    responses={
        status.HTTP_200_OK: {"model": list[MoleculeResponse]},
    },
)
def get_molecules(
    service: Annotated[MoleculeService, Depends(get_molecule_service)],
    pagination: Annotated[PaginationQueryParams, Depends(get_pagination_query_params)],
) -> list[MoleculeResponse]:
    return service.find_all(pagination.page, pagination.page_size)


@app.put(
    "/molecules/{molecule_id}",
    status_code=200,
    responses={
        status.HTTP_200_OK: {"model": MoleculeResponse},
        status.HTTP_404_NOT_FOUND: {
            "model": str,
            "description": "Molecule with the given ID not found",
        },
    },
)
def update_molecule(
    molecule_id: Annotated[
        int, Path(..., description="Unique identifier for the molecule")
    ],
    molecule_request: Annotated[MoleculeRequest, Body(...)],
    service: Annotated[MoleculeService, Depends(get_molecule_service)],
) -> MoleculeResponse:
    return service.update(molecule_id, molecule_request)


@app.delete(
    "/molecules/{molecule_id}",
    status_code=200,
    responses={
        status.HTTP_200_OK: {"description": "Molecule deleted successfully"},
        status.HTTP_404_NOT_FOUND: {
            "model": str,
            "description": "Molecule with the given ID not found",
        },
    },
)
def delete_molecule(
    molecule_id: Annotated[
        int, Path(..., description="Unique identifier for the molecule")
    ],
    service: Annotated[MoleculeService, Depends(get_molecule_service)],
):
    return service.delete(molecule_id)
