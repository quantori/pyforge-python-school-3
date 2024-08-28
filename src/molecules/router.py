from typing import Annotated
from fastapi import Depends, status, Body, Path, Query, UploadFile, APIRouter
from src.molecules.schemas import MoleculeRequest, MoleculeResponse
from src.molecules.dependencies import get_pagination_query_params
from src.molecules.service import get_molecule_service
from src.molecules.schemas import PaginationQueryParams
from src.molecules.service import MoleculeService

router = APIRouter()


@router.post(
    "/",
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


@router.get(
    "/{molecule_id}",
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


@router.get(
    "/",
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


@router.put(
    "/{molecule_id}",
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


@router.delete(
    "/{molecule_id}",
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
) -> bool:
    return service.delete(molecule_id)


@router.get(
    "/search/substructures",
    responses={
        status.HTTP_200_OK: {"model": list[MoleculeResponse]},
        status.HTTP_400_BAD_REQUEST: {
            "model": str,
            "description": "Probably due to Invalid SMILES string",
        },
    },
)
def substructure_search(
    smiles: Annotated[
        str,
        Query(
            ...,
            description="Find substructures of the given SMILES string",
        ),
    ],
    service: Annotated[MoleculeService, Depends(get_molecule_service)],
) -> list[MoleculeResponse]:
    """
    Find all molecules that ARE SUBSTRUCTURES of the given smile, not vice vera.
    """
    return service.get_substructures(smiles)


@router.get(
    "/search/substructure_of",
    responses={
        status.HTTP_200_OK: {"model": list[MoleculeResponse]},
        status.HTTP_400_BAD_REQUEST: {
            "model": str,
            "description": "Probably due to Invalid SMILES string",
        },
    },
)
def substructure_search_of(
    smiles: Annotated[
        str,
        Query(
            ...,
            description="SMILES string that has to be substructure of the found molecules",
        ),
    ],
    service: Annotated[MoleculeService, Depends(get_molecule_service)],
) -> list[MoleculeResponse]:
    """
    Find all molecules that the given smile IS SUBSTRUCTURE OF, not vice vera.
    """
    return service.get_is_substructure_of(smiles)


@router.post("/upload/upload_molecules_csv", status_code=status.HTTP_201_CREATED)
def upload_molecules(
    file: UploadFile,
    service: Annotated[MoleculeService, Depends(get_molecule_service)],
):
    """
    Upload a CSV file containing molecules to the repository.

    The CSV file should have the following columns: smiles,name

    Lines that have incorrect format, missing smiles string or invalid smiles string are ignored.
    """

    # Uploaded CSV file is not stored on the server, only the molecules are extracted and stored in the memory.

    res = service.process_csv_file(file)
    return {"number_of_molecules_added": res}
