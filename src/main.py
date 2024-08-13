from fastapi import FastAPI, Request
from src.api.routes.molecule_routes import router as molecule_router
import src.exceptions as exceptions
import src.exception_handlers as handlers

app = FastAPI()
app.include_router(molecule_router, prefix="/molecules", tags=["molecules"])


app.add_exception_handler(exceptions.InvalidSmilesException, handlers.invalid_smiles_exception_handler)
app.add_exception_handler(exceptions.UnknownIdentifierException, handlers.unknown_identifier_exception_handler)


# @app.post("/upload_molecules_csv", status_code=status.HTTP_201_CREATED)
# async def upload_molecules(file: UploadFile):
#     """
#     Upload a CSV file containing molecules to the repository.
#
#     Uploaded CSV file is not stored on the server, only the molecules are extracted and stored in the memory.
#
#     The CSV file should have the following columns: smiles,name,description
#
#     Lines that have incorrect format, missing smiles string or invalid smiles string are ignored.
#
#     :param file:
#
#     :return: A dictionary containing the key: number_of_molecules_added
#
#     :raises HTTPException: If the file is not a CSV file, or does not have the required columns.
#     """
#
#     if file.content_type != "text/csv":
#         raise HTTPException(status_code=status.HTTP_400_BAD_REQUEST, detail="File must be a CSV file")
#
#     contents = await file.read()
#
#     csv_data = contents.decode("utf-8")
#     # if the first line is not smiles,name,description raise an exception
#     # this is a very naive check, because more correct way would be to check column names, because
#     # the order of the columns might be different, or there might be more columns
#     # TODO: change this check to a better one
#     if not csv_data.startswith("smiles,name,description"):
#         raise HTTPException(status_code=status.HTTP_400_BAD_REQUEST, detail="CSV file should have the following "
#                                                                             "columns: smiles,name,description")
#
#     number_of_molecules_added = 0
#
#     csv_reader = csv.DictReader(csv_data.splitlines())
#     for row in csv_reader:
#         try:
#             validate_smiles(row["smiles"])
#             molecule = Molecule(smiles=row["smiles"], molecule_name=row["name"], description=row["description"])
#             molecules_repository.add(molecule)
#             number_of_molecules_added += 1
#         except exception.InvalidSmilesException:
#             # log somewhere. just printing this time, I have not looked into logging in FastAPI yet.
#             print(f"Invalid SMILES string: {row['smiles']}")
#
#     return {"number_of_molecules_added": number_of_molecules_added}
