from fastapi import FastAPI, HTTPException, status, UploadFile, File
from rdkit import Chem
import json
import logging
from os import getenv
import uvicorn
from molecules.schema import MoleculeAdd, MoleculeResponse
from molecules.dao import MoleculeDAO
from typing import List, Dict
from dotenv import load_dotenv
import redis
import traceback

logging.basicConfig(
    filename='app.log',
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

app = FastAPI()

load_dotenv(".env")

DB_URL = getenv("DB_URL")
REDIS_URL = getenv("REDIS_URL")

redis_client = redis.from_url(REDIS_URL)

try:
    redis_client.ping()
    logger.info("Successfully connected to Redis")
except redis.ConnectionError:
    logger.error("Unable to connect to Redis")


if DB_URL is None:
    logger.error("Database URL is not in env")
    raise ValueError("Database URL is not in env")

if REDIS_URL is None:
    logger.error("Redis URL is not in env")
    raise ValueError("Redis URL is not in env")


@app.get("/")
def get_server():
    server_id = getenv("SERVER_ID", "1")
    logger.info(
        f"Server {server_id} received a request to the root endpoint"
    )
    return {"server_id": server_id, "message": "Welcome to Molecules app!"}


@app.post("/molecules", status_code=status.HTTP_201_CREATED)
async def add_molecule(molecule: MoleculeAdd):
    try:
        logger.info(f"Adding molecule: {molecule.name}")
        existing_molecule = await MoleculeDAO.find_full_data(molecule.id)
        if existing_molecule:
            logger.warning(f"Molecule with ID {molecule.id} already exists")
            raise HTTPException(
                status_code=400,
                detail="Molecule with this ID already exists"
            )
        if not Chem.MolFromSmiles(molecule.name):
            logger.warning(f"Invalid SMILES: {molecule.name}")
            raise HTTPException(
                status_code=400,
                detail=f"Invalid SMILES: {molecule.name}"
            )
        await MoleculeDAO.add_mol(**molecule.model_dump())
        return {"message": "Molecule added successfully"}
    except Exception as e:
        logger.error(f"Error occurred while adding molecule: {str(e)}")
        raise HTTPException(status_code=500, detail="Internal Server Error")


@app.get(
        "/molecules",
        tags=["Molecules"],
        response_model=List[MoleculeResponse]
        )
async def retrieve_molecules(limit: int = 100) -> List[MoleculeResponse]:
    logger.info("Retrieving {limit} molecules")
    try:
        molecule_iterator = MoleculeDAO.find_all_molecules_iterator(limit)
        response = [
            MoleculeResponse(**molecule)
            async for molecule in molecule_iterator
        ]
        return response
    except Exception as e:
        logger.error(f"Error retrieving molecules: {e}")
        raise HTTPException(status_code=500, detail="Internal Server Error")


@app.get("/molecules/{mol_id}", response_model=MoleculeResponse)
async def get_molecule(mol_id: int):
    logging.info(f"Received request for molecule ID: {mol_id}")
    mol_data = await MoleculeDAO.find_full_data(mol_id)
    if mol_data is None:
        logging.warning(f"Molecule with ID {mol_id} not found")
        raise HTTPException(status_code=404, detail="Molecule not found")
    return mol_data


@app.put(
        "/molecules/{mol_id}",
        tags=["Molecules"],
        response_description="Update molecule by ID"
    )
async def update_molecule(mol_id: int, name: str):
    try:
        logger.info(f"Updating molecule with ID {mol_id} to new name {name}")
        await MoleculeDAO.update_mol(mol_id, name)
        return {"message": "Molecule updated successfully"}
    except ValueError as e:
        logger.warning(f"Update failed: {str(e)}")
        raise HTTPException(status_code=404, detail=str(e))
    except Exception as e:
        logger.error(f"Internal server error during molecule update: {str(e)}")
        raise HTTPException(status_code=500, detail="Internal Server Error")


@app.delete(
        "/molecules/{mol_id}",
        tags=["Molecules"],
        response_description="Delete molecule by ID"
    )
async def delete_mol(mol_id: int) -> dict:
    logger.info(f"Deleting molecule with ID {mol_id}")
    molecule = await MoleculeDAO.delete_mol_by_id(mol_id=mol_id)
    if molecule:
        logger.info(f"Molecule with ID {mol_id} deleted")
        return {"message": f"The molecule with id {mol_id} is deleted!"}
    else:
        logger.warning(f"Molecule with ID {mol_id} not found for deletion")
        raise HTTPException(status_code=404, detail="Molecule not found")


def get_cached_result(key: str):
    try:
        result = redis_client.get(key)
        if result:
            decoded_result = json.loads(result)
            logger.info(f"Cache hit for key: {key}")
            return decoded_result
        else:
            logger.info(f"Cache miss for key: {key}")
    except redis.RedisError as e:
        logger.error(f"Redis error during cache retrieval: {e}")
    except json.JSONDecodeError as e:
        logger.error(f"JSON decode error: {e}")
    return None


def set_cache(key: str, value: dict, expiration: int = 60):
    try:
        json_value = json.dumps(value)
        logger.info(f"Caching data for key: {key} with value: {json_value}")
        redis_client.setex(key, expiration, json_value)
        logger.info(f"Cache set for key: {key}")
    except redis.RedisError as e:
        logger.error(f"Redis error during cache setting: {e}")
    except json.JSONDecodeError as e:
        logger.error(f"JSON encode error: {e}")


@app.get("/substructure_search", tags=["Molecules"], response_model=List[Dict])
async def substructure_search(
    substructure_name: str,
    limit: int = 100
) -> List[Dict]:
    logger.info(f"Received request for substructure_search with substructure_name={substructure_name} and limit={limit}")

    if not substructure_name:
        logger.error("Substructure SMILES string cannot be empty")
        raise HTTPException(status_code=400, detail="Substructure SMILES string cannot be empty")

    substructure_mol = Chem.MolFromSmiles(substructure_name)
    if substructure_mol is None:
        logger.error("Invalid SMILES string")
        raise HTTPException(status_code=400, detail="Invalid SMILES string")
    
    try:
        key = f"substructure_search:{substructure_name}:{limit}"
        cached_result = get_cached_result(key)
        if cached_result:
            logger.info(f"Returning cached result for substructure: {substructure_name}")
            return {"source": "cache", "data": cached_result}

        logger.info(f"Searching for molecules with substructure: {substructure_name}")
        iterator = MoleculeDAO.find_by_substructure_iterator(substructure_name, limit)
        matches = [match async for match in iterator]

        if not matches:
            logger.warning(f"No molecules found for substructure: {substructure_name}")
            raise HTTPException(status_code=404, detail="No molecules found matching the substructure")

        logger.info(f"Setting cache for key: {key} with {len(matches)} matches")
        set_cache(key, matches, expiration=300)
        logger.info("Substructure search completed successfully")
        return matches

    except HTTPException as e:
        logger.warning(f"HTTP error during substructure search: {e.detail}")
        raise e
    except ValueError as e:
        logger.warning(f"Bad request for substructure search: {e}")
        raise HTTPException(status_code=400, detail="Bad request")
    except Exception as e:
        logger.error(f"An error occurred: {str(e)}", exc_info=True)
        raise HTTPException(status_code=500, detail="Internal Server Error")


@app.post(
        "/upload_file/",
        status_code=status.HTTP_201_CREATED,
        tags=["File Upload"],
        response_description="File uploaded and molecules parsed successfully"
)
async def upload_file(file: UploadFile = File(...)):
    logger.info(f"Uploading and processing file: {file.filename}")
    try:
        content = await file.read()
        content = content.decode("utf-8")
        molecules = json.loads(content)
        if not isinstance(molecules, list):
            logger.warning(
                "Invalid file format: Expected a list of molecules"
                )
            raise HTTPException(
                status_code=400,
                detail="Invalid file format: expected a list of molecules"
            )
    except (UnicodeDecodeError, json.JSONDecodeError) as e:
        logger.warning(f"Failed to process file {file.filename}: {str(e)}")
        raise HTTPException(status_code=400, detail="Invalid JSON file")
    except ValueError as e:
        logger.warning(f"Failed to process file {file.filename}: {str(e)}")
        raise HTTPException(status_code=400, detail="Invalid JSON file")
    added_count = 0
    for molecule in molecules:
        try:
            mol_id = molecule.get("mol_id")
            name = molecule.get("name")
            logger.debug(f"Processing molecule: ID={mol_id}, Name={name}")
            if not mol_id or not name:
                logger.warning(f"Missing ID or name in file: {molecule}")

            existing_molecule = await MoleculeDAO.find_full_data(mol_id)
            if existing_molecule:
                logger.warning(
                    f"Molecule with ID {mol_id} already exists, skipping."
                )
                continue

            mol = Chem.MolFromSmiles(name)
            if not mol:
                logger.warning(f"Invalid SMILES: {name}, skipping.")
                continue

            await MoleculeDAO.add_mol(id=mol_id, name=name)
            added_count += 1

        except Exception as e:
            logger.error(f"Error processing molecule {molecule}: {str(e)}")
            continue
    logger.info(
        f"File processing complete. Number of molecules added: {added_count}"
        )
    return {
        "message": "File uploaded and molecules parsed successfully",
        "num_molecules": added_count
    }


if __name__ == "__main__":
    port = int(getenv("PORT", 8000))
    logger.info(f"Starting server on port {port}")
    uvicorn.run(app, host="0.0.0.0", port=port)
