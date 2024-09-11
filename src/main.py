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
from tasks import substructure_search_task
from celery.result import AsyncResult
from celery_worker import celery

logging.basicConfig(
    filename='app.log',
    level=logging.DEBUG,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

app = FastAPI()

load_dotenv(".env")

DB_URL = getenv("DB_URL")
REDIS_URL = getenv("REDIS_URL", "redis://redis_cache:6379/0")


if DB_URL is None:
    logger.error("Database URL is not in env")
    raise ValueError("Database URL is not in env")

if REDIS_URL is None:
    logger.error("Redis URL is not in env")
    raise ValueError("Redis URL is not in env")


redis_client = redis.StrictRedis(
    host="redis_cache",
    port=6379,
    db=0
)

try:
    redis_client.ping()
    logger.info("Redis is available")
except redis.ConnectionError as e:
    logger.error(f"Redis connection failed: {e}")
    raise HTTPException(status_code=500, detail="Redis connection failed")


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


def set_cache(key: str, value: dict, expiration: int = 3600):
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
    logger.info(
        f"Request: substructure_name={substructure_name}, limit={limit}"
    )

    if not substructure_name:
        logger.error("Substructure SMILES string cannot be empty")
        raise HTTPException(
            status_code=400,
            detail="Substructure SMILES string cannot be empty"
        )

    substructure_mol = Chem.MolFromSmiles(substructure_name)
    if substructure_mol is None:
        logger.error("Invalid SMILES string")
        raise HTTPException(status_code=400, detail="Invalid SMILES string")

    key = f"substructure_search:{substructure_name}:{limit}"
    cached_result = get_cached_result(key)
    if cached_result:
        logger.info(f"Cache hit: Returning cached result for key: {key}")
        return {"source": "cache", "data": cached_result}

    logger.info(f"Initiating substructure search for {substructure_name}")

    task = substructure_search_task.delay(substructure_name, limit)
    try:
        result = task.get(timeout=300)
    except Exception as e:
        logger.error(f"Task execution failed: {e}")
        raise HTTPException(
            status_code=500,
            detail="Internal server error during task execution"
        )

    if not result:
        logger.warning(f"No molecules found for {substructure_name}")
        raise HTTPException(
            status_code=404,
            detail="No molecules found matching the substructure"
        )

    set_cache(key, result, expiration=3600)
    logger.info(f"Cache set: Stored result for key {key}")
    return {"source": "database", "data": result}


@app.post("/tasks/substructure_search")
async def create_substructure_search_task(
    substructure_name: str,
    limit: int = 100
):
    if not substructure_name:
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail="Substructure SMILES string cannot be empty"
        )
    task = substructure_search_task.delay(substructure_name, limit)
    return {"task_id": task.id, "status": task.status}


@app.get("/tasks/{task_id}")
async def get_task_result(task_id: str):
    task_result = AsyncResult(task_id, app=celery)
    if task_result.state == 'PENDING':
        return {"task_id": task_id, "status": "Task is still processing"}
    elif task_result.state == 'SUCCESS':
        return {
            "task_id": task_id, "status": "Task completed",
            "result": task_result.result
        }
    else:
        return {"task_id": task_id, "status": task_result.state}


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
    cache_key = "uploaded_molecules"
    for molecule in molecules:
        try:
            mol_id = molecule.get("mol_id")
            name = molecule.get("name")
            logger.debug(f"Processing molecule: ID={mol_id}, Name={name}")
            if not mol_id or not name:
                logger.warning(f"Missing ID or name in file: {molecule}")
                continue

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

    cached_data = {
        "num_molecules": added_count
    }
    set_cache(cache_key, cached_data, expiration=300)

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
    uvicorn.run(app, host="0.0.0.0", port=port, debug=True)
