import requests
import json
import redis
import time
from os import getenv


ENDPOINT = "http://localhost:8000"
REDIS_URL = getenv("REDIS_URL")

redis_client = redis.from_url(REDIS_URL)


def test_cache_molecule_by_id():
    """Test that molecule data is cached in Redis."""
    molecule_id = 4
    redis_key = f"molecule:{molecule_id}"

    redis_client.delete(redis_key)
    response = requests.get(ENDPOINT + f"/molecules/{molecule_id}")
    assert response.status_code == 200
    expected_data = response.json()

    cached_data = redis_client.get(redis_key)
    assert cached_data is not None
    assert json.loads(cached_data) == expected_data


def test_use_cached_molecule():
    """Test that a cached molecule is used on subsequent requests."""
    molecule_id = 4
    redis_key = f"molecule:{molecule_id}"

    cached_data = redis_client.get(redis_key)
    assert cached_data is not None

    response = requests.get(ENDPOINT + f"/molecules/{molecule_id}")
    assert response.status_code == 200
    assert response.json() == json.loads(cached_data)


def test_cache_expiration():
    """Test that the cache expires after the specified time."""
    molecule_id = 4
    redis_key = f"molecule:{molecule_id}"

    response = requests.get(ENDPOINT + f"/molecules/{molecule_id}")
    assert response.status_code == 200

    redis_client.expire(redis_key, 2)

    time.sleep(3)

    cached_data = redis_client.get(redis_key)
    assert cached_data is None


def test_cache_invalidation_on_update():
    """Test that the cache is invalidated when a molecule is updated."""
    molecule_id = 2
    redis_key = f"molecule:{molecule_id}"

    response = requests.get(ENDPOINT + f"/molecules/{molecule_id}")
    assert response.status_code == 200

    cached_data = redis_client.get(redis_key)
    assert cached_data is not None

    response = requests.put(
        ENDPOINT + f"/molecules/{molecule_id}",
        params={"name": "CNO"}
    )
    assert response.status_code == 200

    cached_data = redis_client.get(redis_key)
    assert cached_data is None


def test_redis_connection_failure(monkeypatch):
    """Test that the application handles Redis connection failure."""
    def mock_redis_client(*args, **kwargs):
        raise redis.exceptions.ConnectionError

    monkeypatch.setattr(redis.Redis, 'get', mock_redis_client)

    molecule_id = 4
    response = requests.get(ENDPOINT + f"/molecules/{molecule_id}")

    assert response.status_code == 200
    assert response.json() is not None


def upload_molecules_json(filename='molecules.json'):
    """Helper function to upload a molecules JSON file."""
    with open(filename, 'rb') as file:
        files = {'file': ('molecules.json', file, 'application/json')}
        response = requests.post(ENDPOINT + "/upload_file/", files=files)
    return response


def test_upload_file_invalid_json():
    """Test uploading an invalid JSON file."""
    files = {
        'file': (
            'molecules.json',
            '{"mol_id": 5, "name": "C1=CC=CC=C1"',
            'application/json'
        )
    }
    response = requests.post(ENDPOINT + "/upload_file/", files=files)
    assert response.status_code == 400
    assert response.json() == {"detail": "Invalid JSON file"}


def test_upload_file_success():
    """Test successful file upload."""
    response = upload_molecules_json()
    assert response.status_code == 201
    assert response.json() == {
        "message": "File uploaded and molecules parsed successfully",
        "num_molecules": 10
    }


def test_get_server():
    """Test to check if the server is up and running."""
    response = requests.get(ENDPOINT)
    assert response.status_code == 200
    assert "server_id" in response.json()


def test_get_molecule_by_id():
    response = requests.get(ENDPOINT + "/molecules/4")
    assert response.status_code == 200
    assert response.json() == {"id": 4, "name": "CNC"}


def test_update_molecule():
    response = requests.put(
        ENDPOINT + "/molecules/2",
        params={"name": "CNO"}
    )
    assert response.status_code == 200
    assert response.json() == {"message": "Molecule updated successfully"}


def test_delete_molecule():
    response = requests.delete(ENDPOINT + "/molecules/7")
    assert response.status_code == 200
    assert response.json() == {"message": "The molecule with id 7 is deleted!"}
