from fastapi.testclient import TestClient
from fastapi import FastAPI

from src.routers.molecules import router


app = FastAPI()
app.include_router(router)
client = TestClient(app)


def test_read_main():
    response = client.get("/molecules")
    assert response.status_code == 200
