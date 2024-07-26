import uvicorn
from fastapi import FastAPI

from src.routers import molecules

app = FastAPI(title="Python Summer School 2024: SMILES")
app.include_router(molecules.router, prefix="/molecules", tags=["molecules"])


if __name__ == "__main__":
    uvicorn.run(app, host="127.0.0.1", port=8000)
