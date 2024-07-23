from fastapi import FastAPI
from src.routers import molecules

app = FastAPI(title="Python Summer School 2024: SMILES")
app.include_router(molecules.router, prefix="/molecules", tags=["molecules"])
