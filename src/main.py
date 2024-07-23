from fastapi import FastAPI
from src.routers import moleculs

app = FastAPI(title="Python Summer School 2024: SMILES")
app.include_router(moleculs.router, prefix="/molecules", tags=["molecules"])
