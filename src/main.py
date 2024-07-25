from fastapi import FastAPI
from src.routers import molecules

app = FastAPI()

app.include_router(molecules.router)