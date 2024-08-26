from fastapi import FastAPI
from src.users.router import router as user_router
from src.molecules.router import router as molecule_router
from src.handler import register_exception_handlers

app = FastAPI()

# register the routers
app.include_router(user_router, prefix="/users")
app.include_router(molecule_router, prefix="/molecules")
register_exception_handlers(app)
