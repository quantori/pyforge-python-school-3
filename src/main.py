from fastapi import FastAPI
from src.users.router import router as user_router

app = FastAPI()

# register the routers
app.include_router(user_router, prefix="/users")
