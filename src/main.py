from fastapi import FastAPI

app = FastAPI()

@app.get("/m")
def a():
    pass