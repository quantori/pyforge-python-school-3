from fastapi import FastAPI
import logging

# Create an instance of FastAPI
app = FastAPI()

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# Example endpoint using logging
@app.post("/add")
def add_molecule(molecule: dict):  # Assume a dictionary for simplicity
    logger.info("Adding a new molecule")
    # Simulate adding the molecule
    return {"message": "Molecule added", "molecule": molecule}

# Another example for logging
@app.get("/")
def read_root():
    logger.info("Root endpoint accessed")
    return {"message": "Hello World"}
