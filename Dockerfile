# Use Miniconda base image
FROM continuumio/miniconda3

# Set the working directory
WORKDIR /app

# Install FastAPI and Uvicorn
RUN pip install fastapi uvicorn

# Install RDKit from conda-forge
RUN conda install -c conda-forge rdkit

# Copy the application code to the container
COPY ./src /app

# Command to run the FastAPI application
CMD ["uvicorn", "main:app", "--host", "0.0.0.0", "--port", "8000"]
