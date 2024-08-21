# Use an official Miniconda image as the base image
FROM continuumio/miniconda3:latest

# Set the working directory
WORKDIR /app

# Install required system packages for mysqlclient
RUN apt-get update && apt-get install -y \
    pkg-config \
    libmysqlclient-dev \
    && rm -rf /var/lib/apt/lists/*

# Install RDKit via conda
RUN conda install -c conda-forge rdkit

# Copy the requirements file into the container
COPY requirements.txt .

# Install Python dependencies
RUN pip3 install -r requirements.txt

# Copy the current directory contents into the container at /app
COPY . .

# Command to run the application
CMD ["uvicorn", "main:app", "--host", "0.0.0.0", "--port", "8000"]
