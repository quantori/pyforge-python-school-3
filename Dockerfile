# I was using conda environment before but i found it too inconvenient to use with docker
# So i switched to using pip and requirements.txt
# This base image example is from fastapi official documentation
FROM tiangolo/uvicorn-gunicorn-fastapi:python3.10

# Install rdkit library with pip, This is a big library and takes a while to install,
# so  I think it's better to install it before copying the rest of the dependencies,
RUN pip install rdkit

# Copy requirements.txt to the container at /app
COPY requirements.txt /app/

# Install the dependencies
RUN pip install -r requirements.txt

# Copy the content of the local src directory to the container at /app
COPY . /app/

# Set the working directory in the container
WORKDIR /app/

# make migrations
#RUN ["alembic", "upgrade", "head"]

EXPOSE 8000

ENTRYPOINT ["fastapi", "run", "src/main.py"]
