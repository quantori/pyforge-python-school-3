

# I was using conda environment before but i found it too inconvenient to use with docker
# So i switched to using pip and requirements.txt
# This base image example is from fastapi official documentation
FROM tiangolo/uvicorn-gunicorn-fastapi:python3.10
# add maintaner information and other important information

LABEL maintainer="gaiozi tabatadzegaga@gmail.com"
LABEL description="Be careful, migrations are applied automatically:\n\
 Currently best way to do this is to ovverride entrypoint.\n\
 for example: \n\
 entrypoint: ['/bin/sh', '-c', 'alembic revision --autogenerate && alembic upgrade head && fastapi run src/main.py']"

# Install rdkit library with pip, This is a big library and takes a while to install,
# so  I think it's better to install it before copying the rest of the dependencies,

RUN pip install rdkit

WORKDIR /app/

# Copy requirements.txt to the container at /app
COPY requirements.txt .

# Install the dependencies
RUN pip install -r requirements.txt

# Copy the content of the local src directory to the container at /app
COPY . .

EXPOSE 8000

ENTRYPOINT ["fastapi", "run", "src/main.py"]
