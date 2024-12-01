# 1. Using a base image with Python 3.12
FROM python:3.12

# 2. Create a working directory
WORKDIR /fastapi_app

# 3. Copy requirements.txt and install dependencies via pip, including RDKit
COPY requirements.txt /fastapi_app/
RUN pip install --upgrade pip && \
    pip install --no-cache-dir -r requirements.txt

# 4. Copy the application code into the container
COPY . /fastapi_app

# 5. Opening a port for the application
EXPOSE 8010

# 6. Launching FastAPI server with Uvicorn, available on port 8010, with automatic reboot
CMD ["uvicorn", "src.main:app", "--host", "0.0.0.0", "--port", "8010", "--reload"]
