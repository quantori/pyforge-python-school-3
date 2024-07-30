FROM python:3.12-slim
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt
WORKDIR /app
COPY /src/main.py .
CMD uvicorn --host 0.0.0.0 --port 8000 --reload main:app