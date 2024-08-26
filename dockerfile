FROM python:3.10-slim

WORKDIR /app

COPY . /app

RUN pip install --no-cache-dir fastapi uvicorn

EXPOSE 8000

ENV PYTHONUNBUFFERED=1

CMD ["uvicorn", "homework_10:app", "--host", "0.0.0.0", "--port", "8000"]
