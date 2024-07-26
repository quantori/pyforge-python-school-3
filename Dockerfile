FROM python:3.9-slim
WORKDIR /app/src
COPY ./requirements.txt /app
RUN pip install --no-cache-dir -r /app/requirements.txt
COPY ./src /app/src
EXPOSE 8000
ENTRYPOINT [ "uvicorn" ]
CMD [ "main:app", "--host", "0.0.0.0", "--port", "8000"]