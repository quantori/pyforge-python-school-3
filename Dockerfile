FROM continuumio/miniconda3
RUN conda install conda-forge::rdkit


WORKDIR /app

COPY src/ ./app/
COPY tests/ ./app/
COPY requirements.txt ./app/

RUN pip install --no-cache-dir -r requirements.txt


RUN pytest tests/


ENV PORT=8000

EXPOSE 8000

CMD ["uvicorn", "main:app", "--host", "0.0.0.0", "--port", "8000"]
LABEL authors="Nika"