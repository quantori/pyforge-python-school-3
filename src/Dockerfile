FROM continuumio/miniconda3

ENV PATH="/opt/conda/envs/myenv/bin:$PATH"

WORKDIR /app

COPY requirements.txt /app/

RUN conda create --name myenv python=3.9 -y && \
    conda run -n myenv pip install --no-cache-dir -r requirements.txt && \
    conda run -n myenv conda install -y -c conda-forge rdkit && \
    conda run -n myenv pip install alembic

COPY . /app

ENV PYTHONPATH=/app

EXPOSE 8000

CMD ["conda", "run", "-n", "myenv", "uvicorn", "main:app", "--host", "0.0.0.0", "--port", "8000"]