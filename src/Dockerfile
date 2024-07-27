FROM continuumio/miniconda3

RUN conda create --name myenv python=3.11

RUN echo "source activate myenv" > ~/.bashrc

ENV PATH /opt/conda/envs/myenv/bin:$PATH

RUN conda install --name myenv -c conda-forge rdkit

WORKDIR /app

COPY . /app

COPY requirements.txt .

RUN pip install --no-cache-dir -r requirements.txt

EXPOSE 8000

CMD ["python", "app.py"]
