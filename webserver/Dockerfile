FROM continuumio/miniconda3

WORKDIR /home/graph-part
#TODO do we need all these copies?
COPY graph_part graph_part
COPY setup.py .
COPY README.md .
COPY webserver/requirements.txt .
#RUN conda install -c anaconda python=3.8
RUN pip install --no-cache-dir .
RUN conda install -c bioconda emboss

COPY webserver/make_output.py .
COPY webserver/run.sh .
RUN chmod +x run.sh