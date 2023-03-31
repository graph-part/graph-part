FROM continuumio/miniconda3

WORKDIR /home/graph-part
#TODO do we need all these copies?
COPY graph_part graph_part
COPY setup.py .
COPY README.md .
COPY webserver/requirements.txt .
COPY benchmarking/data/protein/NetGPI/netgpi_dataset.fasta sample.fasta
#RUN conda install -c anaconda python=3.8
RUN pip install --no-cache-dir .
RUN pip install -r requirements.txt # server has additional requirements
RUN conda install -c bioconda emboss
RUN conda install -c conda-forge -c bioconda mmseqs2 

COPY webserver/make_output.py .
COPY webserver/run.sh .
RUN chmod +x run.sh