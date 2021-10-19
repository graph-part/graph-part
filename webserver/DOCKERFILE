FROM python:3.8-slim

WORKDIR /home/graph-part
#TODO do we need all these copies?
COPY graph-part graph-part
COPY setup.py .
COPY README.md .
COPY webserver/requirements.txt .
RUN pip3 install --no-cache-dir .


#COPY webserver/generate_output.py .
COPY webserver/run.sh .
RUN chmod +x run.sh