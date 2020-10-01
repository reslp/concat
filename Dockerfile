FROM continuumio/miniconda3:4.7.10
MAINTAINER <philipp.resl@uni-graz.at>

RUN apt-get update && \
	apt-get install -y wget git 

RUN conda config --append channels bioconda && \
	conda config --append channels conda-forge && \
	conda config --append channels anaconda && \
	conda install -c bioconda biopython=1.78 && \
	conda install -c bioconda mafft=7.464 && \
	conda clean -a -y

WORKDIR /software

RUN git clone https://github.com/reslp/concat.git
ENV PATH="/software/concat:$PATH"

CMD ["concat.py"]


