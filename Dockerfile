FROM reslp/mamba:1.4.1
MAINTAINER <philipp.resl@uni-graz.at>

RUN apt-get --allow-releaseinfo-change update && \
	apt-get install -y wget git 

RUN conda config --append channels bioconda && \
	conda config --append channels conda-forge && \
	conda config --append channels anaconda && \
	mamba install -c bioconda biopython=1.77 && \
	mamba install -c bioconda mafft=7.464 && \
	conda clean -a -y

WORKDIR /software

RUN git clone https://github.com/reslp/concat.git
ENV PATH="/software/concat:$PATH"

WORKDIR /data 

CMD ["concat.py"]


