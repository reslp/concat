FROM reslp/mamba:0.25.0
MAINTAINER <philipp.resl@uni-graz.at>

RUN apt-get --allow-releaseinfo-change update && \
	apt-get install -y wget git 

RUN mamba install -c bioconda biopython=1.77 mafft=7.464 && conda clean -a -y

WORKDIR /software

RUN git clone https://github.com/reslp/concat.git
ENV PATH="/software/concat:$PATH"

WORKDIR /data 

CMD ["concat.py"]


