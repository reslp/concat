FROM ubuntu:18.04
MAINTAINER <philipp.resl@uni-graz.at>

RUN apt-get update && \
	apt-get install -y wget git python2.7 && \
	ln -s /usr/bin/python2.7 /usr/bin/python

WORKDIR /software



RUN git clone https://github.com/reslp/concat.git
ENV PATH="/software/concat:$PATH"

CMD ["concat.py"]


