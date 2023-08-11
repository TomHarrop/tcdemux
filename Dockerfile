FROM python:3.10.12

LABEL SOFTWARE_NAME tcdemux
LABEL MAINTAINER "Tom Harrop"
LABEL version=0.0.1

ENV DEBIAN_FRONTEND=noninteractive
ENV LC_ALL=C
ENV PATH=/usr/share/bbmap:${PATH}

RUN     apt-get clean && \
        rm -r /var/lib/apt/lists/*

RUN     apt-get update && apt-get upgrade -y --fix-missing

RUN     apt-get install -y  --no-install-recommends \
            bbmap \
            pigz

COPY    . /source

RUN     /usr/local/bin/python3.10 \
            -m pip install /source

ENTRYPOINT ["/usr/local/bin/tcdemux"]
