# swglh/carp:1.0.0
# docker build --network=host -f Dockerfile_carp -t swglh/carp:1.0.0 .
# DNAnexus file-id: 

# author: Suzy Hocking
# date: 11/12/2025

FROM ubuntu:22.04

WORKDIR /usr

# set up environment
ENV DEBIAN_FRONTEND=noninteractive

# update and install required packages
RUN apt-get update && \
    apt-get install -y \
    sudo git nano python3 python3-dev python3-pip default-jre-headless bedtools \
    && apt-get clean && rm -rf /var/lib/apt/lists/*

# initialise pip
RUN python3 -m pip install --upgrade pip cython numpy Cmake wheel dxpy

# install requirements
COPY requirements.txt /usr/carp/requirements.txt
RUN python3 -m pip install -r /usr/carp/requirements.txt

# copy scripts
COPY Dockerfile_carp /usr/
COPY README.md /usr/carp/
COPY src/allele_fraction.py /usr/carp/
COPY src/savvycnv_dosage.py /usr/carp/
COPY src/generate_plots.py /usr/carp/