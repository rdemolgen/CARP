# swglh/carp:1.1.0
# docker build --network=host -f Dockerfile -t swglh/carp:1.1.0 .
# DNAnexus file-id: project-G729Kkj4fq4Q9X1BPy0807bK:file-J7P5Xj84fq4jx7fP6Pg5Jvz2

# author: Suzy Hocking
# date: 11/12/2025

FROM ubuntu:22.04

WORKDIR /usr/carp

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

# copy Dockerfile
COPY Dockerfile /usr/carp/

# copy scripts
COPY README.md /usr/carp/
COPY src/carp.py /usr/carp/src/
COPY src/baf.py /usr/carp/src/
COPY src/dosage.py /usr/carp/src/
COPY src/plots.py /usr/carp/src/
COPY src/utility.py /usr/carp/src/

# copy resource files
COPY resources/carp /usr/carp/resources/
COPY resources/hg38_cytoBand.txt /usr/carp/resources/
COPY resources/web_ClinGen_region_curation_list_GRCh38_20250425.tsv /usr/carp/resources/

# copy unit test files
COPY tests/test_carp.py /usr/carp/tests/
COPY tests/test_baf.py /usr/carp/tests/
COPY tests/test_dosage.py /usr/carp/tests/
COPY tests/test_plots.py /usr/carp/tests/
COPY tests/test_utility.py /usr/carp/tests/
