# NOT IN USE

FROM ubuntu:20.04

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y \
    openjdk-11-jre-headless \
    wget \
    unzip \
    perl \
    python3 \
    python3-pip \
    && apt-get clean

# FastQC
RUN wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.9.zip \
    && unzip fastqc_v0.11.9.zip \
    && chmod +x FastQC/fastqc \
    && mv FastQC /opt/fastqc

# BTools
RUN wget https://downloads.sourceforge.net/project/bbmap/BBMap_38.94.tar.gz \
    && tar -xzf BBMap_38.94.tar.gz \
    && mv bbmap /opt/bbmap



WORKDIR /data

COPY read_qc.sh /opt/read_qc.sh
RUN chmod +x /opt/read_qc.sh

ENTRYPOINT ["/opt/read_qc.sh"]
