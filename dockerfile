FROM ubuntu:22.04
ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y \
    wget curl git build-essential unzip ca-certificates \
    openjdk-11-jre-headless python3 python3-pip \
    zlib1g-dev libbz2-dev liblzma-dev libcurl4-openssl-dev libssl-dev && \
    rm -rf /var/lib/apt/lists/*

RUN curl -s https://get.nextflow.io | bash && mv nextflow /usr/local/bin/

RUN wget -q https://github.com/OpenGene/fastp/releases/download/v0.23.2/fastp && \
    chmod +x fastp && mv fastp /usr/local/bin/

RUN wget -q https://github.com/bwa-mem2/bwa-mem2/releases/download/v2.2.1/bwa-mem2 && \
    chmod +x bwa-mem2 && mv bwa-mem2 /usr/local/bin/

RUN wget -q https://github.com/samtools/samtools/releases/download/1.17/samtools-1.17.tar.bz2 && \
    tar -xjf samtools-1.17.tar.bz2 && cd samtools-1.17 && \
    ./configure --prefix=/usr && make -j4 && make install && cd .. && rm -rf samtools-1.17*

RUN wget -q https://github.com/samtools/bcftools/releases/download/1.17/bcftools-1.17.tar.bz2 && \
    tar -xjf bcftools-1.17.tar.bz2 && cd bcftools-1.17 && \
    ./configure --prefix=/usr && make -j4 && make install && cd .. && rm -rf bcftools-1.17*

RUN mkdir -p /opt/snpeff && \
    wget -q -O /opt/snpeff/snpEff_latest_core.zip https://snpeff.blob.core.windows.net/versions/snpEff_latest_core.zip && \
    cd /opt/snpeff && unzip snpEff_latest_core.zip && rm snpEff_latest_core.zip

VOLUME ["/opt/annovar"]

WORKDIR /data
ENTRYPOINT ["/bin/bash"]
