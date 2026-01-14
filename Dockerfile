FROM ubuntu:24.04
ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update && apt-get install -y curl bash openjdk-17-jre-headless python3 python3-pip python3-venv git build-essential fasttree && rm -rf /var/lib/apt/lists/*
RUN python3 -m venv /opt/venv
ENV PATH="/opt/venv/bin:${PATH}"
RUN curl -s https://get.nextflow.io | bash && mv nextflow /usr/local/bin/nextflow
RUN /opt/venv/bin/pip install --no-cache-dir biopython pandas networkx pyvis matplotlib seaborn numpy requests nextstrain-augur
COPY . /opt/pipeline
WORKDIR /opt/pipeline
CMD ["bash"]
