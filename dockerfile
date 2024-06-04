FROM continuumio/miniconda3

# Python 3.8, OpenJDK
RUN conda create -n myenv python=3.8 openjdk=11 -y

RUN echo "conda activate myenv" >> ~/.bashrc
SHELL ["conda", "run", "-n", "myenv", "/bin/bash", "-c"]

RUN conda install -n myenv -c conda-forge numpy pandas biopython -y

WORKDIR /workspace

COPY . /workspace

CMD ["bash"]