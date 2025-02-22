FROM ubuntu:24.04
ENV LANG=C.UTF-8
ENV LC_ALL=C.UTF-8
ADD image_development/bioscanflow_env.yml /bioscanflow_env.yml
RUN apt-get update && \
    apt-get install wget -y && \
    apt-get install unzip -y && \
    apt-get install libjpeg9 -y

ENV CONDA_DIR=/opt/conda
RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh && \
    /bin/bash ~/miniconda.sh -b -p /opt/conda

# Put conda in path so we can use conda activate
ENV PATH=$CONDA_DIR/bin:$PATH

RUN conda env create -f /bioscanflow_env.yml    
RUN eval "$(conda shell.bash hook)" && \
    conda activate bioscanflow_env  && \
    Rscript -e "library(devtools); devtools::install_github('tobiasgf/lulu')"

RUN wget https://github.com/mothur/mothur/releases/download/v1.48.2/Mothur.Ubuntu_22.zip -O /tmp/mothur.zip && \
    unzip /tmp/mothur.zip -d /usr/bin/ && rm /tmp/mothur.zip
ENV PATH="$PATH:/usr/bin/mothur"

ENTRYPOINT ["bash", "-c", "eval \"$(conda shell.bash hook)\" && conda activate bioscanflow_env && exec bash"]
