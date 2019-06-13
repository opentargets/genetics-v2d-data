FROM continuumio/miniconda:4.6.14

COPY ./ /v2d
WORKDIR /v2d
RUN conda env create -n v2d_data --file environment.yaml
RUN echo "source activate v2d_data" > ~/.bashrc
ENV PATH /opt/conda/envs/v2d_data/bin:$PATH
