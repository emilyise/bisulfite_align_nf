FROM continuumio/miniconda:4.7.12
LABEL authors="Robin F Chan" \
      description="Docker image containing all env reqs for robinfchan/bisulfite_align Nextflow pipeline"

# Install procps so that Nextflow can poll CPU usage
RUN apt-get update && apt-get install -y procps && apt-get clean -y 

# Install the conda environment
COPY environment.yml /
COPY trimRRBSdiversityAdaptCustomers.py /
RUN conda env create -f /environment.yml && conda clean -a
ENV PATH /opt/conda/envs/bisulfite_align/bin:$PATH
