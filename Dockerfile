# Multi-stage Dockerfile for Amplicon Analysis Pipeline
# Uses micromamba for fast Conda environment creation

FROM mambaorg/micromamba:1.5.6 as base

# Set working directory
WORKDIR /app

# Copy environment file first (better layer caching)
COPY --chown=$MAMBA_USER:$MAMBA_USER analysis/environment.yml /tmp/environment.yml

# Create the conda environment
RUN micromamba create -y -n amplicon_pipeline -f /tmp/environment.yml && \
    micromamba clean --all --yes

# Activate environment by default
ARG MAMBA_DOCKERFILE_ACTIVATE=1
ENV ENV_NAME=amplicon_pipeline

# Copy the rest of the project
COPY --chown=$MAMBA_USER:$MAMBA_USER . /app

# Set the default shell to use the activated environment
SHELL ["micromamba", "run", "-n", "amplicon_pipeline", "/bin/bash", "-c"]

# Verify installations
RUN nextflow -version && \
    bwa 2>&1 | head -1 && \
    samtools --version | head -1 && \
    gatk --version | head -1

# Default command
ENTRYPOINT ["micromamba", "run", "-n", "amplicon_pipeline"]
CMD ["nextflow", "run", "/app/analysis/main.nf", "--help"]
