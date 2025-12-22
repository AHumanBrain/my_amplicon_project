#!/bin/bash
set -e

echo "--- Starting WSL Setup ---"

# 1. Install Miniconda if not present
if ! command -v conda &> /dev/null; then
    echo "Miniconda not found. Installing..."
    mkdir -p ~/miniconda3
    wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda3/miniconda.sh
    bash ~/miniconda3/miniconda.sh -b -u -p ~/miniconda3
    rm -rf ~/miniconda3/miniconda.sh
    
    # Initialize conda for bash
    ~/miniconda3/bin/conda init bash
    eval "$(/home/$USER/miniconda3/bin/conda shell.bash hook)"
    echo "Miniconda installed successfully."
    
    # Remove defaults to avoid ToS issues
    ~/miniconda3/bin/conda config --remove channels defaults || true
    ~/miniconda3/bin/conda config --add channels conda-forge
else
    echo "Miniconda already installed."
    eval "$(conda shell.bash hook)"
    
    # Ensure defaults is removed for existing installs too
    conda config --remove channels defaults || true
    conda config --add channels conda-forge
fi

# 2. Java is required for Nextflow.
# We will install it via Conda (openjdk) in the environment, so system-level install is suspended to avoid sudo prompts.
echo "Skipping system-level Java install (will rely on Conda openjdk)..."


# 3. Create/Update Conda Environment
ENV_FILE=analysis/environment.yml
ENV_NAME=amplicon_pipeline

echo "Creating/Updating Conda environment '$ENV_NAME'..."

# Use direct conda create with --override-channels to STRICTLY avoid defaults/ToS issues
# We explicitly list packages here to mirror environment.yml
PACKAGES="nextflow fastqc cutadapt bwa samtools picard gatk4 deeptools multiqc openjdk=17"

if conda env list | grep -q "$ENV_NAME"; then
    echo "Environment exists. updating..."
    conda install -n $ENV_NAME -c conda-forge -c bioconda --override-channels $PACKAGES -y
else
    echo "Creating environment..."
    conda create -n $ENV_NAME -c conda-forge -c bioconda --override-channels $PACKAGES -y
fi

echo "--- Setup Complete ---"
echo "To run the analysis, use: wsl bash -c 'source ~/.bashrc && conda activate amplicon_pipeline && cd analysis && bash run_analysis.sh'"
