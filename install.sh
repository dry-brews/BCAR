#!/bin/bash

# Function to check if command exists
command_exists() {
    command -v "$1" >/dev/null 2>&1
}

# Check if conda is installed
if ! command_exists conda; then
    echo "Conda is not installed. Please install Conda first."
    echo "Visit: https://docs.conda.io/projects/conda/en/latest/user-guide/install/"
    exit 1
fi

# Create/update conda environment
echo "Setting up conda environment..."
if conda env list | grep -q "bcar2-env"; then
    echo "Environment exists..."
    #conda env update -f environment.yml
else
    conda env create -f environment.yml
fi

# Activate environment and compile
echo "Activating environment and compiling..."
eval "$(conda shell.bash hook)"
conda activate bcar2-env

# Run make
make clean
make

echo "Installation complete!"
