#!/bin/bash

# --- Quartet2RIdeogram Dependency Installation Script ---

echo "Starting dependency installation..."

# Check if conda is available
if ! command -v conda &> /dev/null
then
    echo "Error: Conda is not installed or not in your PATH."
    echo "Please install Miniconda or Anaconda first: https://docs.conda.io/en/latest/miniconda.html"
    exit 1
fi

# Check if a conda environment is active
if [ -z "$CONDA_DEFAULT_ENV" ]; then
    echo "Warning: No Conda environment is active."
    echo "It is highly recommended to run this in a dedicated environment."
    echo "You can create one with: conda create -n q2r python=3.9 r-base=4.2 -y"
    read -p "Do you want to continue anyway? (y/n) " -n 1 -r
    echo
    if [[ ! $REPLY =~ ^[Yy]$ ]]; then
        exit 1
    fi
fi

# --- Install Python Packages ---
echo ""
echo "--- Installing Python packages using pip... ---"
pip install pandas

if [ $? -ne 0 ]; then
    echo "Error: Failed to install Python packages."
    exit 1
fi
echo "Python packages installed successfully."


# --- Install R Packages ---
echo ""
echo "--- Installing R packages... ---"
R -e '
options(repos = c(CRAN = "https://cloud.r-project.org/"))
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("RIdeogram", quietly = TRUE)) BiocManager::install("RIdeogram")
if (!requireNamespace("stringr", quietly = TRUE)) install.packages("stringr")
if (!requireNamespace("getopt", quietly = TRUE)) install.packages("getopt")
'

if [ $? -ne 0 ]; then
    echo "Error: Failed to install R packages."
    exit 1
fi
echo "R packages installed successfully."

echo ""
echo "--- Installation complete! ---"
echo "You can now run the main pipeline using quartet2rideogram.sh"
