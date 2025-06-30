#!/bin/bash

# Foolproof conda environment installation script
# This script installs all required packages in the correct order to avoid conflicts

set -e  # Exit on any error

ENV_NAME="graphbp_complete"
PYTHON_VERSION="3.10"

echo "============================================"
echo "Creating conda environment: $ENV_NAME"
echo "============================================"

# Remove existing environment if it exists
conda env remove -n $ENV_NAME -y 2>/dev/null || true

# Create new environment with Python
echo "Step 1: Creating base environment with Python $PYTHON_VERSION"
conda create -n $ENV_NAME python=$PYTHON_VERSION -y

# Activate environment
echo "Step 2: Activating environment"
source $(conda info --base)/etc/profile.d/conda.sh
conda activate $ENV_NAME

# Install mamba first for faster package resolution
echo "Step 3: Installing mamba"
conda install -c conda-forge mamba -y

# Install conda-only packages first (these must be installed via conda)
echo "Step 4: Installing conda-only packages"
mamba install -c conda-forge -c bioconda \
    rdkit \
    openbabel \
    snakemake \
    -y

# Install core scientific packages
echo "Step 5: Installing core scientific packages"
mamba install -c conda-forge \
    numpy \
    pandas \
    matplotlib \
    scikit-learn \
    networkx \
    -y

# Install PyTorch (CPU version first, then we'll handle CUDA if needed)
echo "Step 6: Installing PyTorch"
# Check if CUDA is available
if command -v nvidia-smi &> /dev/null; then
    echo "CUDA detected, installing PyTorch with CUDA support"
    mamba install pytorch torchvision torchaudio pytorch-cuda=11.8 -c pytorch -c nvidia -y
else
    echo "No CUDA detected, installing CPU-only PyTorch"
    mamba install pytorch torchvision torchaudio cpuonly -c pytorch -y
fi

# Install PyTorch Geometric and related packages
echo "Step 7: Installing PyTorch Geometric ecosystem"
# Get PyTorch and CUDA versions for compatibility
TORCH_VERSION=$(python -c "import torch; print(torch.__version__)")
CUDA_VERSION=$(python -c "import torch; print(torch.version.cuda if torch.cuda.is_available() else 'cpu')")
echo "Detected PyTorch version: $TORCH_VERSION"
echo "Detected CUDA version: $CUDA_VERSION"

# Install PyTorch Geometric first
pip install torch-geometric

# Install torch-scatter, torch-sparse, torch-spline-conv with pre-built wheels
echo "Installing PyTorch Geometric extensions..."

if [[ "$CUDA_VERSION" != "cpu" ]]; then
    echo "Attempting to install CUDA-enabled extensions..."
    # Try multiple common CUDA versions in order of preference
    CUDA_VERSIONS=("cu118" "cu117" "cu116" "cu113" "cu111")
    
    for cuda_ver in "${CUDA_VERSIONS[@]}"; do
        echo "Trying with $cuda_ver..."
        if pip install torch-scatter torch-sparse torch-spline-conv -f https://data.pyg.org/whl/torch-$(echo $TORCH_VERSION | cut -d'+' -f1)+${cuda_ver}.html; then
            echo "✓ Successfully installed with $cuda_ver"
            break
        else
            echo "Failed with $cuda_ver, trying next..."
        fi
    done
else
    echo "Installing CPU-only extensions..."
    pip install torch-scatter torch-sparse torch-spline-conv -f https://data.pyg.org/whl/torch-$(echo $TORCH_VERSION | cut -d'+' -f1)+cpu.html
fi

# Fallback: try conda-forge if pip fails
if ! python -c "import torch_scatter" 2>/dev/null; then
    echo "Pip installation failed, trying conda-forge as fallback..."
    mamba install -c conda-forge -c pyg pytorch-scatter pytorch-sparse pytorch-spline-conv -y
fi

# Install remaining packages via pip
echo "Step 8: Installing remaining packages"
pip install biopython

# Verify installation
echo "Step 9: Verifying installation"
python -c "
import sys
failed_imports = []

packages = [
    ('rdkit', 'rdkit'),
    ('openbabel', 'openbabel'),
    ('torch', 'torch'),
    ('torch_geometric', 'torch_geometric'),
    ('torch_scatter', 'torch_scatter'),
    ('torch_sparse', 'torch_sparse'),
    ('networkx', 'networkx'),
    ('numpy', 'numpy'),
    ('pandas', 'pandas'),
    ('matplotlib', 'matplotlib'),
    ('sklearn', 'scikit-learn'),
    ('Bio', 'biopython')
]

for module, package_name in packages:
    try:
        __import__(module)
        print(f'✓ {package_name}')
    except ImportError as e:
        print(f'✗ {package_name}: {e}')
        failed_imports.append(package_name)

if failed_imports:
    print(f'\nFailed to import: {failed_imports}')
    print('Please check the installation logs above for errors.')
    sys.exit(1)
else:
    print('\n✓ All packages imported successfully!')
    
    # Print versions for key packages
    import rdkit
    import torch
    import torch_geometric
    print(f'RDKit version: {rdkit.__version__}')
    print(f'PyTorch version: {torch.__version__}')
    print(f'PyTorch Geometric version: {torch_geometric.__version__}')
    print(f'CUDA available: {torch.cuda.is_available()}')
    if torch.cuda.is_available():
        print(f'CUDA version: {torch.version.cuda}')
"

echo "============================================"
echo "Environment '$ENV_NAME' created successfully!"
echo "============================================"
echo ""
echo "To activate the environment, use:"
echo "conda activate $ENV_NAME"
echo ""
echo "To deactivate, use:"
echo "conda deactivate"
