#!/bin/bash

# Quick fix script for PyTorch Geometric installation issues
# Run this if you already have PyTorch installed but torch-scatter/sparse/spline-conv failed

echo "============================================"
echo "Fixing PyTorch Geometric Extensions"
echo "============================================"

# Activate your environment (replace with your actual environment name)
ENV_NAME="graphbp_complete"  # Change this to your environment name
source $(conda info --base)/etc/profile.d/conda.sh
conda activate $ENV_NAME

# Get PyTorch and CUDA versions
echo "Detecting PyTorch and CUDA versions..."
TORCH_VERSION=$(python -c "import torch; print(torch.__version__)")
CUDA_AVAILABLE=$(python -c "import torch; print(torch.cuda.is_available())")
echo "PyTorch version: $TORCH_VERSION"
echo "CUDA available: $CUDA_AVAILABLE"

# Remove any failed installations
echo "Cleaning up failed installations..."
pip uninstall torch-scatter torch-sparse torch-spline-conv -y 2>/dev/null || true

# Install PyTorch Geometric first if not already installed
pip install torch-geometric

# Install extensions with pre-built wheels
echo "Installing PyTorch Geometric extensions..."

if [[ "$CUDA_AVAILABLE" == "True" ]]; then
    echo "Installing CUDA-enabled extensions..."
    
    # Try different CUDA versions
    CUDA_VERSIONS=("cu118" "cu117" "cu116" "cu113" "cu111")
    SUCCESS=false
    
    for cuda_ver in "${CUDA_VERSIONS[@]}"; do
        echo "Trying $cuda_ver..."
        if pip install torch-scatter torch-sparse torch-spline-conv -f https://data.pyg.org/whl/torch-$(echo $TORCH_VERSION | cut -d'+' -f1)+${cuda_ver}.html --no-cache-dir; then
            echo "✓ Success with $cuda_ver!"
            SUCCESS=true
            break
        else
            echo "Failed with $cuda_ver"
        fi
    done
    
    if [[ "$SUCCESS" == "false" ]]; then
        echo "All CUDA versions failed, trying conda-forge..."
        mamba install -c conda-forge -c pyg pytorch-scatter pytorch-sparse pytorch-spline-conv -y
    fi
else
    echo "Installing CPU-only extensions..."
    pip install torch-scatter torch-sparse torch-spline-conv -f https://data.pyg.org/whl/torch-$(echo $TORCH_VERSION | cut -d'+' -f1)+cpu.html --no-cache-dir
fi

# Test installation
echo "Testing installation..."
python -c "
try:
    import torch_scatter
    import torch_sparse
    print('✓ PyTorch Geometric extensions installed successfully!')
except ImportError as e:
    print(f'✗ Installation failed: {e}')
    exit(1)
"

echo "✓ Fix completed successfully!"
