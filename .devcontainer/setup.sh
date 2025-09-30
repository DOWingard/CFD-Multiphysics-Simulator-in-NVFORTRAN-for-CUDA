#!/bin/bash

export NVHPC_ROOT=/opt/nvidia/hpc_sdk
export NVHPC_VERSION=25.7
export NVHPC_HOME=${NVHPC_ROOT}/Linux_x86_64/${NVHPC_VERSION}

export CUDA_HOME=${NVHPC_HOME}/cuda
export CUDA_ROOT=${CUDA_HOME}
export CUDACXX=${CUDA_HOME}/bin/nvcc

export NVHPC_COMPILERS=${NVHPC_HOME}/compilers

export PATH=${NVHPC_COMPILERS}/bin:${CUDA_HOME}/bin:${PATH}
export LD_LIBRARY_PATH=${NVHPC_COMPILERS}/lib:${CUDA_HOME}/lib64:${LD_LIBRARY_PATH}


chmod +x scripts/*
chmod +x pyinterface/plot.py
chmod +x kernel/scripts/*


export PATH="$PATH:$(pwd)/scripts"
export PATH="$PATH:$(pwd)/pyinterface"
export PATH="$PATH:$(pwd)/kernel/scripts"

for dir in scripts pyinterface kernel/scripts; do
    grep -qxF "export PATH=\"\$PATH:$(pwd)/$dir\"" ~/.bashrc || \
        echo "export PATH=\"\$PATH:$(pwd)/$dir\"" >> ~/.bashrc
done

source ~/.bashrc

echo "Setup complete. Scripts and plot.py are executable and PATH is updated."


pip3 install numpy matplotlib jupyter ipython imageio scipy

WORKSPACE_DIR=$(pwd)
mkdir -p "${WORKSPACE_DIR}/build"

# if [ -f "${WORKSPACE_DIR}/test_cuda.cu" ]; then
#     /opt/nvidia/hpc_sdk/Linux_x86_64/25.7/cuda/bin/nvcc \
#         -allow-unsupported-compiler \
#         -arch=sm_75 \
#         -o build/test_cuda \
#         test_cuda.cu
# fi
