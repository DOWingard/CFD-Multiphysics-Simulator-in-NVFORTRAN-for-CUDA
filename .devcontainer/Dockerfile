FROM ubuntu:22.04

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y \
    build-essential \
    cmake \
    git \
    wget \
    curl \
    vim \
    nano \
    gdb \
    valgrind \
    pkg-config \
    software-properties-common \
    python3 \
    python3-pip \
    sudo \
    ca-certificates \
    gnupg \
    lsb-release \
    && rm -rf /var/lib/apt/lists/*

ARG USERNAME=vscode
ARG USER_UID=1000
ARG USER_GID=$USER_UID

RUN groupadd --gid $USER_GID $USERNAME \
    && useradd --uid $USER_UID --gid $USER_GID -m $USERNAME -s /bin/bash \
    && echo $USERNAME ALL=\(root\) NOPASSWD:ALL > /etc/sudoers.d/$USERNAME \
    && chmod 0440 /etc/sudoers.d/$USERNAME

RUN mkdir -p /opt/nvidia/hpc_sdk

ENV NVHPC_ROOT=/opt/nvidia/hpc_sdk
ENV CUDA_HOME=/opt/nvidia/hpc_sdk/Linux_x86_64/25.7/cuda
ENV PATH=/opt/nvidia/hpc_sdk/Linux_x86_64/25.7/compilers/bin:/opt/nvidia/hpc_sdk/Linux_x86_64/25.7/cuda/bin:${PATH}
ENV LD_LIBRARY_PATH=/opt/nvidia/hpc_sdk/Linux_x86_64/25.7/compilers/lib:/opt/nvidia/hpc_sdk/Linux_x86_64/25.7/cuda/lib64:${LD_LIBRARY_PATH}
ENV MANPATH=/opt/nvidia/hpc_sdk/Linux_x86_64/25.7/compilers/man:${MANPATH}

RUN pip3 install numpy matplotlib jupyter ipython
RUN chmox +x scripts/*
RUN export PATH="$PATH:$(pwd)/scripts"
RUN source ~/.bashrc

USER $USERNAME

WORKDIR /workspace