FROM nvidia/cuda:12.6.2-devel-ubuntu22.04

RUN apt-get update && apt-get install -y \
    build-essential \
    cmake \
    git \
    gdb \
    vim \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /workspace
