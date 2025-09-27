FROM ubuntu:22.04

RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    gfortran \
    cmake \
    git \
    wget \
    curl \
    python3 \
    python3-pip \
    nvidia-cuda-toolkit \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /workspace

SHELL ["/bin/bash", "-c"]
