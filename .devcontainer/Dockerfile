FROM ubuntu:22.04

ENV DEBIAN_FRONTEND=noninteractive
ENV PATH=/home/VOID/.local/bin:$PATH

RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        adduser \
        bash \
        bash-completion \
        build-essential \
        cmake \
        curl \
        g++ \
        gcc \
        gfortran \
        git \
        libblas3 \
        liblapack3 \
        make \
        nano \
        python3 \
        python3-pip \
        sudo \
        wget \
        zlib1g \
    && rm -rf /var/lib/apt/lists/*

RUN git config --global user.name "Derek Wingard" && \
    git config --global user.email "derekwingard@gmail.com"

RUN useradd -ms /bin/bash VOID && \
    echo "VOID ALL=(ALL) NOPASSWD:ALL" >> /etc/sudoers && \
    mkdir -p /home/VOID/app && \
    chown -R VOID:VOID /home/VOID/app

USER VOID
WORKDIR /home/VOID/app

RUN [ -d scripts ] && chmod +x scripts/* || echo "No scripts folder yet"
ENV PATH=scripts:$PATH



RUN echo "source /etc/bash_completion" >> /home/VOID/.bashrc && \
    echo "alias lr='ls -lart'" >> /home/VOID/.bashrc && \
    echo "alias cd..='cd ..'" >> /home/VOID/.bashrc

RUN pip install --user \
        fortls  \
        json5  \
        packaging 
