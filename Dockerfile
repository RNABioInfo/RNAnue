############################################################
# Dockerfile to build RNAnue container
# Based on Debian
# v1.0.0
############################################################

# set the base image to debian
FROM ubuntu:23.04
# tag version
ARG VERSION=v1.0.0
# file author
LABEL authors="Christopher Adelmann and Richard A. Schaefer"

# update sources list
RUN apt-get update && apt-get -y upgrade
RUN apt-get install -y curl build-essential cmake git pkg-config
RUN apt-get install -y libboost-program-options-dev libbz2-dev zlib1g-dev libncurses5-dev liblzma-dev

# install htslib
WORKDIR /
RUN curl -L https://github.com/samtools/htslib/releases/download/1.20/htslib-1.20.tar.bz2 -o htslib-1.20.tar.bz2
RUN tar -xvf htslib-1.20.tar.bz2 && rm htslib-1.20.tar.bz2
WORKDIR /htslib-1.20
RUN ./configure
RUN make
RUN make install

# install segemehl
WORKDIR /
RUN curl -L http://legacy.bioinf.uni-leipzig.de/Software/segemehl/downloads/segemehl-0.3.4.tar.gz -o segemehl-0.3.4.tar.gz
RUN tar -xvf segemehl-0.3.4.tar.gz && rm segemehl-0.3.4.tar.gz
WORKDIR /segemehl-0.3.4
RUN export PKG_CONFIG_PATH=/usr/local/lib/pkgconfig
RUN make all
RUN cp segemehl.x /usr/local/bin
RUN echo 'alias segemehl="segemehl.x"' >> ~/.bashrc

# TODO Change to main when release is available
# retrieve RNAnue
WORKDIR /
RUN git clone -b develop --recurse-submodules https://github.com/ChristopherAdelmann/RNAnue.git
WORKDIR /RNAnue

# install RNAnue
WORKDIR /RNAnue
RUN cmake --preset release
RUN cmake --build --preset release --parallel 10
RUN echo 'alias RNAnue="/RNAnue/build/release/RNAnue"' >> ~/.bashrc
WORKDIR /
