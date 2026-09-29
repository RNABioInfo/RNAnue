# Build exactly this checkout, including its pinned submodules.
FROM ubuntu:24.04
ARG VERSION=v1.0.0
ARG BUILD_JOBS=2
LABEL authors="Christopher Adelmann and Richard A. Schaefer"
RUN apt-get update && apt-get install -y --no-install-recommends \
    gcc-14 g++-14 make cmake ninja-build git pkg-config xxd patchelf \
    autoconf automake libtool python3 gnuplot-nox libboost-program-options-dev \
    libbz2-dev zlib1g-dev liblzma-dev libhts-dev libtbb-dev libpng-dev \
    libncurses-dev ca-certificates && rm -rf /var/lib/apt/lists/*
WORKDIR /src/RNAnue
# Keep Git metadata: CMake verifies the STAR gitlink against the supported pin.
COPY . .
RUN cmake --preset release -G Ninja \
    -DCMAKE_C_COMPILER=gcc-14 -DCMAKE_CXX_COMPILER=g++-14 \
    -DCMAKE_INSTALL_PREFIX=/opt/rnanue -DRNANUE_BUILD_STAR=ON \
    && cmake --build --preset release --parallel ${BUILD_JOBS} \
    && cmake --install build/release \
    && /opt/rnanue/libexec/rnanue/STAR --version
ENV PATH="/opt/rnanue/bin:${PATH}"
WORKDIR /data
CMD ["RNAnue", "--help"]
