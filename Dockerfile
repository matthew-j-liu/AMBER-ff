# syntax=docker/dockerfile:1.6
#
# Ready-to-use container for the AMBER-ff project.
# Bundles AmberTools 25 + Armadillo (conda-forge), a C++20 toolchain, and the
# compiled `amber_ff` / `dihedral_scan` binaries.
#
# Build:  docker build -t amber-ff .
# Run:    docker run --rm -it -v "$PWD/input_molecules_copy:/data" amber-ff \
#             amber_ff /data/ethane.xyz

# ---------- stage 1: build ----------------------------------------------------
FROM condaforge/miniforge3:24.9.0-0 AS build

ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update && apt-get install -y --no-install-recommends \
        build-essential cmake git ca-certificates \
    && rm -rf /var/lib/apt/lists/*

# AmberTools 25 + Armadillo from conda-forge (prebuilt; no source compile).
RUN mamba install -n base -c conda-forge -y ambertools=24 armadillo openbabel h5py \
    && mamba clean -afy
ENV AMBERHOME=/opt/conda
ENV PATH=/opt/conda/bin:$PATH
# Make conda's Armadillo discoverable by CMake.
ENV CMAKE_PREFIX_PATH=/opt/conda

WORKDIR /src
COPY CMakeLists.txt ./
COPY src/ ./src/
RUN cmake -B build -S . -DCMAKE_BUILD_TYPE=Release \
    && cmake --build build -j"$(nproc)"

# ---------- stage 2: runtime --------------------------------------------------
FROM condaforge/miniforge3:24.9.0-0 AS runtime

RUN mamba install -n base -c conda-forge -y ambertools=24 armadillo openbabel h5py \
    && mamba clean -afy
ENV AMBERHOME=/opt/conda
ENV PATH=/opt/conda/bin:/usr/local/bin:$PATH

# Project binaries + data needed by the benchmark scripts.
COPY --from=build /src/build/amber_ff       /usr/local/bin/amber_ff
COPY --from=build /src/build/dihedral_scan  /usr/local/bin/dihedral_scan
WORKDIR /work
COPY params/              /work/params/
COPY input_molecules_copy/ /work/input_molecules_copy/
COPY benchmarks/          /work/benchmarks/
COPY run_AMBER_benchmarking.sh run_AMBER_custom.sh test.sh /work/

CMD ["bash"]
