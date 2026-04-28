FROM rocker/r-ver:4.4.2

RUN apt-get update && apt-get install -y --no-install-recommends \
    make git libcurl4-openssl-dev libssl-dev libxml2-dev \
    libcairo2-dev libxt-dev libfontconfig1-dev \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /bs2
RUN mkdir -p scripts
COPY DESCRIPTION README.md /bs2/
COPY scripts/00_install.R /bs2/scripts/00_install.R
RUN Rscript scripts/00_install.R

COPY . /bs2

ENV OMP_NUM_THREADS=1 \
    OPENBLAS_NUM_THREADS=1 \
    MKL_NUM_THREADS=1 \
    BLAS_NUM_THREADS=1 \
    VECLIB_MAXIMUM_THREADS=1

CMD ["Rscript", "scripts/00_run_all_quick.R"]
