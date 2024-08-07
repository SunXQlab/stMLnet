# hash:sha256:7677ead8cdf7806b1a4148916cbdc12c67940c3eee8d3add2e58602faa84621c
FROM registry.codeocean.com/codeocean/r-studio:1.2.5019-r4.0.3-ubuntu18.04

ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        libgmp-dev=2:6.1.2+dfsg-2ubuntu0.1 \
        libgsl0-dev \
        libmagick++-dev=8:6.9.7.4+dfsg-16ubuntu6.15 \
    && rm -rf /var/lib/apt/lists/*

RUN Rscript -e 'remotes::install_version("Matrix", "1.5-4.1")' \
    && Rscript -e 'remotes::install_version("R.utils", "2.12.2")' \
    && Rscript -e 'remotes::install_version("ROCR", "1.0-11")' \
    && Rscript -e 'remotes::install_version("Rfast2", "0.1.4")' \
    && Rscript -e 'remotes::install_version("ape", "5.7-1")' \
    && Rscript -e 'remotes::install_version("caret", "6.0-94")' \
    && Rscript -e 'remotes::install_version("doParallel", "1.0.17")' \
    && Rscript -e 'remotes::install_version("doSNOW", "1.0.20")' \
    && Rscript -e 'remotes::install_version("dplyr", "1.1.2")' \
    && Rscript -e 'remotes::install_version("foreach", "1.5.2")' \
    && Rscript -e 'remotes::install_version("igraph", "1.4.3")' \
    && Rscript -e 'remotes::install_version("plotrix", "3.8-2")' \
    && Rscript -e 'remotes::install_version("ranger", "0.12.1")' \
    && Rscript -e 'remotes::install_version("rvcheck", "0.1.8")' \
    && Rscript -e 'remotes::install_version("scatterpie", "0.1.8")'

RUN Rscript -e 'options(warn=2); install.packages("BiocManager")'
RUN Rscript -e 'options(warn=2); BiocManager::install(c( \
        "ggalluvial", \
        "ggsci", \
        "org.Hs.eg.db" \
    ))' # Original versions: 0.12.5 3.0.0 3.12.0

COPY postInstall /
RUN /postInstall
