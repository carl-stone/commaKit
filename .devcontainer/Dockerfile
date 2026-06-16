FROM bioconductor/bioconductor_docker:RELEASE_3_22

ENV DEBIAN_FRONTEND=noninteractive
ENV RENV_CONFIG_SANDBOX_ENABLED=FALSE
ENV RENV_CONFIG_SYNCHRONIZED_CHECK=FALSE

USER root

RUN apt-get update && apt-get install -y --no-install-recommends \
    bash-completion \
    bubblewrap \
    build-essential \
    ca-certificates \
    cmake \
    curl \
    gdb \
    git \
    git-lfs \
    less \
    libbz2-dev \
    libcairo2-dev \
    libcurl4-openssl-dev \
    libfontconfig1-dev \
    libfreetype6-dev \
    libfribidi-dev \
    libgit2-dev \
    libglpk-dev \
    libharfbuzz-dev \
    libhdf5-dev \
    libicu-dev \
    libjpeg-dev \
    liblzma-dev \
    libpng-dev \
    libssl-dev \
    libtiff-dev \
    libxml2-dev \
    libxt-dev \
    make \
    pandoc \
    python3-venv \
    qpdf \
    wget \
    xz-utils \
    zlib1g-dev \
  && rm -rf /var/lib/apt/lists/*

RUN python3 -m venv /opt/radian \
  && /opt/radian/bin/pip install --no-cache-dir --upgrade pip \
  && /opt/radian/bin/pip install --no-cache-dir radian \
  && ln -s /opt/radian/bin/radian /usr/local/bin/radian

RUN Rscript --vanilla -e "install.packages('https://cran.r-project.org/src/contrib/Archive/renv/renv_1.1.8.tar.gz', repos = NULL, type = 'source')" \
  && Rscript --vanilla -e "install.packages(c('devtools', 'roxygen2', 'testthat', 'pkgdown', 'languageserver', 'httpgd', 'lintr', 'styler'), repos = 'https://cloud.r-project.org')"

RUN apt-get update && apt-get install -y --no-install-recommends ripgrep \
  && rm -rf /var/lib/apt/lists/*

WORKDIR /workspaces/commaKit
