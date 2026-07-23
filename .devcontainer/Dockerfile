FROM bioconductor/bioconductor_docker:RELEASE_3_22@sha256:e34177896c1579f66d48b756e4229b58237499183a089da084af72d7aba1d40a

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

RUN Rscript --vanilla -e "install.packages('https://cran.r-project.org/src/contrib/Archive/renv/renv_1.1.8.tar.gz', repos = NULL, type = 'source')"

RUN apt-get update && apt-get install -y --no-install-recommends ripgrep \
  && rm -rf /var/lib/apt/lists/*

ARG AIR_VERSION=0.10.0
ARG TARGETARCH
RUN set -eux; \
  case "${TARGETARCH}" in \
    amd64) \
      air_arch="x86_64"; \
      air_sha256="abfb6db0f5f45701a17b9c49d5323eecde6586cba909bcd7392d938716821493" \
      ;; \
    arm64) \
      air_arch="aarch64"; \
      air_sha256="4873b28b29e242061d4604d5c6d004d0902cc7a0b59fa79e59b521bfc44b13e8" \
      ;; \
    *) \
      echo "Unsupported architecture: ${TARGETARCH}" >&2; \
      exit 1 \
      ;; \
  esac; \
  air_archive="/tmp/air-${AIR_VERSION}.tar.gz"; \
  curl --proto '=https' --tlsv1.2 -fsSL \
    "https://github.com/posit-dev/air/releases/download/${AIR_VERSION}/air-${air_arch}-unknown-linux-gnu.tar.gz" \
    -o "${air_archive}"; \
  echo "${air_sha256}  ${air_archive}" | sha256sum -c -; \
  tar -xzf "${air_archive}" \
    --strip-components=1 \
    -C /usr/local/bin \
    "air-${air_arch}-unknown-linux-gnu/air"; \
  chmod 0755 /usr/local/bin/air; \
  rm "${air_archive}"; \
  air --version

RUN id rstudio \
  && mkdir -p /workspaces/commaKit /home/rstudio/.codex \
  && chown -R rstudio:rstudio /workspaces/commaKit /home/rstudio/.codex

USER rstudio

WORKDIR /workspaces/commaKit
