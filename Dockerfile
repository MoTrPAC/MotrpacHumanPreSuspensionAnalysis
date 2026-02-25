# Dockerfile for testing MotrpacHumanPreSuspensionAnalysis installation
FROM rocker/r-ver:4.4.0

LABEL maintainer="MoTrPAC"
LABEL description="Test installation of MotrpacHumanPreSuspensionAnalysis R package"

# System dependencies required by R packages
RUN apt-get update && apt-get install -y --no-install-recommends \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libfontconfig1-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    cmake \
    libgit2-dev \
    libglpk-dev \
    git \
    tcl-dev \
    tk-dev \
    libx11-dev \
    libxt-dev \
    xvfb \
    pandoc \
    && rm -rf /var/lib/apt/lists/*

# Install BiocManager first
RUN R -e "install.packages('BiocManager', repos='https://cloud.r-project.org')"

# Install Bioconductor dependencies (from Remotes)
# Mfuzz depends on e1071 and Biobase; install its deps first
RUN R -e "install.packages('e1071', repos='https://cloud.r-project.org')"
RUN R -e "BiocManager::install(c('ComplexHeatmap', 'Biobase', 'Mfuzz', 'TMSig', 'variancePartition'), ask=FALSE, update=TRUE, force=TRUE)"

# Verify all Bioconductor packages installed
RUN R -e " \
  pkgs <- c('ComplexHeatmap', 'Biobase', 'Mfuzz', 'TMSig', 'variancePartition'); \
  missing <- pkgs[!sapply(pkgs, requireNamespace, quietly=TRUE)]; \
  if (length(missing)) stop('Failed to install: ', paste(missing, collapse=', ')) \
"

# Install CRAN Imports
RUN R -e "install.packages(c( \
    'data.table', 'dplyr', 'magrittr', 'tibble', 'tidyr', \
    'ggplot2', 'ggpubr', 'latex2exp', 'forcats', 'scales', 'stringr' \
  ), repos='https://cloud.r-project.org')"

# Install CRAN Suggests
RUN R -e "install.packages(c( \
    'here', 'knitr', 'pheatmap', 'rmarkdown', 'readr', \
    'roxygen2', 'table.glue', 'testthat' \
  ), repos='https://cloud.r-project.org')"

# Install GitHub Remotes
RUN R -e "install.packages('remotes', repos='https://cloud.r-project.org')"
RUN R -e "remotes::install_github('MoTrPAC/MotrpacBicQC', upgrade='never')"

# Copy the package source
COPY . /tmp/package
WORKDIR /tmp/package

# Install all package dependencies declared in DESCRIPTION (Imports/Suggests/Remotes)
RUN R -e "install.packages('pak', repos='https://cloud.r-project.org')"
RUN R -e "pak::local_install_dev_deps(dependencies = TRUE)"

# Verify a known required dependency is available before package install
RUN R -e "if (!requireNamespace('ggpubr', quietly = TRUE)) stop('ggpubr is not installed')"

# Build and install the package with R CMD
RUN R CMD build --no-build-vignettes . && \
    R CMD INSTALL MotrpacHumanPreSuspensionAnalysis_*.tar.gz

# Verify: load the package and print session info
RUN R -e " \
  library(MotrpacHumanPreSuspensionAnalysis); \
  cat('\n=== Package loaded successfully ===\n\n'); \
  print(sessionInfo()) \
"

CMD ["R", "-e", "library(MotrpacHumanPreSuspensionAnalysis); cat('Package is ready.\\n')"]
