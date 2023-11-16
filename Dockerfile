FROM rocker/r-ver:4.3.1

RUN apt-get update && apt-get install -y -qq --no-install-recommends \
  libgsl0-dev libxml2-dev libboost-all-dev libssl-dev libhdf5-dev unzip curl \
  libudunits2-dev libgdal-dev libgeos-dev libproj-dev libgmp3-dev libreadline-dev \
  build-essential libglpk-dev xorg-dev libreadline-dev libc6-dev libclang-dev \
  libbz2-dev liblzma-dev libcurl4-openssl-dev libcairo2-dev libatlas-base-dev \
  libunwind-dev libpango1.0-dev tcl-dev tk-dev openjdk-17-jdk gfortran gnupg2 \
  lsb-release software-properties-common libpng-dev zlib1g-dev libicu-dev \
  libfftw3-dev libgfortran5 libgraphviz-dev libgtk-3-dev libgtkmm-3.0-dev \
  libxt-dev pandoc

RUN apt-get install -y -qq --no-install-recommends python3 python3-dev python3-pip

RUN pip3 install umap-learn leidenalg igraph

RUN Rscript -e "install.packages(c('BiocManager', 'devtools', 'randomcoloR', 'Seurat', 'readr', 'plyr'))"

RUN Rscript -e "devtools::install_github('iaaka/visutils')"

RUN Rscript -e "BiocManager::install('DropletUtils')"

RUN Rscript -e "install.packages('SeuratObject')"

RUN Rscript -e "devtools::install_github('mojaveazure/seurat-disk')"

ENV LD_LIBRARY_PATH=/usr/local/lib/R/lib:/usr/local/lib:/usr/lib/x86_64-linux-gnu:/usr/lib/jvm/java-17-openjdk-amd64/lib/server
