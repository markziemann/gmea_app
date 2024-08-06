FROM rocker/shiny:4.4.1

# Install system requirements for index.R as needed
RUN apt-get update && apt-get install -y \
    --no-install-recommends \
    git-core \
    libssl-dev \
    libcurl4-gnutls-dev \
    curl \
    libncurses-dev \
    libsodium-dev \
    libxml2-dev \
    libicu-dev \
    liblzma-dev \
    libbz2-dev \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Install CRAN packages
RUN Rscript -e 'install.packages(c("pkgload","markdown","rmarkdown","shiny","BiocManager","RhpcBLASctl","shinybusy"))'

# Install bioconductor packages
RUN Rscript -e 'BiocManager::install(c("mitch"))'

# get a clone of the code
RUN git clone https://github.com/markziemann/gmea_app.git
ENV DIRPATH /gmea_app
WORKDIR $DIRPATH
COPY report.Rmd intro.md app.R 450K.rds EPIC.rds /srv/shiny-server/
COPY favicon.ico /www/favicon.ico

# Expose the application port
USER shiny
EXPOSE 3838

CMD ["/usr/bin/shiny-server"]
