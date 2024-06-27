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
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Install CRAN packages
RUN Rscript -e 'install.packages(c("markdown","rmarkdown","shiny","HGNChelper","BiocManager"))'

# Install bioconductor packages
RUN Rscript -e 'BiocManager::install(c("IlluminaHumanMethylation450kanno.ilmn12.hg19","IlluminaHumanMethylationEPICanno.ilm10b4.hg19" ,"mitch"))'

# get a clone of the code
RUN git clone https://github.com/markziemann/gmea_app.git
ENV DIRPATH /gmea_app
WORKDIR $DIRPATH
RUN cp -r gmea_app/* /srv/shiny-server/

# Expose the application port
USER shiny
EXPOSE 3838

CMD ["/usr/bin/shiny-server"]
