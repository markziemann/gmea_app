FROM bioconductor/bioconductor_docker:RELEASE_3_19

# Update apt-get
RUN apt-get update \
        && apt-get upgrade -y \
        && apt-get install -y nano git  libncurses-dev \
        ## Install the python package magic wormhole to send files
        && pip install magic-wormhole           \
        ## Remove packages in '/var/cache/' and 'var/lib'
        ## to remove side-effects of apt-get update
        && apt-get clean \
        && rm -rf /var/lib/apt/lists/*

# Install CRAN packages
RUN Rscript -e 'install.packages(c("markdown","rmarkdown","shiny","HGNChelper"))'

# Install bioconductor packages
RUN Rscript -e 'BiocManager::install(c("IlluminaHumanMethylation450kanno.ilmn12.hg19","IlluminaHumanMethylationEPICanno.ilm10b4.hg19" ,"mitch"))'

# get a clone of the code
RUN git clone https://github.com/markziemann/gmea_app.git

# Set the container working directory
ENV DIRPATH /gmea_app
WORKDIR $DIRPATH

# Expose the application port
EXPOSE 3838

# Run the R Shiny app
CMD Rscript /gmea_app/gmea_app/app.R
