FROM ghcr.io/web-mev/base-sctk-docker:sha-9210ec50158783a7a98f4a1a0acbab204aae2e76

# Install the optparse package, which allows better commandline arg parsing
RUN R -e "install.packages('optparse',dependencies=TRUE, repos='http://cran.rstudio.com/')"

RUN mkdir -p /opt/software
ADD decontx.R /opt/software/

ENTRYPOINT ["/bin/bash"]
