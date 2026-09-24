FROM sneumann/spcca-base:4.5

LABEL maintainer="Steffen Neumann <sneumann@ipb-halle.de>"
LABEL Description="spCCA Web Application"

WORKDIR /srv/shiny-server

# Remove all Demo Apps
RUN rm -rf *

## Now install real MetFamily
ADD . /tmp/spCCA
RUN R CMD INSTALL /tmp/spCCA

WORKDIR /srv/shiny-server
RUN rm -rf *
ADD inst/spCCA /srv/shiny-server/

# Make readable by shiny user
RUN chmod -R 755 /srv/shiny-server/*



