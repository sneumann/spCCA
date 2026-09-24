# Notes for developpers

## Docker images
The image `sneumann/spcca-base` contains all dependencies
for the spCCA web application.

### Building the container(s)
Build via `docker build -t sneumann/spcca-base -f Dockerfile-base . `or with the correct tagging:
```
docker build -t sneumann/spcca-base:`grep ^FROM Dockerfile-base | cut -d: -f 2 | tr -d " "` -f Dockerfile-base .
```

The image `sneumann/spcca` is built on top and contains
the actual spCCA package and web application.

Build via `docker build -t sneumann/spcca . ` or with the correct tagging:
```
docker build -t sneumann/spcca:`grep ^FROM Dockerfile | cut -d: -f 2 | tr -d " "`-`grep Version DESCRIPTION | cut -d: -f 2 | tr -d " "`-`grep spccaAppVersion inst/spCCA/version.R | cut -d'"' -f2` .
```

### Running from a container

To run the resulting container, start with
`docker run --rm -p 3838:3838 sneumann/spcca:latest`

and point your browser to http://localhost:3838/

### Developing and debugging in a container

If you want to develop and debug stuff, you can build a container
on top of `spcca:latest` that has an added rstudio server.
First build using `docker build -t spcca-rstudio -f Dockerfile-rstudio .`
and then run via `docker run -it --rm -p 8787:8787 spcca-rstudio:latest`.
CAVEAT: the `Dockerfile-rstudio` specifies a fixed user/password combo
of `rstudio:rstudio`. Do not use in Production !

You can also pass a local directory with checked out Metabolite GWAS git tree
via the `docker run -v` argument.

