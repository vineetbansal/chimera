FROM continuumio/miniconda3:4.6.14

# ----------------------------------
# Env vars needed during build
# These are passed to the running container
# ----------------------------------
ENV CHIMERA_DIRS_DATA_PFAM=/pfam
# ----------------------------------

ARG APP_DIR=/app

RUN set -x && \
  apt-get update && \
  apt-get install --no-install-recommends --no-install-suggests -y \
    linux-headers-amd64 build-essential libc-dev gcc ghostscript hmmer liblpsolve55-dev lp-solve && \
  apt-get clean && \
rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# lp-solve installs the shared libraries, liblpsolve55-dev installs the headers and static libraries
# At runtime, the static libraries are preferred over the .so unless we do the following.
RUN ln -s /usr/lib/lp_solve/liblpsolve55.so /usr/lib/liblpsolve55.so

RUN mkdir $APP_DIR/ && mkdir $APP_DIR/src && mkdir $APP_DIR/tests

COPY environment.yml $APP_DIR
WORKDIR $APP_DIR
RUN ["conda", "env", "create", "-f", "environment.yml"]

# dPuc2 installation
# Configure cpan for first use
RUN ["/bin/bash", "-c", "echo | cpan"]
RUN ["/bin/bash", "-c", "cpan install Inline::C"]
RUN ["git", "clone", "https://github.com/alexviiia/dpuc2.git", "/dpuc2"]
RUN ["git", "clone", "https://github.com/alexviiia/dpuc2-data.git", "/dpuc2-data"]

# DomStratStats installation
RUN ["git", "clone", "https://github.com/alexviiia/DomStratStats.git", "/DomStratStats"]

COPY setup.py $APP_DIR
COPY app.py $APP_DIR
COPY MANIFEST.in $APP_DIR
COPY src $APP_DIR/src/
COPY tests $APP_DIR/tests/

RUN ["/bin/bash", "-c", "source activate chimera && python setup.py install"]

RUN echo "source activate chimera" > ~/.bashrc
ENV PATH /opt/conda/envs/chimera/bin:$PATH

# ----------------------------------
# Env vars passed by this Dockerfile to the running container
# ----------------------------------
ENV CHIMERA_DIRS_BIN_DPUC2=/dpuc2
ENV CHIMERA_FILES_DATA_DPUC2_NET=/dpuc2-data/dpucNet.pfam32.txt.gz
ENV CHIMERA_DIRS_BIN_DOMSTRATSTATS=/DomStratStats
# Where apt-get installs hmmr by default
ENV CHIMERA_DIRS_BIN_HMMR=/usr/bin
# ----------------------------------

CMD ["/bin/bash", "-c", "python app.py"]
