FROM continuumio/miniconda3

RUN set -x && \
  apt-get update && \
  apt-get install --no-install-recommends --no-install-suggests -y \
    linux-headers-amd64 build-essential libc-dev gcc ghostscript hmmer && \
  apt-get clean && \
rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

EXPOSE 5000

RUN mkdir /app/ && mkdir /app/src && mkdir /app/tests

COPY environment.yml /app
COPY setup.py /app
COPY MANIFEST.in /app
COPY src /app/src/
COPY tests /app/tests/
COPY init.sh /app

WORKDIR /app

RUN chmod u+x init.sh
RUN ["conda", "env", "create", "-f", "environment.yml"]
RUN ["/bin/bash", "-c", "source activate chimera && python setup.py install"]

RUN echo "source activate chimera" > ~/.bashrc
ENV PATH /opt/conda/envs/chimera/bin:$PATH

ENV FLASK_APP=chimera.flask:create_app
ENTRYPOINT ["/bin/bash", "init.sh"]
