FROM continuumio/miniconda3

EXPOSE 5000

RUN mkdir /app/ && mkdir /app/src && mkdir /app/tests

COPY environment.yml /app
COPY setup.py /app
COPY src /app/src/
COPY tests /app/tests/

WORKDIR /app

RUN ["conda", "env", "create", "-f", "environment.yml"]
RUN ["/bin/bash", "-c", "source activate chimera && python setup.py install"]

RUN echo "source activate chimera" > ~/.bashrc
ENV PATH /opt/conda/envs/chimera/bin:$PATH

ENV FLASK_APP=chimera.flask:create_app
ENTRYPOINT ["/bin/bash", "-c"]
CMD ["source activate chimera && exec python -m flask run --host=0.0.0.0"]
