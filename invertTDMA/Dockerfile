FROM julia:1.10

RUN apt-get update
RUN apt-get -y install fonts-dejavu-core xfonts-base zip  

ENV USER pluto
ENV USER_HOME_DIR /home/${USER}
ENV JULIA_DEPOT_PATH ${USER_HOME_DIR}/.julia
ENV NOTEBOOK_DIR ${USER_HOME_DIR}/notebooks
ENV JULIA_NUM_THREADS 100

RUN useradd -m -d ${USER_HOME_DIR} ${USER} \
    && mkdir ${NOTEBOOK_DIR}

COPY prestartup.jl ${USER_HOME_DIR}/
COPY startup.jl ${USER_HOME_DIR}/
COPY notebooks/*.* ${NOTEBOOK_DIR}/

RUN julia ${USER_HOME_DIR}/prestartup.jl \
    && chown -R ${USER} ${USER_HOME_DIR}

USER ${USER}
EXPOSE 1234
VOLUME ${NOTEBOOK_DIR}
WORKDIR ${NOTEBOOK_DIR}

CMD [ "julia", "/home/pluto/startup.jl" ]
