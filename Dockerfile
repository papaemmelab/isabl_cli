FROM leukgen/docker-pcapcore:v0.1.1

# set mantainers
LABEL maintainers="github.com/jsmedmar, github.com/juanesarango"

# install python 3
RUN \
    apt-get update -yq && \
    apt-get -yq install software-properties-common && \
    add-apt-repository ppa:deadsnakes/ppa && \
    apt update -yq && \
    apt install -yq python3.6 && \
    curl https://bootstrap.pypa.io/get-pip.py | python3.6

# # install repo
COPY . /code
RUN pip install /code && rm -rf /code

# set entrypoint
ENTRYPOINT [ "isabl" ]
