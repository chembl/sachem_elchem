FROM docker.elastic.co/elasticsearch/elasticsearch:7.5.2
ARG USERNAME=${USERNAME:-elastic}
ARG UID=${UID:-1000}
ARG GID=${GID:-1000}
ARG WORKDIR=${WORKDIR:-/sachem}
# setup user and app root directory
RUN mkdir -p ${WORKDIR}
RUN chown -R ${UID}:${GID} ${WORKDIR}
WORKDIR ${WORKDIR}

RUN yum -y update
RUN yum -y install ant
RUN yum -y install gcc
RUN yum -y install make

# Prepare ANT for cpptask build
ENV JAVA_HOME /usr/share/elasticsearch/jdk
COPY ./ ${WORKDIR}/

# Build elchem.jar and elchem.zip
RUN ant

# Install plugin in elasticsearch
RUN elasticsearch-plugin install -b file:///${WORKDIR}/elchem.zip

ENTRYPOINT ["/usr/local/bin/docker-entrypoint.sh"]

