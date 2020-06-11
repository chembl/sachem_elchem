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

COPY ./ ${WORKDIR}/
RUN mkdir build
ENV JAVA_HOME /usr/lib/jvm/java-1.8.0-openjdk-1.8.0.242.b08-1.el7.x86_64
#RUN ant

ENTRYPOINT ["/usr/local/bin/docker-entrypoint.sh"]

