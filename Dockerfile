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

# Prepare ANT for cpptask build
COPY ./ ${WORKDIR}/
RUN mkdir build
ENV JAVA_HOME /usr/share/elasticsearch/jdk
ENV CLASSPATH ${WORKDIR}/ant-ext/ant-contrib-0.6-bin/lib:${WORKDIR}/ant-ext/cpptasks-1.0b4/lib
RUN cd ant-ext
RUN cd ant-ext && tar -xvf ant-contrib-0.6-bin.tar.gz
RUN cd ant-ext && tar -xvf cpptasks-1.0b4.tar.gz
RUN cd ant-ext/cpptasks-1.0b4 && mkdir lib && mv cpptasks.jar ./lib

# build elchem
RUN ant

ENTRYPOINT ["/usr/local/bin/docker-entrypoint.sh"]

