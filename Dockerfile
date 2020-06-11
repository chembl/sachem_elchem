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
RUN cd ant-ext && mkdir ant-contrib-0.6-bin && mv ant-contrib-0.6-bin.tar.gz ant-contrib-0.6-bin
RUN cd ant-ext/ant-contrib-0.6-bin && tar -xvf ant-contrib-0.6-bin.tar.gz
RUN cd ant-ext && tar -xvf cpptasks-1.0b4.tar.gz
RUN cd ant-ext && mkdir lib
RUN ls -la ./
RUN ls -la ./ant-ext/
RUN ls -la ./ant-ext/ant-contrib-0.6-bin/
RUN ls -la ./ant-ext/ant-contrib-0.6-bin/lib/
RUN cp ./ant-ext/ant-contrib-0.6-bin/lib/ant-contrib-0.6.jar /usr/share/ant/lib
RUN cp ./ant-ext/cpptasks-1.0b4/cpptasks.jar /usr/share/ant/lib
RUN ls -la /usr/share/ant/lib

# build elchem
RUN ant

ENTRYPOINT ["/usr/local/bin/docker-entrypoint.sh"]

