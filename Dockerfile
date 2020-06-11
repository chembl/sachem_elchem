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

# Prepare ANT for cpptask build
ENV JAVA_HOME /usr/share/elasticsearch/jdk
COPY ./ ${WORKDIR}/
RUN mkdir jni-build
RUN mkdir -p target/META-INF/isomorphism/2.0/LINUX-AMD64
RUN gcc -c -fPIC -I${JAVA_HOME}/include -I${JAVA_HOME}/include/linux -std=c99 -o jni-build/libisomorphism.o jni/native.c
RUN gcc -shared -fPIC -o target/META-INF/isomorphism/2.0/LINUX-AMD64/libisomorphism.so jni-build/libisomorphism.o -lc

# Build elchem.jar and elchem.zip
RUN ant

# Install plugin in elasticsearch
RUN elasticsearch-plugin install -b file:///${WORKDIR}/elchem.zip

ENTRYPOINT ["/usr/local/bin/docker-entrypoint.sh"]

