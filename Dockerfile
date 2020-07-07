FROM docker.elastic.co/elasticsearch/elasticsearch:7.5.2
ARG ELCHEM_DIR=${ELCHEM_DIR:-/usr/share/sachem-elchem}
# setup base directory
RUN mkdir -p ${ELCHEM_DIR}

RUN yum -y update
RUN yum -y install ant
RUN yum -y install gcc
RUN yum -y install make

# Prepare ANT for cpptask build
ENV JAVA_HOME /usr/share/elasticsearch/jdk
COPY ./ ${ELCHEM_DIR}/

# Build elchem.jar and elchem.zip
RUN cd ${ELCHEM_DIR} && ant

# Install plugin in elasticsearch
RUN elasticsearch-plugin install -b file:///${ELCHEM_DIR}/elchem.zip

ENTRYPOINT ["/usr/local/bin/docker-entrypoint.sh"]

