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
RUN yum -y install automake
RUN yum -y install libtool
RUN yum -y install ant
RUN yum -y install wget
RUN wget https://download.postgresql.org/pub/repos/yum/9.6/redhat/rhel-7-x86_64/pgdg-redhat-repo-42.0-9.noarch.rpm
RUN yum -y install pgdg-redhat-repo-42.0-9.noarch.rpm epel-release
RUN yum -y update
RUN yum -y install postgresql96-server postgresql96-contrib postgresql96-devel

COPY ./ ${WORKDIR}/
RUN mkdir build
ENV PATH /usr/pgsql-9.6/bin:$PATH
ENV JAVA_HOME /usr/lib/jvm/java-1.8.0-openjdk-1.8.0.242.b08-1.el7.x86_64
RUN aclocal
RUN autoconf
RUN autoreconf --install
RUN automake -a --foreign
RUN cd build && ../configure

ENTRYPOINT ["/usr/local/bin/docker-entrypoint.sh"]

