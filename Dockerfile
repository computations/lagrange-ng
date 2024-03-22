FROM ubuntu:22.04
RUN apt-get update && \
apt-get upgrade --assume-yes && \
DEBIAN_FRONTEND=noninteractive TZ=Etc/UTC \
apt-get install --assume-yes build-essential gfortran git cmake && mkdir lagrange-ng
RUN git clone https://github.com/computations/lagrange-ng/ --depth=1 --recursive
RUN cd lagrange-ng && cmake -Bbuild -H. -DCMAKE_BUILD_TYPE=Release && cd build && make && cp ../bin/lagrange-ng /usr/local/bin/
ENTRYPOINT ["lagrange-ng"]
