FROM ubuntu:16.04

# Update base image ubuntu
RUN apt-get update && apt-get install -y --no-install-recommends apt-utils
RUN apt-get upgrade -y
RUN apt-get install -y build-essential
RUN apt-get install -y git
RUN apt-get install -y wget
RUN apt-get install -y autoconf
RUN apt-get install -y gfortran
RUN apt-get install -y libblas-dev liblapack-dev  libpthread-stubs0-dev
RUN apt-get install -y pkg-config
RUN apt-get install -y libtool
RUN apt-get install -y libboost-all-dev

# CMAKE 3.10
RUN wget https://cmake.org/files/v3.10/cmake-3.10.3-Linux-x86_64.sh; chmod +x cmake-3.10.3-Linux-x86_64.sh; ./cmake-3.10.3-Linux-x86_64.sh --skip-license; ln -s ./cmake-3.10.3-Linux-x86_64/bin/cmake /usr/bin/cmake;

# Get dependencies
RUN git clone https://github.com/google/glog; cd glog; cmake -H. -Bbuild -G "Unix Makefiles"; cmake --build build --target install;
