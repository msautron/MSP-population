#Basic layer> Ubuntu
# Use phusion base image to get a basic Ubuntu envorinemnt and other additional stuff to mantain the
# container running
FROM ubuntu:jammy
MAINTAINER giuliano.iorio.astro@gmail.com
ENV WORKDIR=/home/sevnuser
WORKDIR $WORKDIR

ARG USERNAME=sevnuser
ARG USER_UID=1000
ARG USER_GID=$USER_UID

# Create the user
RUN groupadd --gid $USER_GID $USERNAME \
    && useradd --uid $USER_UID --gid $USER_GID -m $USERNAME
# [Optional] Add sudo support. Enable this at your own risk.
#RUN    apt-get update \
    	#&& apt-get install -y sudo \
    	#&& echo $USERNAME ALL=\(root\) NOPASSWD:ALL > /etc/sudoers.d/$USERNAME \
    	#&& chmod 0440 /etc/sudoers.d/$USERNAME


#1 layer: install g++, git, cmake
RUN   apt update &&   apt -y dist-upgrade \
&&   apt update \
&&   apt install -y \
g++ \
cmake \
git \
wget \
vim \
nano \
sudo \
emacs \
libgsl-dev

RUN  chown -R $USERNAME $WORKDIR
USER $USERNAME

#2 layer: SEVN \
# clone and build SEVN fork,
# the variable SEVN_VER is used to breack the Docker cache and alreay rebuilt this layer so
# that we install always the last version
# clone and build SEVN fork,
# the variable SEVN_VER is used to breack the Docker cache and alreay rebuilt this layer so
# that we install always the last version.
# It is rebuilt everytime the SEVN_VER change, so to rebuild it everytime it is enough
# just to use a timestamp, e.g. --build-arg SEVN_VER=$(date +%Y%m%d-%H%M%S)
ARG SEVN_VER=unknown
RUN  git clone https://gitlab.com/sevncodes/sevn.git \
&& cd sevn && mkdir build &&  cd build && cmake .. && make \
&& cd  ../.. && cp sevn/run_scripts/run.sh .



CMD ["tail", "-f", "/dev/null"]
