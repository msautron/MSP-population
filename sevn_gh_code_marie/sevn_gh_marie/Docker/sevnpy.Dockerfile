#Basic layer> Ubuntu
FROM ubuntu:jammy
MAINTAINER giuliano.iorio.astro@gmail.com
ENV WORKDIR=/home/sevnuser
WORKDIR $WORKDIR
ENV PATH="$WORKDIR/miniconda3/bin:$PATH"
ARG TARGETPLATFORM

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
RUN  apt update && apt -y dist-upgrade \
&& apt update \
&& apt install -y \
g++ \
cmake \
git \
wget \
sudo \
vim \
nano \
emacs \
libgsl-dev

RUN sudo chown -R $USERNAME $WORKDIR
USER $USERNAME

#2 layer: miniconda \
RUN case "${TARGETPLATFORM}" in \
         "linux/amd64")  export MINICONDAS=Miniconda3-py310_23.5.2-0-Linux-x86_64.sh ;; \
         "linux/arm64")  export MINICONDAS=Miniconda3-py310_23.5.2-0-Linux-aarch64.sh  ;; \
    esac; \
     cd $WORKDIR \
     && wget \
     https://repo.anaconda.com/miniconda/${MINICONDAS} \
     &&  bash ${MINICONDAS}  -b \
     && rm -f ${MINICONDAS}


#3 create env and basic python packages \
RUN  pip install numpy pandas matplotlib scipy scikit-learn jupyter

#4 layer: SEVN \
# clone and build SEVN fork,
# the variable SEVN_VER is used to breack the Docker cache and alreay rebuilt this layer so
# that we install always the last version.
# It is rebuilt everytime the SEVN_VER change, so to rebuild it everytime it is enough
# just to use a timestamp, e.g. --build-arg SEVN_VER=$(date +%Y%m%d-%H%M%S)
ARG SEVN_VER=unknown
RUN  git clone https://gitlab.com/sevncodes/sevn.git \
&& cd sevn && mkdir build &&  cd build && cmake .. && make \
&& cd  ../.. && cp sevn/run_scripts/run.sh .

#5 Install SEVNpy
RUN cd sevn/SEVNpy && pip install .

#6 Setup easy jupyter access
RUN mkdir -p $HOME/.jupyter/ \
    && echo "c.NotebookApp.token=''" >> $HOME/.jupyter/jupyter_notebook_config.py;


CMD ["jupyter", "notebook", "--port=8888", "--no-browser", "--ip=0.0.0.0", "--allow-root"]

EXPOSE 8888