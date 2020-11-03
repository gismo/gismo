# This file creates docker containers of the G+Smo library
#
# G+Smo (pronounced gismo or gizmo) is a C++ library for isogeometric
# analysis (IGA). Geometry plus simulation modules aims at the
# seamless integration of Computer-aided Design (CAD) and Finite
# Element Analysis (FEA).
#
# https://github.com/gismo
#
# The docker containers provide all prerequisite software tools,
# compilers and libraries and a complete build of the G+Smo library.
FROM ubuntu:20.04
MAINTAINER Matthias Moller <m.moller@tudelft.nl>

ARG BUILD_RFC3339="1970-01-01T00:00:00Z"
ARG COMMIT="local"
ARG VERSION="dirty"

STOPSIGNAL SIGKILL

LABEL org.opencontainers.image.ref.name="gismo/gismo" \
      org.opencontainers.image.created=$BUILD_RFC3339 \
      org.opencontainers.image.authors="Matthias Moller <m.moller@tudelft.nl>" \
      org.opencontainers.image.documentation="https://github.com/gismo/gismo/README.md" \
      org.opencontainers.image.description="G+Smo docker image" \
      org.opencontainers.image.licenses="MPL-2.0" \
      org.opencontainers.image.source="https://github.com/gismo/gismo" \
      org.opencontainers.image.revision=$COMMIT \
      org.opencontainers.image.version=$VERSION \
      org.opencontainers.image.url="https://hub.docker.com/r/mmoelle1/gismo/"

ENV BUILD_RFC3339 "$BUILD_RFC3339"
ENV COMMIT "$COMMIT"
ENV VERSION "$VERSION"

# This file accepts the following build-time arguments:

# C compiler              : {clang-6.0,clang-7,clang-8,clang-9,clang-10,gcc-7,gcc-8,gcc-9,gcc-10}
ARG CC=gcc-9
# C++ compiler            : {clang++-6.0,clang++-7,clang++-8,clang++-9,clang++-10,g++-7,g++-8,g++-9,g++-10}
ARG CXX=g++-9

# Fortran compiler        : {gfortran-7,gfortran-8,gfortran-9,gfortran-10}
ARG FC=na

# CMAKE_CXX_STANDARD      : {98,11,14,17,20}
ARG CMAKE_CXX_STANDARD=14

# CMAKE_BUILD_TYPE        : {Release,Debug,RelWithDebInfo,MinSizeRel}
ARG CMAKE_BUILD_TYPE=Release

# CMAKE_INSTALL_PREFIX
ARG CMAKE_INSTALL_PREFIX=/home/gismo/gismo

# GISMO_COEFF_TYPE        : {double,float}
ARG GISMO_COEFF_TYPE=double

# GISMO_VERSION           : {v0.8.1,v0.8.2,v0.8.3,v0.8.4,HEAD}
ARG GISMO_VERSION=HEAD

# GISMO_WITH_MPI          : {ON,OFF}
ARG GISMO_WITH_MPI=OFF

# GISMO_WITH_OPENMP       : {ON,OFF}
ARG GISMO_WITH_OPENMP=ON

# GISMO_WITH_SUPERLU      : {ON,OFF}
ARG GISMO_WITH_SUPERLU=ON

# GISMO_WITH_TRILINOS     : {ON,OFF}
ARG GISMO_WITH_TRILINOS=OFF

# GISMO_WITH_UMFPACK      : {ON,OFF}
ARG GISMO_WITH_UMFPACK=ON

# Target architecture     : auto by default
ARG TARGET_ARCHITECTURE=auto

# Set the variable DEBIAN_FRONTEND to noninteractive. Don't use ENV
# since this sets the variable DEBIAN_FRONTEND also in the container.
ARG DEBIAN_FRONTEND=noninteractive

# Install prerequisite software
RUN apt-get update -q &&                        \
    apt-get install --no-install-recommends -yq \
    ca-certificates                             \
    cmake                                       \
    git                                         \
    libboost-all-dev                            \
    locales                                     \
    make                                        \
    software-properties-common                  \
    sudo                                        \
    wget                                        \
    zlib1g-dev      &&                          \
    apt-get clean   &&                          \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# Install MPI
RUN if [ "$GISMO_WITH_MPI" = "ON" ] ;              \
    then                                           \
    apt-get update -q &&                           \
    apt-get install --no-install-recommends -yq    \
        libopenmpi-dev &&                          \
    apt-get clean &&                               \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*; \
    fi

# Install SUPERLU
RUN if [ "$GISMO_WITH_SUPERLU" = "ON" ] ;          \
    then                                           \
    apt-get update -q &&                           \
    apt-get install --no-install-recommends -yq    \
        libopenblas-base                           \
        libopenblas-dev                            \
        liblapacke-dev                             \
        libsuperlu-dev &&                          \
    apt-get clean &&                               \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*; \
    fi
ENV CMAKE_INCLUDE_PATH=$CMAKE_INCLUDE_PATH:/usr/include/superlu

# Install TRILINOS
RUN if [ "$GISMO_WITH_TRILINOS" = "ON" ] ;         \
    then                                           \
    apt-get update -q &&                           \
    apt-get install --no-install-recommends -yq    \
        bc                                         \
        libopenblas-base                           \
        libopenblas-dev                            \
        liblapacke-dev &&                          \
    apt-get clean &&                               \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*; \
    fi

# Install UMFPACK
RUN if [ "$GISMO_WITH_UMFPACK" = "ON" ] ;          \
    then                                           \
    apt-get update -q &&                           \
    apt-get install --no-install-recommends -yq    \
        libopenblas-base                           \
        libopenblas-dev                            \
        liblapacke-dev                             \
        libsuitesparse-dev &&                      \
    apt-get clean &&                               \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*; \
    fi
ENV CMAKE_INCLUDE_PATH=$CMAKE_INCLUDE_PATH:/usr/include/suitesparse

# Install Clang C/C++ compiler version 6.0 (if required)
RUN if [ "$CC" = "clang-6.0" ] || [ "$CXX" = "clang++-6.0" ] ; \
    then                                               \
    apt-get update -q &&                               \
    apt-get install --no-install-recommends -yq        \
        clang-6.0 &&                                   \
    apt-get clean &&                                   \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*;     \
    fi

# Install Clang C/C++ compiler version 7 (if required)
RUN if [ "$CC" = "clang-7" ] || [ "$CXX" = "clang++-7" ] ; \
    then                                               \
    apt-get update -q &&                               \
    apt-get install --no-install-recommends -yq        \
        clang-7 &&                                     \
    apt-get clean &&                                   \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*;     \
    fi

# Install Clang C/C++ compiler version 8 (if required)
RUN if [ "$CC" = "clang-8" ] || [ "$CXX" = "clang++-8" ] ; \
    then                                               \
    apt-get update -q &&                               \
    apt-get install --no-install-recommends -yq        \
        clang-8 &&                                     \
    apt-get clean &&                                   \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*;     \
    fi

# Install Clang C/C++ compiler version 9 (if required)
RUN if [ "$CC" = "clang-9" ] || [ "$CXX" = "clang++-9" ] ; \
    then                                               \
    apt-get update -q &&                               \
    apt-get install --no-install-recommends -yq        \
        clang-9 &&                                     \
    apt-get clean &&                                   \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*;     \
    fi

# Install Clang C/C++ compiler version 10 (if required)
RUN if [ "$CC" = "clang-10" ] || [ "$CXX" = "clang++-10" ] ; \
    then                                               \
    apt-get update -q &&                               \
    apt-get install --no-install-recommends -yq        \
        clang-10 &&                                    \
    apt-get clean &&                                   \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*;     \
    fi

# Install GCC C compiler version 7 (if required)
RUN if [ "$CC" = "gcc-7" ] ; then                  \
    apt-get update -q &&                           \
    apt-get install --no-install-recommends -yq    \
        gcc-7 &&                                   \
    apt-get clean &&                               \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*; \
    fi

# Install GCC C compiler version 8 (if required)
RUN if [ "$CC" = "gcc-8" ] ; then                  \
    apt-get update -q &&                           \
    apt-get install --no-install-recommends -yq    \
        gcc-8 &&                                   \
    apt-get clean &&                               \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*; \
    fi

# Install GCC C compiler version 9 (if required)
RUN if [ "$CC" = "gcc-9" ] ; then                  \
    apt-get update -q &&                           \
    apt-get install --no-install-recommends -yq    \
        gcc-9 &&                                   \
    apt-get clean &&                               \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*; \
    fi

# Install GCC C compiler version 10 (if required)
RUN if [ "$CC" = "gcc-10" ] ; then                 \
    apt-get update -q &&                           \
    apt-get install --no-install-recommends -yq    \
        gcc-10 &&                                  \
    apt-get clean &&                               \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*; \
    fi

# Install GCC C++ compiler version 7 (if required)
RUN if [ "$CXX" = "g++-7" ] ; then                 \
    apt-get update -q &&                           \
    apt-get install --no-install-recommends -yq    \
        g++-7 &&                                   \
    apt-get clean &&                               \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*; \
    fi

# Install GCC C++ compiler version 8 (if required)
RUN if [ "$CXX" = "g++-8" ] ; then                 \
    apt-get update -q &&                           \
    apt-get install --no-install-recommends -yq    \
        g++-8 &&                                   \
    apt-get clean &&                               \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*; \
    fi

# Install GCC C++ compiler version 9 (if required)
RUN if [ "$CXX" = "g++-9" ] ; then                 \
    apt-get update -q &&                           \
    apt-get install --no-install-recommends -yq    \
        g++-9 &&                                   \
    apt-get clean &&                               \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*; \
    fi

# Install GCC C++ compiler version 10 (if required)
RUN if [ "$CXX" = "g++-10" ] ; then                \
    apt-get update -q &&                           \
    apt-get install --no-install-recommends -yq    \
        g++-10 &&                                  \
    apt-get clean &&                               \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*; \
    fi

# Install GCC Fortran compiler version 7 (if required)
RUN if [ "$FC" = "gfortran-7" ] ; then             \
    apt-get update -q &&                           \
    apt-get install --no-install-recommends -yq    \
        gfortran-7 &&                              \
    ln -nsf /usr/bin/gfortran-7 /usr/bin/gfortran && \
    apt-get clean &&                               \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*; \
    fi

# Install GCC Fortran compiler version 8 (if required)
RUN if [ "$FC" = "gfortran-8" ] ; then             \
    apt-get update -q &&                           \
    apt-get install --no-install-recommends -yq    \
        gfortran-8 &&                              \
    ln -nsf /usr/bin/gfortran-8 /usr/bin/gfortran && \
    apt-get clean &&                               \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*; \
    fi

# Install GCC Fortran compiler version 9 (if required)
RUN if [ "$FC" = "gfortran-9" ] ; then             \
    apt-get update -q &&                           \
    apt-get install --no-install-recommends -yq    \
        gfortran-9 &&                              \
    ln -nsf /usr/bin/gfortran-9 /usr/bin/gfortran && \
    apt-get clean &&                               \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*; \
    fi

# Install GCC Fortran compiler version 10 (if required)
RUN if [ "$FC" = "gfortran-10" ] ; then            \
    apt-get update -q &&                           \
    apt-get install --no-install-recommends -yq    \
        gfortran-10 &&                             \
    apt-get clean &&                               \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*; \
    fi

# Set the locale
RUN locale-gen en_US.UTF-8  
ENV LANG en_US.UTF-8  
ENV LANGUAGE en_US:en  
ENV LC_ALL en_US.UTF-8

# Add and enable the default user
ARG USER=gismo
RUN adduser --disabled-password --gecos '' $USER
RUN adduser $USER sudo; echo '%sudo ALL=(ALL) NOPASSWD:ALL' >> /etc/sudoers

# Make sure everything is in place
RUN chown -R $USER:$USER /home/$USER
USER $USER
ENV HOME /home/$USER
ENV USER $USER
WORKDIR $HOME

# Check out the open-source stable version of G+Smo
RUN if [ "$GISMO_VERSION" = "HEAD" ] ; \
    then \
    git clone https://github.com/gismo/gismo.git gismo; \
    else \
    git clone https://github.com/gismo/gismo.git gismo && cd gismo && git fetch --tags && git checkout $GISMO_VERSION; \
    fi

# Configure and build G+Smo library
RUN cd gismo    &&                                  \
    mkdir build &&                                  \
    cd build    &&                                  \
    CC=$CC                                          \
    CXX=$CXX                                        \
    FC=$FC                                          \
    cmake ..                                        \
        -DTARGET_ARCHITECTURE=$TARGET_ARCHITECTURE  \
        -DCMAKE_BUILD_TYPE=$CMAKE_BUILD_TYPE        \
        -DCMAKE_CXX_STANDARD=$CMAKE_CXX_STANDARD    \
        -DCMAKE_INSTALL_PREFIX=$CMAKE_INSTALL_PREFIX\        
        -DGISMO_COEFF_TYPE=$GISMO_COEFF_TYPE        \
        -DGISMO_WITH_MPI=$GISMO_WITH_MPI            \
        -DGISMO_WITH_OPENMP=$GISMO_WITH_OPENMP      \
        -DGISMO_WITH_SUPERLU=$GISMO_WITH_SUPERLU    \
        -DGISMO_WITH_TRILINOS=$GISMO_WITH_TRILINOS  \
        -DGISMO_WITH_UMFPACK=$GISMO_WITH_UMFPACK && \
    make

# Add directory of executables to search path
ENV PATH $PATH:$CMAKE_INSTALL_PREFIX:$HOME/gismo/build/bin
