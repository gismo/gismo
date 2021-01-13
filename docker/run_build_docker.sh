#!/bin/bash

# Options management
usage() { echo "Usage: $0 imagename

An open-source script to facilitate the building of G+Smo docker images

Options:

imagename: An image name following the naming convention
           os-compiler-cxx-buildtype-cpu[-option0][-option1]...[:arch]
           where
                   os : is the operating system, e.g. ubuntu16.04
             compiler : is the compiler, e.g. gcc5
                  cxx : is the C++ standard, e.g. cxx11
           buildstype : is the build type, e.g. release
                  cpu : is the cpu type, e.g. skylake
                 arch : is the architectore, e.g. amd64" 1>&2; exit 1; }

shift $((OPTIND-1))
IMAGE_NAME=$1
BASE_PATH=$(dirname $(readlink -f $0))

if [ -z "${IMAGE_NAME}" ]; then
    usage
fi

shopt -s extglob;
case "$IMAGE_NAME" in
    ?(*:)?(*\-)ubuntu16.04?(\-*))
        (DOCKER_REPO="mmoelle1/gismo"
         DOCKERFILE_PATH="${BASE_PATH}/ubuntu16.04/Dockerfile"
         source ${BASE_PATH}/ubuntu16.04/hooks/build)
        ;;
    ?(*:)?(*\-)ubuntu18.04?(\-*))
        (DOCKER_REPO="mmoelle1/gismo"
         DOCKERFILE_PATH="${BASE_PATH}/ubuntu18.04/Dockerfile"
         source ${BASE_PATH}/ubuntu18.04/hooks/build)
        ;;
    ?(*:)?(*\-)ubuntu20.04?(\-*))
        (DOCKER_REPO="mmoelle1/gismo"
         DOCKERFILE_PATH="${BASE_PATH}/ubuntu20.04/Dockerfile"
         source ${BASE_PATH}/ubuntu20.04/hooks/build)
        ;;
    *)
        echo "Unsupported OS"
        exit 1
        ;;
esac
shopt -u extglob;
