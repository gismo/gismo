#!/bin/bash

# Options management
usage() { echo "Usage: $0 imagename

An open-source script to facilitate the building of G+Smo docker images

Options:

imagename: An image name following the naming convention
           os-compiler-cxx-buildtype-cpu[-option0][-option1]...
           where
                   os : is the operating system, e.g. ubuntu20.04
             compiler : is the compiler, e.g. gcc10
                  cxx : is the C++ standard, e.g. cxx11
           buildstype : is the build type, e.g. release
                  cpu : is the cpu type, e.g. skylake
              optionX : is an optional option, e.g. openmp" 1>&2; exit 1; }

realpath() {
    [[ $1 = /* ]] && echo "$1" || echo "$PWD/${1#./}"
}

shift $((OPTIND-1))
IMAGE_NAME=$1
BASE_DIR=$(dirname $(realpath $0))

if [ -z "${IMAGE_NAME}" ]; then
    usage
fi

shopt -s extglob;
case "$IMAGE_NAME" in
    ?(*:)?(*\-)ubuntu16.04?(\-*))
        (DOCKER_REPO="mmoelle1/gismo"
         DOCKERFILE="${BASE_DIR}/ubuntu16.04/Dockerfile"
         source ${BASE_DIR}/ubuntu16.04/hooks/build)
        ;;
    ?(*:)?(*\-)ubuntu18.04?(\-*))
        (DOCKER_REPO="mmoelle1/gismo"
         DOCKERFILE="${BASE_DIR}/ubuntu18.04/Dockerfile"
         source ${BASE_DIR}/ubuntu18.04/hooks/build)
        ;;
    ?(*:)?(*\-)ubuntu20.04?(\-*))
        (DOCKER_REPO="mmoelle1/gismo"
         DOCKERFILE="${BASE_DIR}/ubuntu20.04/Dockerfile"
         source ${BASE_DIR}/ubuntu20.04/hooks/build)
        ;;
    *)
        echo "Unsupported OS"
        exit 1
        ;;
esac
shopt -u extglob;
