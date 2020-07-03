# G+Smo Docker images

Scripts for building G+Smo Docker images

|||
|--:|---|
|License|[MPL 2.0](https://www.mozilla.org/en-US/MPL/2.0/)|
|OS support|Linux, Windows, macOS|
|Docker pulls| ![Docker pulls](https://img.shields.io/docker/pulls/mmoelle1/gismo) |
|Repository|[mmoelle1/gismo](https://hub.docker.com/repository/docker/mmoelle1/gismo)|
|Status|completed|
|Developer|Matthias Möller|
|Maintainer|M.Moller@tudelft.nl|
|Last checked|03-07-2020|

***
__Table of content__
1. [Prebuilt Docker images](#prebuilt-docker-images)
2. [Using Docker images](#using-docker-images)
3. [Building Docker images](#using-docker-images)

***

## Prebuilt Docker images
The following Docker images are built nightly from the G+Smo stable branch:
- `latest-ubuntu16.04-gcc5-release` ![Docker Image Size (tag)](https://img.shields.io/docker/image-size/mmoelle1/gismo/latest-ubuntu16.04-gcc5-release) ([GCC](https://gcc.gnu.org/) 5.x, C++98, Release mode with OpenMP)
- `latest-ubuntu16.04-gcc5-debug` ![Docker Image Size (tag)](https://img.shields.io/docker/image-size/mmoelle1/gismo/latest-ubuntu16.04-gcc5-debug) ([GCC](https://gcc.gnu.org/) 5.x, C++98, Debug mode with OpenMP)
- `latest-ubuntu16.04-clang3.8-release` ![Docker Image Size (tag)](https://img.shields.io/docker/image-size/mmoelle1/gismo/latest-ubuntu16.04-clang3.8-release) ([Clang](http://llvm.org/) 3.8, C++98, Release mode)
- `latest-ubuntu16.04-clang3.8-debug` ![Docker Image Size (tag)](https://img.shields.io/docker/image-size/mmoelle1/gismo/latest-ubuntu16.04-clang3.8-debug) ([Clang](http://llvm.org/) 3.8, C++98, Debug mode)

- `latest-ubuntu18.04-gcc7-cxx11-release` ![Docker Image Size (tag)](https://img.shields.io/docker/image-size/mmoelle1/gismo/latest-ubuntu18.04-gcc7-cxx11-release) ([GCC](https://gcc.gnu.org/) 7.x, C++11, Release mode with OpenMP)
- `latest-ubuntu18.04-gcc7-cxx11-debug` ![Docker Image Size (tag)](https://img.shields.io/docker/image-size/mmoelle1/gismo/latest-ubuntu18.04-gcc7-cxx11-debug) ([GCC](https://gcc.gnu.org/) 7.x, C++11, Debug mode with OpenMP)
- `latest-ubuntu18.04-clang6.0-cxx11-release` ![Docker Image Size (tag)](https://img.shields.io/docker/image-size/mmoelle1/gismo/latest-ubuntu18.04-clang6.0-cxx11-release) ([Clang](http://llvm.org/) 6.x, C++11, Release mode)
- `latest-ubuntu18.04-clang6.0-cxx11-debug` ![Docker Image Size (tag)](https://img.shields.io/docker/image-size/mmoelle1/gismo/latest-ubuntu18.04-clang6.0-cxx11-debug) ([Clang](http://llvm.org/) 6.x, C++11, Debug mode)

- `latest-ubuntu20.04-gcc9-cxx14-release` ![Docker Image Size (tag)](https://img.shields.io/docker/image-size/mmoelle1/gismo/latest-ubuntu20.04-gcc9-cxx14-release) ([GCC](https://gcc.gnu.org/) 9.x, C++14, Release mode with OpenMP)
- `latest-ubuntu20.04-gcc9-cxx14-debug` ![Docker Image Size (tag)](https://img.shields.io/docker/image-size/mmoelle1/gismo/latest-ubuntu20.04-gcc9-cxx14-debug) ([GCC](https://gcc.gnu.org/) 9.x, C++14, Debug mode with OpenMP)
- `latest-ubuntu20.04-clang10-cxx14-release` ![Docker Image Size (tag)](https://img.shields.io/docker/image-size/mmoelle1/gismo/latest-ubuntu20.04-clang10-cxx14-release) ([Clang](http://llvm.org/) 10.x, C++14, Release mode)
- `latest-ubuntu20.04-clang10-cxx14-debug` ![Docker Image Size (tag)](https://img.shields.io/docker/image-size/mmoelle1/gismo/latest-ubuntu20.04-clang10-cxx14-debug) ([Clang](http://llvm.org/) 10.x, C++14, Debug mode)

Two additional Docker images are built nightly with MPI and [Trilinos](http://trilinos.org) support enabled:
- `latest-ubuntu16.04-gcc5-cxx11-trilinos-release` ![Docker Image Size (tag)](https://img.shields.io/docker/image-size/mmoelle1/gismo/latest-ubuntu16.04-gcc5-cxx11-trilinos-release) ([GCC](https://gcc.gnu.org/) 5.x, C++11, Release mode with OpenMP)
- `latest-ubuntu16.04-gcc5-cxx11-trilinos-debug` ![Docker Image Size (tag)](https://img.shields.io/docker/image-size/mmoelle1/gismo/latest-ubuntu16.04-gcc5-cxx11-trilinos-debug) ([GCC](https://gcc.gnu.org/) 5.x, C++11, Debug mode with OpenMP)

The detailed build configuration is as follows (see the [G+Smo Wiki](https://www.gs.jku.at/trac/gismo/wiki/public/Compiling#Configuringandbuilding) for a description of the different options):
```
Configuration:
  CMAKE_BUILD_TYPE        {Release,Debug}
  CMAKE_INSTALL_PREFIX    /usr/local
  GISMO_COEFF_TYPE        double
  GISMO_EXTRA_DEBUG       OFF
  GISMO_BUILD_LIB         ON
  GISMO_BUILD_EXAMPLES    ON
  GISMO_WITH_OPENMP       {ON,OFF}
  GISMO_WITH_SUPERLU      ON
  GISMO_WITH_TRILINOS     {ON,OFF}
  GISMO_WITH_UMFPACK      ON
```
&nbsp;
The directory structure of the G+Smo library is as follows:
- `$HOME/gismo`: G+Smo main directory
- `$HOME/gismo/build/bin`: pre-compiled executables of example and tutorial applications
- `$HOME/gismo/examples`: source code of example and tutorial applications
- `$HOME/gismo/filedata`: parameterisations and other configuration files
- `$HOME/gismo/src`: source code of the G+Smo library

For a full description of the source code directory structure see the [G+Smo Wiki](https://www.gs.jku.at/trac/gismo/wiki/public/Compiling#Sourcedirectorytree).

## Using Docker images
To install Docker for your platform (Windows, macOS, Linux, cloud platforms, etc.) follow the instructions at [docker.com](https://docs.docker.com/get-started/).

Once you have Docker installed, you can run one of the above images by executing the following command:
```
docker run --rm -ti mmoelle1/gismo:latest-ubuntu16.04-gcc5-release
```
&nbsp;
You will be logged in under username `gismo` with standard password `gismo`. The `--rm` flag to the `docker run` command ensures that the container is removed upon termination. The directory `$HOME/gismo/build/bin` is added to the `PATH` variable so that example and tutorial applications can be run directly, e.g.:
```
docker run --rm -ti mmoelle1/gismo:latest-ubuntu16.04-gcc5-release geometry_example
```

&nbsp;
If you want to access data that is stored in a directory on your host computer from within your Docker container you need to mount (`-v`) the directory on your host computer when starting the Docker container. Assume that you want to make your current working directory on your host computer (`$(pwd)`) accessible under `/home/gismo/shared` in the docker container you need start the docker image as follows:
```
docker run --rm -ti -v $(pwd):/home/gismo/shared mmoelle1/gismo:latest-ubuntu16.04-gcc5-release
```
&nbsp;
It is also possible to mount a directory on your host computer to a directory of the Docker container and make it the working (`-w`) so that files generated inside the Docker container are accessible on your host computer:
```
docker run --rm -ti -w/home/gismo/shared -v $(pwd):/home/gismo/shared mmoelle1/gismo:latest-ubuntu16.04-gcc5-release
```
&nbsp;
A good starting point for learning the G+Smo library is to walk through the [example](https://www.gs.jku.at/trac/gismo/wiki/public/Doxygen/Examples) and [tutorial](https://www.gs.jku.at/trac/gismo/wiki/public/Doxygen/Tutorials) applications.

## Building your own Docker images

Scripts for building G+Smo Docker images are provided in the directories
```
docker
+- ubuntu16.04
   +- Dockerfile
+- ubuntu18.04
   +- Dockerfile
+- ubuntu20.04
   +- Dockerfile
```

To build a G+Smo Docker image with the default parameters run
```
docker build -t repo/image:tag ubuntuXX.YY
```
where `repo` is the name of the repository at Dockerhub (can be left blank), `image` is the name of the image and `tag` is the tag. `XX.YY` can be `16.04`, `18.04` or `20.04`. All default parameters can be changed using the `--build-arg` flag. For instance, building G+Smo with the C++ standard 14 can be achieved by running
```
docker build --build-arg --build-arg CMAKE_CXX_STANDARD=14 -t repo/image:tag ubuntuXX.YY
```
For a complete list of parameters have a look at the `Dockerfile`. The admissible values are given as comments. Other values might work as well but have not been tested. Please note that the configuration `--build-arg GISMO_WITH_TRILINOS=ON` also requires MPI to be enabled (`--build-arg GISMO_WITH_MPI=ON`) since otherwise the compilation of the Trilinos libraries will fail.
