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
1. [Building Docker images](#building-docker-images)
2. [Prebuilt Docker images](#prebuilt-docker-images)
3. [Using Docker images](#using-docker-images)
4. [Building Docker images](#using-docker-images)

***

## Building Docker images

To build G+Smo Docker images, use the `run_build_gismo.sh` script
providing as parameter the name of the docker image following the
naming convention

`os-compiler-cxx-buildtype[-arch][-coeff_type][-index_type][-option0][-option1]-...`

_Mandatory configuration_

-  `os` is the name of the operating system including specific version number. Supported values are: `ubuntu16.04`, `ubuntu18.04`, and `ubuntu20.04`.
-  `compiler` is the name of the compiler including a specific version number. Supported values are: `clang3.5`-`clang10` and `gcc4.7`-`gcc10`. Note that not all compiler versions are available on all operating systems. Newer compiler versions typically require also newer operating systems.
-  `cxx` is the C++ standard. Supported values are: `cxx98`, `cxx11`, `cxx14`, `cxx17`, and `cxx20`
-  `buildtype` is the build type. Supported values are: `debug`, `release`, `debinfo`, and `minsize`.

_Optional configuration_

-  `arch` is the name of the architecture. Supported values are `auto`, `non`, `generic`, `core`, `merom`, `penryn`, `nehalem`, `westmere`, `sandybridge`, `ivybridge`, `haswell`, `broadwell`, `skylake`, `skylake-xeon`, `kabylake`, `cannonlake`, `scacadelake`, `cooperlake`, `icelake`, `icelake-xeon`, `tigerlake`, `alderlake`, `sapphirerapids`, `bonnell`, `silvermont`, `goldmont`, `goldmont-plus`, `tremont`, `knl`, `knm`, `atom`, `k8`, `k8-sse3`, `barcelona`, `istanbul`, `magny-cours`, `bulldozer`, `interlagos`, `piledriver`, `steamroller`, `excavator`, `amd14h`, `amd16h`, `zen`, `zen2`, and `zen3`. If not given then `generic` is used as default.
-  `coeff_type` is the global coefficient type. Supported values are: `float`, `double`, `longdouble`, `mpfr::mpreal`, `mpq_class`, `posit_2_0`, `posit_3_0`, `posit_3_1`, `posit_4_0`, `posit_4_1`, `posit_8_0`, `posit_8_1`, `posit_16_1`, `posit_32_2`, `posit_64_3`, `posit_128_4`, and `posit_256_5`. If not given then `double` is used as default.
-  `index_type` is the global index type. Supported values are: `int`, `int32_t`, `int64_t`, `long`, and `longlong`. If not given then `int` is used as default.

_Optional features_

-  `codipack`/`nocodipack` to enable/disable building G+Smo with CoDiPack support. Default is disabled.
-  `examples`/`noexamples` to enable/disable building of examples. Default is enabled.
-  `gmp`/`nogmp` to enable/disable building G+Smo with GMP support. Default is disabled.
-  `ipopt`/`noipopt` to enable/disable building G+Smo with IPOpt support. Default is disabled.
-  `lib`/`nolib` to enable/disable building G+Smo as library or and header-only mode. Default is enabled.
-  `mpfr`/`nompfr` to enable/disable building G+Smo with MPFR support. Default is disabled.
-  `mpi`/`nompi` to enable/disable building G+Smo with MPI support. Default is disabled. Not that Trilinos requires MPI to be enabled and will overwrite this setting.
-  `occ`/`noocc` to enable/disable building G+Smo with OpenCascade support. Default is disabled.
-  `omp`/`noomp` to enable/disable building G+Smo with OpenMP support. Default is disabled for Clang compiler and enabled for all other compilers.
-  `onurbs`/`noonurbs` to enable/disable building G+Smo with OpenNURBS support. Default is disabled.
-  `pardiso`/`nopardiso` to enable/disable building G+Smo with Pardiso support. Default is disabled.
-  `pastix`/`nopastix` to enable/disable building G+Smo with Pastix support. Default is disabled.
-  `pch`/`nopch` to enable/disable building G+Smo with support for pre-compiled headers. Default is disabled.
-  `psolid`/`nopsolid` to enable/disable building G+Smo with PSOLID support. Default is disabled.
-  `smesh`/`nosmesh` to enable/disable building G+Smo with Surface Mesh support. Default is disabled.
-  `spectra`/`nospectra` to enable/disable building G+Smo with Spectra support. Default is disabled.
-  `superlu`/`nosuperlu` to enable/disable building G+Smo with SuperLU support. Default is disabled.
-  `taucs`/`notaucs` to enable/disable building G+Smo with Taucs support. Default is disabled.
-  `trilinos`/`notrilinos` to enable/disable building G+Smo with Trilinos support. Default is disabled.
-  `umfpack`/`noumfpack` to enable/disable building G+Smo with UMFPACK support. Default is disabled.
-  `unittests`/`nounittests` to enable/disable building of unittests. Default is disabled.
-  `unum`/`nounum` to enable/disable building G+Smo with UNUM support. Default is disabled.

## Prebuilt Docker images
The following Docker images are built automatically for commits to the G+Smo stable branch:

- `ubuntu16.04-gcc5-cxx98-release` ![Docker Image Size (tag)](https://img.shields.io/docker/image-size/mmoelle1/gismo/latest-ubuntu16.04-gcc5-release) ([GCC](https://gcc.gnu.org/) 5.x, C++98, Release mode with OpenMP)

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
