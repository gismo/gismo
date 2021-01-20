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
|Last checked|09-12-2020|

***
__Table of content__
1. [Building Docker images](#building-docker-images)
2. [Prebuilt Docker images](#prebuilt-docker-images)
3. [Using Docker images](#using-docker-images)

***

## Building Docker images

The `run_build_docker.sh` script provides an easy mechanism to create
customized Docker images. We assume that you have installed Docker as
described at [Docker.com](https://www.docker.com/get-started) and can
execute it with your user privileges. To build a customized Docker
image run

```bash
./run_build_docker.sh os-compiler-cxx-buildtype[-cpu][-coeff_type][-index_type][-option0][-option1][...][:arch][/branch]
```
| | Mandatory configuration |
|-|-|
| `os` | _Name of the operating system including specific version number._ Supported values are: `ubuntu16.04`, `ubuntu18.04`, and `ubuntu20.04` |
| `compiler` | _Name of the compiler including a specific version number._ Supported values are: `clang3.5`-`clang10` and `gcc4.7`-`gcc10`. Note that not all compiler versions are available on all operating systems. Newer compiler versions typically require also newer operating systems. |
| `cxx` | _C++ standard._ Supported values are: `cxx98`, `cxx11`, `cxx14`, `cxx17`, and `cxx20` |
| `buildtype` | _Build type._ Supported values are: `debug`, `release`, `debinfo`, and `minsize` |

| | Optional configuration |
|-|-|
| `cpu`| _Name of the cpu architecture._ Supported values are `auto`, `none`, `generic`, `core`, `merom`, `penryn`, `nehalem`, `westmere`, `sandybridge`, `ivybridge`, `haswell`, `broadwell`, `skylake`, `skylake-xeon`, `kabylake`, `cannonlake`, `scacadelake`, `cooperlake`, `icelake`, `icelake-xeon`, `tigerlake`, `alderlake`, `sapphirerapids`, `bonnell`, `silvermont`, `goldmont`, `goldmont-plus`, `tremont`, `knl`, `knm`, `atom`, `k8`, `k8-sse3`, `barcelona`, `istanbul`, `magny-cours`, `bulldozer`, `interlagos`, `piledriver`, `steamroller`, `excavator`, `amd14h`, `amd16h`, `zen`, `zen2`, and `zen3`. If not given then `generic` is used as default. |
| `coeff_type`| _Global coefficient type._ Supported values are: `float`, `double`, `longdouble`, `mpfr::mpreal`, `mpq_class`, `posit_2_0`, `posit_3_0`, `posit_3_1`, `posit_4_0`, `posit_4_1`, `posit_8_0`, `posit_8_1`, `posit_16_1`, `posit_32_2`, `posit_64_3`, `posit_128_4`, and `posit_256_5`. If not given then `double` is used as default. |
| `index_type` | _Global index type._ Supported values are: `int`, `int32_t`, `int64_t`, `long`, and `longlong`. If not given then `int` is used as default. |

__Disclaimer__: Optional features are not yet fully supported. 

| | Optional features |
|-|-|
| `codipack`/<br>`nocodipack` | _Enable/disable building G+Smo with CoDiPack support._ Default is disabled. |
| `examples`/<br>`noexamples` | _Enable/disable building G+Smo examples._ Default is enabled. |
| `gmp`/<br>`nogmp` | _Enable/disable building G+Smo with GMP support._ Default is disabled. |
|`ipopt`/<br>`noipopt` | _Enable/disable building G+Smo with IPOpt support._ Default is disabled. |
| `lib`/<br>`nolib` | _Enable/disable building G+Smo as library._ Default is enabled. |
| `mpfr`/<br>`nompfr` | _Enable/disable building G+Smo with MPFR support._ Default is disabled. |
| `mpi`/<br>`nompi` | _Enable/disable building G+Smo with MPI support._ Default is disabled. Not that Trilinos requires MPI to be enabled and will overwrite this setting. |
| `occ`/<br>`noocc` | _Enable/disable building G+Smo with OpenCascade support._ Default is disabled. |
| `omp`/<br>`noomp` | _Enable/disable building G+Smo with OpenMP support._ Default is disabled for Clang compiler and enabled for all other compilers. |
| `onurbs`/<br>`noonurbs` | _Enable/disable building G+Smo with OpenNURBS support._ Default is disabled. |
| `pardiso`/<br>`nopardiso` | _Enable/disable building G+Smo with Pardiso support._ Default is disabled. |
| `pastix`/<br>`nopastix` | _Enable/disable building G+Smo with Pastix support._ Default is disabled. |
| `pch`/<br>`nopch` | _Enable/disable building G+Smo with support for pre-compiled headers._ Default is disabled. |
| `psolid`/<br>`nopsolid` | _Enable/disable building G+Smo with PSOLID support._ Default is disabled. |
| `smesh`/<br>`nosmesh` | _Enable/disable building G+Smo with Surface Mesh support._ Default is disabled. |
| `spectra`/<br>`nospectra` | _Enable/disable building G+Smo with Spectra support._ Default is disabled. |
| `superlu`/<br>`nosuperlu` | _Enable/disable building G+Smo with SuperLU support._ Default is disabled. |
| `taucs`/<br>`notaucs` | _Enable/disable building G+Smo with Taucs support._ Default is disabled. |
| `trilinos`/<br>`notrilinos` | _Enable/disable building G+Smo with Trilinos support._ Default is disabled. |
| `umfpack`/<br>`noumfpack` | _Enable/disable building G+Smo with UMFPACK support._ Default is disabled. |
| `unittests`/<br>`nounittests` | _Enable/disable building G+Smo unittests._ Default is disabled. |
| `unum`/<br>`nounum` | _Enable/disable building G+Smo with UNUM support._ Default is disabled. |

## Prebuilt Docker images
The following Docker images are built automatically for commits to the G+Smo stable branch:

- `ubuntu16.04-gcc5-cxx98-release` ![Docker Image Size (tag)](https://img.shields.io/docker/image-size/mmoelle1/gismo/ubuntu16.04-gcc5-cxx98-release) ([GCC](https://gcc.gnu.org/) 5.x, C++98, Release mode with OpenMP)
- `ubuntu18.04-gcc7-cxx11-release` ![Docker Image Size (tag)](https://img.shields.io/docker/image-size/mmoelle1/gismo/ubuntu18.04-gcc7-cxx11-release) ([GCC](https://gcc.gnu.org/) 7.x, C++11, Release mode with OpenMP)
- `ubuntu20.04-gcc10-cxx14-release` ![Docker Image Size (tag)](https://img.shields.io/docker/image-size/mmoelle1/gismo/ubuntu20.04-gcc10-cxx14-release) ([GCC](https://gcc.gnu.org/) 10.x, C++14, Release mode with OpenMP)

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
docker run --rm -ti mmoelle1/gismo:ubuntu18.04-gcc7-cxx11-release
```
&nbsp;
You will be logged in under username `gismo` with standard password `gismo`. The `--rm` flag to the `docker run` command ensures that the container is removed upon termination. The directory `$HOME/gismo/build/bin` is added to the `PATH` variable so that example and tutorial applications can be run directly, e.g.:
```
docker run --rm -ti mmoelle1/gismo:ubuntu18.04-gcc7-cxx11-release geometry_example
```

&nbsp;
If you want to access data that is stored in a directory on your host computer from within your Docker container you need to mount (`-v`) the directory on your host computer when starting the Docker container. Assume that you want to make your current working directory on your host computer (`$(pwd)`) accessible under `/home/gismo/shared` in the docker container you need start the docker image as follows:
```
docker run --rm -ti -v $(pwd):/home/gismo/shared mmoelle1/gismo:ubuntu18.04-gcc7-cxx11-release
```
&nbsp;
It is also possible to mount a directory on your host computer to a directory of the Docker container and make it the working (`-w`) so that files generated inside the Docker container are accessible on your host computer:
```
docker run --rm -ti -w/home/gismo/shared -v $(pwd):/home/gismo/shared mmoelle1/gismo:ubuntu18.04-gcc7-cxx11-release
```
&nbsp;
A good starting point for learning the G+Smo library is to walk through the [example](https://www.gs.jku.at/trac/gismo/wiki/public/Doxygen/Examples) and [tutorial](https://www.gs.jku.at/trac/gismo/wiki/public/Doxygen/Tutorials) applications.
