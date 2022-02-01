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
|Last checked|21-01-2022|

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
| `os` | _Name of the operating system including specific version number._ Supported values are: `ubuntu16.04`, `ubuntu18.04`, and `ubuntu20.04`. |
| `compiler` | _Name of the compiler including a specific version number._ Supported values are: `clang3.5` to `clang11` and `gcc4.7` to `gcc11`. Note that not all compiler versions are available on all operating systems. Newer compiler versions typically require also newer operating systems. |
| `cxx` | _C++ standard._ Supported values are: `cxx11`, `cxx14`, `cxx17`, `cxx20`, and `cxx23`. |
| `buildtype` | _Build type._ Supported values are: `debug`, `release`, `debinfo`, and `minsize`. |

| | Optional configuration |
|-|-|
| `cpu`| _Name of the cpu architecture._ Supported values are `auto`, `none`, `generic`, `core`, `merom`, `penryn`, `nehalem`, `westmere`, `sandybridge`, `ivybridge`, `haswell`, `broadwell`, `skylake`, `skylake-xeon`, `kabylake`, `cannonlake`, `cacadelake`, `cooperlake`, `icelake`, `icelake-xeon`, `tigerlake`, `alderlake`, `sapphirerapids`, `bonnell`, `silvermont`, `goldmont`, `goldmont-plus`, `tremont`, `knl`, `knm`, `atom`, `k8`, `k8-sse3`, `barcelona`, `istanbul`, `magny-cours`, `bulldozer`, `interlagos`, `piledriver`, `steamroller`, `excavator`, `amd14h`, `amd16h`, `zen`, `zen2`, and `zen3`. <br> If not given then `generic` is used as default. |
| `coeff_type`| _Global coefficient type._ Supported values are: `float`, `double`, `longdouble`, `mpfr::mpreal`, `mpq_class`, `posit_2_0`, `posit_3_0`, `posit_3_1`, `posit_4_0`, `posit_4_1`, `posit_8_0`, `posit_8_1`, `posit_16_1`, `posit_32_2`, `posit_64_3`, `posit_128_4`, and `posit_256_5`. <br> If not given then `double` is used as default. |
| `index_type` | _Global index type._ Supported values are: `int`, `int32_t`, `int64_t`, `long`, and `longlong`. <br> If not given then `int` is used as default. |

__Disclaimer__: Not all optional features are supported on all architectures and/or OS versions. Please report errors to the maintainer.

| | Optional features |
|-|-|
| `adiff`/<br>`noadiff` | _Enable/disable building G+Smo with [AutoDiff](https://eigen.tuxfamily.org/dox/unsupported/group__AutoDiff__Module.html) support._ Default is disabled. |
| `examples`/<br>`noexamples` | _Enable/disable building G+Smo examples._ Default is enabled. |
| `lib`/<br>`nolib` | _Enable/disable building G+Smo as library._ Default is enabled. |
| `mpi`/<br>`nompi` | _Enable/disable building G+Smo with [MPI](https://www.mpi-forum.org) support._ Default is disabled. Not that Trilinos requires MPI to be enabled and will overwrite this setting. |
| `omp`/<br>`noomp` | _Enable/disable building G+Smo with [OpenMP](https://www.openmp.org) support._ Default is disabled for Clang compiler and enabled for all other compilers. |
| `pardiso`/<br>`nopardiso` | _Enable/disable building G+Smo with [Pardiso](https://www.pardiso-project.org) support._ Default is disabled. |
| `pastix`/<br>`nopastix` | _Enable/disable building G+Smo with [Pastix](https://solverstack.gitlabpages.inria.fr/pastix/index.html) support._ Default is disabled. |
| `pch`/<br>`nopch` | _Enable/disable building G+Smo with support for pre-compiled headers._ Default is disabled. |
| `psolid`/<br>`nopsolid` | _Enable/disable building G+Smo with [PSOLID](https://www.plm.automation.siemens.com/global/de/products/plm-components/parasolid.html) support._ Default is disabled. |
| `smesh`/<br>`nosmesh` | _Enable/disable building G+Smo with [Surface Mesh]() support._ Default is disabled. |
| `spectra`/<br>`nospectra` | _Enable/disable building G+Smo with [Spectra](https://spectralib.org) support._ Default is disabled. |
| `superlu`/<br>`nosuperlu` | _Enable/disable building G+Smo with [SuperLU](https://portal.nersc.gov/project/sparse/superlu/) support._ Default is disabled. |
| `taucs`/<br>`notaucs` | _Enable/disable building G+Smo with [Taucs](https://www.tau.ac.il/~stoledo/taucs/) support._ Default is disabled. |
| `umfpack`/<br>`noumfpack` | _Enable/disable building G+Smo with [UMFPACK](https://people.engr.tamu.edu/davis/suitesparse.html) support._ Default is disabled. |
| `unittests`/<br>`nounittests` | _Enable/disable building G+Smo unittests._ Default is disabled. |

| | Optional G+Smo extensions |
|-|-|
| `codipack`/<br>`nocodipack` | _Enable/disable [gsCoDiPack](https://github.com/SciCompKL/CoDiPack) extension._ Default is disabled. |
| `compflow`/<br>`nocompflow` | _Enable/disable [gsCompFlow](https://github.com/gismo/gsCompFlow) extension (**private repository!**)._ Default is disabled.
| `elasticity`/<br>`noelasticity` | _Enable/disable [gsElasticity](https://github.com/gismo/gsElasticity) extension._ Default is disabled.
| `exastencils`/<br>`noexastencils` | _Enable/disable [gsExastencils](https://github.com/gismo/gsExastencils) extension (**private repository!**)._ Default is disabled.
| `gmp`/<br>`nogmp` | _Enable/disable [gsGmp](https://gmplib.org) extension._ Default is disabled. This G+Smo extension is enabled automatically if `coeff_type` is set to `mpfr::mpreal` or `mpq_class`. |
|`ipopt`/<br>`noipopt` | _Enable/disable [gsIpOpt](https://github.com/coin-or/Ipopt) extension._ Default is disabled. |
| `klshell`/<br>`noklshell` | _Enable/disable [gsKLShell](https://github.com/gismo/gsKLShell) extension._ Default is disabled.
| `mpfr`/<br>`nompfr` | _Enable/disable [gsMpfr](https://www.mpfr.org) extension._ Default is disabled. This G+Smo extension is enabled automatically if `coeff_type` is set to `mpfr::mpreal`. |
| `motor`/<br>`nomotor` | _Enable/disable [motor](https://github.com/gismo/motor) extension (**private repository!**)._ Default is disabled.
| `opencascade`/<br>`noopencascade` | _Enable/disable [gsOpenCascade](https://www.opencascade.com) extension._ Default is disabled. |
| `opennurbs`/<br>`noopennurbs` | _Enable/disable [gsOpennurbs](https://www.rhino3d.com/opennurbs/) extension._ Default is disabled. |
| `trilinos`/<br>`notrilinos` | _Enable/disable [gsTrilinos](https://github.com/trilinos/Trilinos) extension._ Default is disabled. |
| `structuralanalysis`/<br>`nostructuralanalysis` | _Enable/disable [gsStructuralAnalysis](https://github.com/gismo/gsStructuralAnalysis) extension._ Default is disabled. |
| `universal`/<br>`nouniversal` | _Enable/disable [gsUniversal](https://github.com/stillwater-sc/universal) extension._ Default is disabled. This G+Smo extension is enabled automatically if `coeff_type` is set to any of the supported `posit_x_y` . |
| `unsupported`/<br>`nounsupported` | _Enable/disable [unsupported](https://github.com/gismo/unsupported) extension._ Default is disabled. |

__Disclaimer__: Not all G+Smo extensions are supported on all architectures and/or OS versions. Please report errors to the maintainer. Private repositories are not accessible without a password.

## Prebuilt Docker images
The following Docker images are built automatically for commits to the G+Smo stable branch. Unspecified options are set to their default values:

- `ubuntu20.04-gcc11-cxx11-release` ([GCC](https://gcc.gnu.org/) 11.x, C++11, Release mode with OpenMP) ![Docker Image Size (tag)](https://img.shields.io/docker/image-size/mmoelle1/gismo/ubuntu20.04-gcc11-cxx11-release) 
- `ubuntu20.04-clang11-cxx11-release` ([Clang](https://gcc.gnu.org/) 11.x, C++11, Release mode with OpenMP) ![Docker Image Size (tag)](https://img.shields.io/docker/image-size/mmoelle1/gismo/ubuntu20.04-clang11-cxx11-release)

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
docker run --rm -ti mmoelle1/gismo:ubuntu20.04-gcc11-cxx11-release
```
&nbsp;
You will be logged in under username `gismo` with standard password `gismo`. The `--rm` flag to the `docker run` command ensures that the container is removed upon termination. The directory `$HOME/gismo/build/bin` is added to the `PATH` variable so that example and tutorial applications can be run directly, e.g.:
```
docker run --rm -ti mmoelle1/gismo:ubuntu20.04-gcc11-cxx11-release geometry_example
```

&nbsp;
If you want to access data that is stored in a directory on your host computer from within your Docker container you need to mount (`-v`) the directory on your host computer when starting the Docker container. Assume that you want to make your current working directory on your host computer (`$(pwd)`) accessible under `/home/gismo/shared` in the docker container you need start the docker image as follows:
```
docker run --rm -ti -v $(pwd):/home/gismo/shared mmoelle1/gismo:ubuntu20.04-gcc11-cxx11-release
```
&nbsp;
It is also possible to mount a directory on your host computer to a directory of the Docker container and make it the working (`-w`) so that files generated inside the Docker container are accessible on your host computer:
```
docker run --rm -ti -w/home/gismo/shared -v $(pwd):/home/gismo/shared mmoelle1/gismo:ubuntu20.04-gcc11-cxx11-release
```
&nbsp;
A good starting point for learning the G+Smo library is to walk through the [example](https://www.gs.jku.at/trac/gismo/wiki/public/Doxygen/Examples) and [tutorial](https://www.gs.jku.at/trac/gismo/wiki/public/Doxygen/Tutorials) applications.

## For G+Smo developers

### Adding G+Smo extensions

Open the file `hooks/build` and add a code block of the form
```bash
# gsKLShell extension
shopt -s extglob;
case "$IMAGE_NAME" in
    ?(*\-)noklshell?(\-*))
        if [[ ${GISMO_SUBMODULES} =~ "gsKLShell" ]]; then
            GISMO_SUBMODULES="${GISMO_SUBMODULES//gsKLShell/}"
        fi
        ;;
    ?(*\-)klshell?(\-*))
        if [[ ${GISMO_SUBMODULES} ]]; then
            GISMO_SUBMODULES="${GISMO_SUBMODULES};gsKLShell"
        else
            GISMO_SUBMODULES="gsKLShell"
        fi
        ;;
esac
shopt -u extglob;
```

That's it for most extensions. 

If your extension requires certain options to be set make sure you add the above code snippet *before* the handling of these options, e.g., the **gsTrilinos** extension requires `GISMO_WITH_MPI=ON` and is therefore handled *before*. 
```bash
# gsTrilinos extension (must be before checking for MPI support)
shopt -s extglob;
case "$IMAGE_NAME" in
    ?(*\-)notrilinos?(\-*))
        if [[ ${GISMO_SUBMODULES} =~ "gsTrilinos" ]]; then
            GISMO_SUBMODULES="${GISMO_SUBMODULES//gsTrilinos/}"
        fi
        ;;
    ?(*\-)trilinos?(\-*))
        if [[ ${GISMO_SUBMODULES} ]]; then
            GISMO_SUBMODULES="${GISMO_SUBMODULES};gsTrilinos"
        else
            GISMO_SUBMODULES="gsTrilinos"
        fi
        GISMO_WITH_MPI=ON
        ;;
esac
shopt -u extglob;

[...]

# MPI support
shopt -s extglob;
case "$IMAGE_NAME" in
    ?(*\-)nompi?(\-*))
        GISMO_WITH_MPI=OFF
        ;;
    ?(*\-)mpi?(\-*))
        GISMO_WITH_MPI=ON
        ;;
    *)
        if [ -z "$GISMO_WITH_MPI" ]; then
            GISMO_WITH_MPI=OFF
        fi
        ;;
esac
shopt -u extglob;
```

If your extension requires additional software packages (e.g., a Fortran compiler or additional libraries and header files) to be installed in the Docker image open the file `ubuntu20.04/Dockerfile` and add a code block of the form
```bash
# Install prerequisites for gsOpenCascade extension
RUN if [ ${GISMO_SUBMODULES} =~ "gsOpenCascade" ] ; then \
    apt-get update -q && \
    apt-get install --no-install-recommends -yq \
        libgl-dev \
        libxi-dev \
        libxmu-dev \
        mesa-common-dev \
        tk-dev && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*; \
    fi
```
Repeat the above step for *all* `Dockerfile`s.

__Disclaimer__: Different OSes might have different naming conventions for libraries and/or header files and even use different commands for installing software packages. Consult your OS documentation for further details.

### Adding G+Smo options

The first step in adding a G+Smo option is similar to adding a G+Smo extension (see above). Open the file `hooks/build` and add a code block of the form
```bash
# MPI support
shopt -s extglob;
case "$IMAGE_NAME" in
    ?(*\-)nompi?(\-*))
        GISMO_WITH_MPI=OFF
        ;;
    ?(*\-)mpi?(\-*))
        GISMO_WITH_MPI=ON
        ;;
    *)
        if [ -z "$GISMO_WITH_MPI" ]; then
            GISMO_WITH_MPI=OFF
        fi
        ;;
esac
shopt -u extglob;
```

Next, open the file `ubuntu20.04/Dockerfile` and add a code block of the form
```bash
# GISMO_WITH_MPI          : {ON,OFF}
ARG GISMO_WITH_MPI=OFF
```
Additionally, add this flag to the `cmake` command at the bottom of the file
```bash
# Configure and build G+Smo library
RUN cd gismo    && \
    mkdir build && \
    cd build    && \
    CC=$CC \
    CXX=$CXX \
    FC=$FC \
    cmake .. \
        [...]
        -DGISMO_WITH_MPI=$GISMO_WITH_MPI \
        [...]
    make

```
Repeat the above step for all `Dockerfile`s.

### Adding Docker images

Adding a new Docker image '`newimage`' (replace with the desired name) requires some more work. 

First, open the file `run_build_docker.sh` and add a code block of the form
```bash
?(*:)?(*\-)newimage?(\-*))
        (DOCKER_REPO="mmoelle1/gismo"
         DOCKERFILE="${BASE_DIR}/newimage/Dockerfile"
         source ${BASE_DIR}/newimage/hooks/build)
        ;;
```

Next, create the directory `newimage` and create a `Dockerfile`. We suggest to copy an existing `Dockerfile`, e.g., from `ubuntu20.04` and update it to your needs. 

Next, create a new directory `newimage/hooks` and create a link named `build` to the file `hooks/build`. 

Finally, create a file `_compiler_settings` in the directory `newimage/hooks` with the following content (adapted to your needs):
```bash
#!/bin/bash
# hooks/_compiler_settings

# Compiler setting
shopt -s extglob;
case "$IMAGE_NAME" in

    # Clang
    
    ?(*\-)clang6.0?(\-*))
        CC=clang-6.0
        CXX=clang++-6.0
        ;;

    [...]

    # GCC

    ?(*\-)gcc7?(\-*))
        CC=gcc-7
        CXX=g++-7
        FC=gfortran-7
        ;;
    
    [...]
                      
    *)
        echo "[***] Unsupported compiler"
        exit 1
        ;;
esac
shopt -u extglob;

echo "[---] CC:                     ${CC}"
echo "[---] CXX:                    ${CXX}"
echo "[---] FC:                     ${FC}"
```

In summary, you should have created the following directory structure and files
```
newimage
|
+-- Dockerfile
|
+-- hooks
    |
    +-- _compiler_settings
    |
    +-- build -> ../../hooks/build
```