# XBraid extension

G+Smo extension for the [XBraid - Parallel-in-time Solver Package](https://github.com/XBraid/xbraid).

|CMake flags|```-DGISMO_OPTIONAL="gsXBraid"```|
|--:|---|
|Required additional CMake flags|```-DGISMO_WITH_MPI=ON``` (recommended)<br>```-DGISMO_WITH_OPENMP=ON``` (optionally)|
|License|[MPL 2.0](https://www.mozilla.org/en-US/MPL/2.0/)|
|OS support|Linux, Windows, macOS|
|Status|completed|
|Developer|Matthias Möller|
|Maintainer|M.Moller@tudelft.nl|
|Last checked|05-05-2021|

***
__Table of content__
1. [Introduction](#introduction)
2. [Usage example](#usage_example)
***

__Introdution__

The XBraid extension builds on the open-source
[XBraid](https://github.com/XBraid/xbraid) package developed at [
Lawrence Livermore National
Laboratory](https://computation.llnl.gov/projects/parallel-time-integration-multigrid/),
and at collaborating [academic
institutions](https://github.com/XBraid/xbraid/wiki/Team). XBraid is a
non-intrusive, optimal-scaling parallel-in-time solver that builds on
multigrid reduction techniques (multigrid-reduction-in-time or MGRIT).

The XBraid extension provides a generic wrapper to XBraid's C++
interface that can be easily customized by deriving an application
from the class `gsXBraid<T>` and overriding some or all virtual methods:

```cpp
virtual braid_Int Access(braid_Vector, BraidAccessStatus&);
virtual braid_Int BufPack(braid_Vector, void*, BraidBufferStatus&);
virtual braid_Int BufSize(braid_Int*, BraidBufferStatus&);
virtual braid_Int BufUnpack(void*, braid_Vector*, BraidBufferStatus&);
virtual braid_Int Clone(braid_Vector, braid_Vector*);
virtual braid_Int Coarsen(braid_Vector, braid_Vector*, BraidCoarsenRefStatus&);
virtual braid_Int Free(braid_Vector);
virtual braid_Int Init(braid_Real, braid_Vector*);
virtual braid_Int Refine(braid_Vector, braid_Vector*, BraidCoarsenRefStatus&);
virtual braid_Int Residual(braid_Vector, braid_Vector, BraidStepStatus&);
virtual braid_Int SpatialNorm(braid_Vector, braid_Real*);
virtual braid_Int Step(braid_Vector, braid_Vector, braid_Vector, BraidStepStatus&);
virtual braid_Int Sum(braid_Real, braid_Vector, braid_Real, braid_Vector);
```

__Usage example__

The file ```xbraid_heatEquation_example.cpp``` illustrates the basic usage of the gsXBraid extension.

1.  Configuration and compilation (MPI-only mode)

    ```bash
    mkdir build
    cd build
    cmake .. -DGISMO_OPTIONAL="gsXBraid" -DGISMO_WITH_MPI=ON
    make xbraid_heatEquation_example -j4
    ```
    
2.  Execution (MPI-only mode)

    ```bash
    mpirun -np <NPROC> --hostfile <HOSTFILE> ./bin/xbraid_heatEquation_example -n 250 -r 6 -i 3
    ```

    This will solve the two-dimensional heat equation on a unit square
    with 250 time steps in the time interval [0, 0.1] using <NPROC>
    MPI processes. The `hostfile` should have the following structure
    
    ```text
    node0 slots=#slots max_slots=#maximum slots
    node1 slots=#slots max_slots=#maximum slots
    ...
    ```
    
    The spatial domain is 6 times regularly refined in space (h-refinement) 
    and the approximation order is increased 3 times (p-refinement). 
    Order elevation instead of order increase can be achieved by replacing 
    the switch`-i` by `-e`.

    For a complete list of command-line argument run
    ```bash
    ./bin/xbraid_heatEquation_example -h
    ```
    
3.  Configuration and compilation (MPI-OpenMP mode)

    ```bash
    mkdir build
    cd build
    cmake .. -DGISMO_OPTIONAL="gsXBraid" -DGISMO_WITH_MPI=ON -DGISMO_WITH_OPENMP=ON
    make xbraid_heatEquation_example -j4
    ```
    
4.  Execution (MPI-OpenMP mode)

    ```bash
    mpirun -np <NPROC> --hostfile <HOSTFILE> -x OMP_NUM_THREADS=<NTHREAD> ./bin/xbraid_heatEquation_example -n 250 -r 6 -i 3
    ```

    The additional parameter `-x OMP_NUM_THREADS=<NTHREAD>` ensures that
    each MPI process executes `NTHREAD` OpenMP threads in parallel. The `-x` 
    flag is not supported by all MPI implementations. If it does not work
    try
    
    ```bash
    mpirun -np <NPROC> --hostfile <HOSTFILE> -env OMP_NUM_THREADS <NTHREAD> ./bin/xbraid_heatEquation_example -n 250 -r 6 -i 3
    ```
