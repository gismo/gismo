/** @file gsMpiTest.cpp

    @brief Testing MPI with G+Smo

    Execute (eg. with 10 processes):
       mpirun -np 10 ./bin/mpi_example

    or provide a hosts file on a cluster:
       mpirun -hostfile hosts.txt ./bin/mpi_example

    If your cluster is using srun:
       srun -N 10 ./bin/mpi_example

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris, C. Hofer
*/

#include <gismo.h>

using namespace gismo;


int main(int argc, char **argv)
{
    gsCmdLine cmd("An example for testing MPI with G+Smo.\n"
        "Execute (eg. with 10 processes):                                      "
        "  *  mpirun -np 10 ./bin/mpi_example\n"
        "or provide a hosts file on a cluster:                                 "
        "  *  mpirun -hostfile hosts.txt ./bin/mpi_example\n"
        "If your cluster is using srun:                                        "
        "  *  srun -N 10 ./bin/mpi_example"
    );
    cmd.getValues(argc,argv);

    // Conditional compilation
#ifdef GISMO_WITH_MPI
    gsInfo << "Gismo was compiled with MPI support.\n";
#else
    gsInfo << "Gismo was compiled without MPI support.\n";
#endif

    // Initialize the MPI environment
    const gsMpi & mpi = gsMpi::init(argc, argv);

    // Get current wall time
    double wtime = mpi.wallTime();

    // Get the world communicator
    gsMpiComm comm = mpi.worldComm();

    //Get size and rank of the processor
    int _size = comm.size();
    int _rank = comm.rank();

    if (0==_rank)
        gsInfo<<"Running on "<<_size<<" processes.\n";
    comm.barrier();

    gsInfo <<"MPI is "<< (mpi.initialized() ? "" : "NOT ")
           <<"initialized on process "<< _rank <<"\n";
    comm.barrier();

    std::string cpuname = mpi.getProcessorName();

    // Print off a hello world message
    gsInfo << "Hello G+Smo, from process " << _rank <<" on "
           << cpuname <<", elapsed time is "<< mpi.wallTime()-wtime<< "\n";

    return 0;
}
