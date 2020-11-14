/** @file mpi_example.cpp

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

    Author(s): A. Mantzaflaris, C. Hofer, R. Schneckenleitner
*/

#include <gismo.h>

using namespace gismo;


// Parallel computation of an approximation of Pi
void approximatePI(const gsMpi & mpi, const gsMpiComm & comm);

// Parallel computation of a dot product
void computeDotProduct(const gsMpi & mpi, const gsMpiComm & comm);

// Subroutine to check whether the sampled position is valid, i.e., if
// the coordinate belongs to the quadrant then count is increased by 1.
void locateDot(const gsMatrix<> & loc, real_t & count);

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
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

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

    // Computes an approximation of PI
    approximatePI(mpi, comm);

    // Computes the dot product of two vectors
    computeDotProduct(mpi, comm);

    comm.barrier();
    gsInfo << "Good bye G+Smo, from process " << _rank << " on "
           << cpuname << "\n";

    return 0;
}

void approximatePI(const gsMpi & mpi, const gsMpiComm & comm)
{
    const int N  = 10000;
    real_t count = 0;

    //Get size and rank of the processor
    int _size = comm.size();
    int _rank = comm.rank();

    gsMatrix<> coord(N, 2);
    gsVector<> result(_size);

    comm.barrier();

    // Get current wall time
    double wtime = mpi.wallTime();

    // Define 0 as the master process
    if(_rank == 0)
    {
        // initialize the random number generator
        std::srand((unsigned)wtime);

        gsInfo<<"I am the master process of "<<_size<<" processes.\n";
        gsInfo << "I will finally compute an approximation of PI! " << "\n";

        // Create the N sample points for each processor
        for(index_t p = 0; p < _size; p++)
        {

            for (index_t i = 0; i < coord.rows(); i++)
            {
                coord(i, 0) = (real_t)(rand()) / RAND_MAX;
                coord(i, 1) = (real_t)(rand()) / RAND_MAX;
            }

            // Send the N sample points to the other processors
            if(p > 0)
                comm.send(coord.data(), coord.size(), p, 0);
                //comm.send(coord.data(), coord.size(), p); // tag = 0
            else
                locateDot(coord, count);
        }

    }
        // If the process is not the master process then wait for the sample points and check if the position of the
        // received coordinates is valid
    else
    {
        comm.recv(coord.data(), coord.size(), 0, 0);
        locateDot(coord, count);
    }

    // The master process gathers all variables count from each other process and stores the values in the vector
    // result
    comm.gather(&count, result.data(), 1, 0);

    // The master process finally counts the values of each processor up and approximates Pi
    if(_rank == 0)
    {
        for(index_t i = 1; i< _size; i++)
            count += result(i);

        gsInfo << "PI is approximately: " << (4. * count) / (N * _size) << "\n";
        gsInfo << "Time needed for parallel computation: " << mpi.wallTime() - wtime << "\n";
    }

}

void computeDotProduct(const gsMpi & mpi, const gsMpiComm & comm)
{
    const int N = 1000000;
    real_t res = 0;

    //Get size and rank of the processor
    int _size = comm.size();
    int _rank = comm.rank();

    // Assumes that the global vector can be split into smaller vectors of the same length
    if(N % _size != 0)
    {
        gsInfo << "Assume that the vectors can be split equivalently \n";
        gsInfo << N << " is not divisible by " << _size << "\n";
        gsInfo << "Dot product will not be computed! \n";
        return;
    }

    gsMatrix<> loc(N / _size, 2);
    gsVector<> result(_size);

    comm.barrier();

    // Get current wall time
    double wtime = mpi.wallTime();

    // Define 0 as the master process
    if (_rank == 0)
    {
        gsInfo << "I am the master process of " << _size << " processes.\n";
        gsInfo << "I will finally compute the dot product of two vectors consisting of ones! " << "\n";

        // The columns of the matrix are the vectors for computing the inner product
        gsMatrix<> colVecs = gsMatrix<>::Ones(N, 2);

        // Split the global vectors equivalently into smaller local vectors
        for (index_t p = 0; p < _size; p++)
        {

            for (index_t i = 0; i < loc.rows(); i++)
            {
                loc(i, 0) = colVecs(loc.rows() * p + i, 0);
                loc(i, 1) = colVecs(loc.rows() * p + i, 1);
            }

            // Distribute the smaller vectors to the other processes
            if (p > 0)
                comm.send(loc.data(), loc.size(), p, 0);
            else
                res = loc.col(0).dot(loc.col(1));
        }
    }
    // The other processes wait for the local vectors and compute the inner product of them
    else
    {
        comm.recv(loc.data(), loc.size(), 0, 0);
        res = loc.col(0).dot(loc.col(1));
    }

    // The inner products of the local vectors are collected and stored in result
    comm.gather(&res, result.data(), 1, 0);

    // The master process finally computes the inner product of the global vector
    if (_rank == 0)
    {
        for (index_t i = 1; i < _size; i++)
            res += result(i);

        gsInfo << "The result is: " << res << "\n";
        gsInfo << "Time needed for parallel computation: " << mpi.wallTime() - wtime << "\n";
    }

}

void locateDot(const gsMatrix<> & loc, real_t & count)
{
    count = 0;

    for(index_t row = 0; row < loc.rows(); row++)
    {
        const real_t pos = loc.row(row).squaredNorm();

        if(pos <= 1.)
            count += 1;
    }
}
