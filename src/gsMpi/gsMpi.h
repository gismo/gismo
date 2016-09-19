/** @file gsMpi.h
    
    @brief Helpers for dealing with MPI.
    
    Basically there are two helpers available:
    <dl>
    <dt>gsFakeMpicomm</dt>
    <dd>A class adhering to the interface of gsMpicomm
    that does not need MPI at all. This can be used
    to create a sequential program even if MPI is
    used to compile it.
    </dd>
    <dt>gsMpicomm</dt>
    <dd>A real MPI helper. When the singleton
    gets instantiated MPI_Init will be
    called and before the program exits
    MPI_Finalize will be called.
    </dd>
    </dl>
   
    Example of who to use these classes:
   
    A program that is parallel if compiled with MPI
    and sequential otherwise:
    \code
    int main(int argc, char** argv){
    typedef gismo::gsMpicomm gsMpicomm;
    gsMpicomm::instance(argc, argv);
    typename gsMpicomm::MPICommunicator world =
    gsMpicomm::getCommunicator();
    ...
    \endcode
   
    If one wants to have sequential program even if the code is
    compiled with mpi then one simply has to exchange the typedef
    with \code typedef gismo::gsMpicomm gsFakeMpicomm; \endcode.
   
    For checking whether we really use MPI or just fake please use
    gsMpicomm::isFake (this is also possible at compile time!)


    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): C. Hofer -- Code based on ideas from the DUNE library
*/


#pragma once

#include <gsCore/gsForwardDeclarations.h>

#ifdef GISMO_WITH_MPI
#include <mpi.h>
#include <gsMpi/gsMpiTraits.h>
#include <gsMpi/gsBinaryFunctions.h>
#endif

#include <gsMpi/gsMpiComm.h>


namespace gismo
{

class gsNoMpi;

#ifndef GISMO_WITH_MPI
/**
 * @brief If no MPI is available gsNoMpi becomes the gsMpi
 * @ingroup Mpi
 */
typedef gsNoMpi gsMpi;
 #else
class gsMpi;
#endif

/// Singleton function returning the gsMpi helper object
GISMO_EXPORT gsMpi & gsMpiSingleton(const int& argc = 0, char** argv = 0);

/// Singleton function returning the gsNoMpi helper object
GISMO_EXPORT gsNoMpi & gsNoMpiSingleton(const int& argc = 0, char** argv = 0);

/**
 * @brief A helper object for use when no MPI is present
 *
 * This helper can be used if no MPI is available or one wants to run
 * sequentially even if MPI is available and used.
 */
class gsNoMpi
{
    friend GISMO_EXPORT gsMpi   & gsMpiSingleton  (const int& argc, char** argv);
    friend GISMO_EXPORT gsNoMpi & gsNoMpiSingleton(const int& argc, char** argv);
    
public:
    enum {
        /**
         * @brief Are we fake (i.e. pretend to have MPI support but are compiled
         * without.)
         */
        isFake = true
    };

    /**
     * @brief The type of the mpi communicator.
     */
    typedef Serial_Comm Communicator;

    /** \brief get the default communicator
     *
     *  Return a communicator to exchange data with all processes
     *
     *  \returns a fake communicator
     */
    static Communicator worldComm()
    {
        return Serial_Comm();
    }

    /** \brief get a local communicator
     *
     *  Returns a communicator to communicate with the local process only
     *
     *  \returns a fake communicator
     */
    static Communicator localComm()
    {
        return worldComm();
    }
    
    /**
     * @brief Get the singleton instance of gsNoMpi
     *
     * This method has to be called with the same arguments
     * that the main method of the program was called:
     * \code
     * int main(int argc, char** argv){
     *   gsMpiComm::instance(argc, argv);
     *   // program code comes here
     *   ...
     * }
     * \endcode
     * @param argc The number of arguments provided to main.
     * @param argv The arguments provided to main.
     */
    static gsNoMpi& instance(int argc = 0, char** argv = NULL)
    {
        return gsNoMpiSingleton(argc,argv);
    }

    /**
     * @brief return rank of process, i.e. zero
     */
    static int worldRank () { return 0; }
    /**
     * @brief return rank of process, i.e. one
     */
    static int worldSize () { return 1; }

private:
    gsNoMpi();
    gsNoMpi(const gsNoMpi&);
    gsNoMpi(const int& argc, char** argv)
    {
        (void)argc;
        (void)argv;
    }

    gsNoMpi& operator=(const gsNoMpi&);

};

#ifdef GISMO_WITH_MPI

/**
 * @brief A parallel communicator class based on MPI
 *
 * @ingroup Mpi
 *
 */
class gsMpi
{
public:
      
    friend GISMO_EXPORT gsMpi & gsMpiSingleton(const int& argc, char** argv);
    
    enum {
        /**
         * @brief Are we fake (i. e. pretend to have MPI support but are compiled
         * without.
         */
        isFake = false
    };

    /**
     * @brief The type of the mpi communicator.
     */
    typedef MPI_Comm Communicator;

    /** \brief get the default communicator
     *
     *  Return a communicator to exchange data with all processes
     *
     *  \returns MPI_COMM_WORLD
     */
    static Communicator worldComm()
    {
        return MPI_COMM_WORLD;
    }

    /** \brief get a local communicator
     *
     *  Returns a communicator to exchange data with the local process only
     *
     *  \returns MPI_COMM_SELF
     */
    static Communicator localComm()
    {
        return MPI_COMM_SELF;
    }
    
    /**
     * @brief Returns the singleton instance of gsMpi
     *
     * This method has to be called with the same arguments
     * that the main method of the program was called:
     * \code
     * int main(int argc, char** argv){
     *   gsMpi::instance(argc, argv);
     *   // program code comes here
     *   ...
     * }
     * \endcode
     * @param argc The number of arguments provided to main.
     * @param argv The arguments provided to main.
     */
    static gsMpi& instance(const int& argc = 0, char** argv = NULL)
    {
        return gsMpiSingleton(argc,argv);
    }
    
    /**
     * @brief Returns the rank of process
     */
    static int worldRank () { return gsMpiComm(worldComm()).rank(); }
    /**
     * @brief Returns the number of processes
     */
    static int worldSize () { return gsMpiComm(worldComm()).size(); }

private:
    gsMpi();
    
    /// \brief calls MPI_Init with argc and argv as parameters
    gsMpi(const int& argc, char** argv)
    {
        init(const_cast<int*>(&argc), argv);
    }

    void init(int * argc = NULL, char** argv = NULL)
    {
        int initialized = -1;
        MPI_Initialized( &initialized );
        if( 0 == initialized )
        {
            //Note: valgrind false positive here, see
            // https://www.open-mpi.org/faq/?category=debugging#valgrind_clean
            initialized = MPI_Init(argc, &argv);
            GISMO_ENSURE(MPI_SUCCESS==initialized, "MPI failed to initialize");
        }
        //gsDebug << "Called  MPI_Init on p=" << rank_ << "!" << std::endl;
    }
    
    /// \brief calls MPI_Finalize
    ~gsMpi()
    {
        int wasFinalized = -1;
        MPI_Finalized( &wasFinalized );
        if( 0 == wasFinalized)
        {
            MPI_Finalize();
        }
        //gsDebug << "Called MPI_Finalize on p=" << rank_ << "!" <<std::endl;
    }
    
    gsMpi(const gsMpi&);
    gsMpi& operator=(const gsMpi&);

};

#endif

}

