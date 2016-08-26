/** @file gsMpiHelper.h
    
    @brief Helpers for dealing with MPI.
 
    @ingroup ParallelCommunication
   
    Basically there are two helpers available:
    <dl>
    <dt>gsFakeMPIHelper</dt>
    <dd>A class adhering to the interface of gsMPIHelper
    that does not need MPI at all. This can be used
    to create a sequential program even if MPI is
    used to compile it.
    </dd>
    <dt>gsMPIHelper</dt>
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
    typedef gismo::gsMPIHelper gsMPIHelper;
    gsMPIHelper::instance(argc, argv);
    typename gsMPIHelper::MPICommunicator world =
    gsMPIHelper::getCommunicator();
    ...
    \endcode
   
    If one wants to have sequential program even if the code is
    compiled with mpi then one simply has to exchange the typedef
    with \code typedef gismo::gsMPIHelper gsFakeMPIHelper; \endcode.
   
    For checking whether we really use MPI or just fake please use
    gsMPIHelper::isFake (this is also possible at compile time!)


    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): C. Hofer -- taken from DUNE
    Created on: 2016-03-23
*/


#pragma once

#include <gsCore/gsConfig.h>
#include <gsCore/gsForwardDeclarations.h>
#include <gsCore/gsMemory.h>

#ifdef GISMO_WITH_MPI
#include <gsMpi/gsMpiTraits.h>
#include <gsMpi/gsBinaryFunctions.h>
#endif

#include <gsMpi/gsCollectiveCommunication.h>

namespace gismo
{

#ifndef GISMO_WITH_MPI
/**
 * @brief A fake mpi helper.
 *
 * This helper can be used if no MPI is available
 * or one wants to run sequentially even if MPI is
 * available and used.
 */
class gsFakeMPIHelper
{
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
    typedef No_Comm MPICommunicator;

    /** \brief get the default communicator
     *
     *  Return a communicator to exchange data with all processes
     *
     *  \returns a fake communicator
     */
    static MPICommunicator getCommunicator ()
    {
        static MPICommunicator comm;
        return comm;
    }

    /** \brief get a local communicator
     *
     *  Returns a communicator to communicate with the local process only
     *
     *  \returns a fake communicator
     */
    static MPICommunicator getLocalCommunicator ()
    {
        return getCommunicator();
    }



    static gsCollectiveCommunication<MPICommunicator>
    getCollectiveCommunication()
    {
        return gsCollectiveCommunication<MPICommunicator>(getCommunicator());
    }

    /**
     * @brief Get the singleton instance of the helper.
     *
     * This method has to be called with the same arguments
     * that the main method of the program was called:
     * \code
     * int main(int argc, char** argv){
     *   gsMPIHelper::instance(argc, argv);
     *   // program code comes here
     *   ...
     * }
     * \endcode
     * @param argc The number of arguments provided to main.
     * @param argv The arguments provided to main.
     */
    static gsFakeMPIHelper& instance(int argc, char** argv)
    {
        (void)argc; (void)argv;
        // create singleton instance
        static gsFakeMPIHelper singleton;
        return singleton;
    }

    static gsFakeMPIHelper& instance()
    {
        // create singleton instance
        static gsFakeMPIHelper singleton;
        return singleton;
    }

    /**
     * @brief return rank of process, i.e. zero
     */
    int rank () const { return 0; }
    /**
     * @brief return rank of process, i.e. one
     */
    int size () const { return 1; }

private:
    gsFakeMPIHelper() {}
    gsFakeMPIHelper(const gsFakeMPIHelper&);
    gsFakeMPIHelper& operator=(const gsFakeMPIHelper);
};
// We do not have MPI therefore gsFakeMPIHelper
// is the gsMPIHelper
/**
 * @brief If no MPI is available gsFakeMPIHelper becomes the gsMPIHelper
 * @ingroup ParallelCommunication
 */
typedef gsFakeMPIHelper gsMPIHelper;

#else

/**
 * @brief A real mpi helper.
 * @ingroup ParallelCommunication
 *
 * This helper should be used for parallel programs.
 */
class gsMPIHelper
{
public:
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
    typedef MPI_Comm MPICommunicator;

    /** \brief get the default communicator
     *
     *  Return a communicator to exchange data with all processes
     *
     *  \returns MPI_COMM_WORLD
     */
    static MPICommunicator getCommunicator ()
    {
        return MPI_COMM_WORLD;
    }

    /** \brief get a local communicator
     *
     *  Returns a communicator to exchange data with the local process only
     *
     *  \returns MPI_COMM_SELF
     */
    static MPICommunicator getLocalCommunicator ()
    {
        return MPI_COMM_SELF;
    }

    static gsCollectiveCommunication<MPICommunicator>
    getCollectiveCommunication()
    {
        return gsCollectiveCommunication<MPICommunicator>(getCommunicator());
    }
    /**
     * @brief Get the singleton instance of the helper.
     *
     * This method has to be called with the same arguments
     * that the main method of the program was called:
     * \code
     * int main(int argc, char** argv){
     *   gsMPIHelper::instance(argc, argv);
     *   // program code comes here
     *   ...
     * }
     * \endcode
     * @param argc The number of arguments provided to main.
     * @param argv The arguments provided to main.
     */
    static gsMPIHelper& instance(const int& argc, char**& argv)
    {
        // create singleton instance
        static gsMPIHelper singleton (argc, argv);
        return singleton;
    }

    static gsMPIHelper& instance()
    {
        // create singleton instance
        static gsMPIHelper singleton;
        return singleton;
    }

    /**
     * @brief return rank of process
     */
    int rank () const { return rank_; }
    /**
     * @brief return number of processes
     */
    int size () const { return size_; }

private:
    int rank_;
    int size_;
    static void prevent_warning(int){}
      
    gsMPIHelper()
    {
        init();
    }
      
    //! \brief calls MPI_Init with argc and argv as parameters
    gsMPIHelper(const int& argc, char**& argv)
    {
        init(&argc, argv);
    }

    void init(const int * argc = NULL, char** argv = NULL)
    {
        int wasInitialized = -1;
        MPI_Initialized( &wasInitialized );
        if(!wasInitialized)
        {
            rank_ = -1;
            size_ = -1;
            static int is_initialized = MPI_Init(const_cast<int*>(argc), &argv);
            prevent_warning(is_initialized);
        }

        MPI_Comm_rank(MPI_COMM_WORLD,&rank_);
        MPI_Comm_size(MPI_COMM_WORLD,&size_);

        GISMO_ASSERT( rank_ >= 0, "Invalid processor rank");
        GISMO_ASSERT( size_ >= 1, "Invalid size");

        gsDebug << "Called  MPI_Init on p=" << rank_ << "!" << std::endl;
    }
    
    //! \brief calls MPI_Finalize
    ~gsMPIHelper()
    {
        int wasFinalized = -1;
        MPI_Finalized( &wasFinalized );
        if(!wasFinalized) {
            MPI_Finalize();
            gsDebug << "Called MPI_Finalize on p=" << rank_ << "!" <<std::endl;
        }

    }
    gsMPIHelper(const gsMPIHelper&);
    gsMPIHelper& operator=(const gsMPIHelper);
};

#endif

}

