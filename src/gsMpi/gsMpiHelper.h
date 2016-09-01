/** @file gsMpiHelper.h
    
    @brief Helpers for dealing with MPI.
    
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

#define gsMPIHelper gsMpiComm

namespace gismo
{

class gsSerialComm;

#ifndef GISMO_WITH_MPI
/**
 * @brief If no MPI is available gsSerialComm becomes the gsMpiComm
 * @ingroup Mpi
 */
typedef gsSerialComm gsMpiComm;
#else
class gsMpiComm;
#endif

/// Singleton function returning the gsMPiHelper
GISMO_EXPORT gsMpiComm & gsMpiCommSingleton(const int& argc = 0, char** argv = 0);

/// Singleton function returning the gsSerialComm
GISMO_EXPORT gsSerialComm & gsSerialCommSingleton(const int& argc = 0, char** argv = 0);

/** \brief Dummy communication type for serial (ie. no) communication */
struct Serial_Comm {};

/**
 * @brief A serial communication class
 *
 * This helper can be used if no MPI is available or one wants to run
 * sequentially even if MPI is available and used.
 */
class gsSerialComm
{
    friend GISMO_EXPORT gsMpiComm    & gsMpiCommSingleton   (const int& argc, char** argv);
    friend GISMO_EXPORT gsSerialComm & gsSerialCommSingleton(const int& argc, char** argv);
    
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
     * @brief Get the singleton instance of gsSerialComm
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
    static gsSerialComm& instance(int argc = 0, char** argv = NULL)
    {
        return gsSerialCommSingleton(argc,argv);
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
    gsSerialComm();
    gsSerialComm(const gsSerialComm&);
    gsSerialComm(const int& argc, char** argv)
    {
        (void)argc;
        (void)argv;
    }

    gsSerialComm& operator=(const gsSerialComm&);

public:

    /** @brief  Compute the sum of the argument over all processes and
        return the result in every process. Assumes that T has an operator+
    */
    template<typename T>
    T sum (T& in) const
    {
        return in;
    }

    /** @brief Compute the sum over all processes for each component
        of an array and return the result in every process. Assumes
        that T has an operator+
    */
    template<typename T>
    int sum (T* inout, int len) const
    {
        return 0;
    }

    /** @brief Compute the product of the argument over all processes
        and return the result in every process. Assumes that T has an
        operator*
    */
    template<typename T>
    T prod (T& in) const
    {
        return in;
    }

    /** @brief Compute the product over all processes for each
        component of an array and return the result in every
        process. Assumes that T has an operator*
    */
    template<typename T>
    int prod (T* inout, int len) const
    {
        return 0;
    }

    /** @brief Compute the minimum of the argument over all processes
        and return the result in every process. Assumes that T has an
        operator<
    */
    template<typename T>
    T min (T& in) const
    {
        return in;
    }

    /** @brief Compute the minimum over all processes
        for each component of an array and return the result
        in every process. Assumes that T has an operator<
    */
    template<typename T>
    int min (T* inout, int len) const
    {
        return 0;
    }

    /** @brief Compute the maximum of the argument over all processes
        and return the result in every process. Assumes that T has an
        operator<
    */
    template<typename T>
    T max (T& in) const
    {
        return in;
    }

    /** @brief Compute the maximum over all processes for each
        component of an array and return the result in every
        process. Assumes that T has an operator<
    */
    template<typename T>
    int max (T* inout, int len) const
    {
        return 0;
    }

    /** @brief Wait until all processes have arrived at this point in
     * the program.
     */
    int barrier () const
    {
        return 0;
    }

    /** @brief Distribute an array from the process with rank root to
     * all other processes
     */
    template<typename T>
    int broadcast (T* inout, int len, int root) const
    {
        return 0;
    }

    /** @brief  Gather arrays on root task.
     *
     * Each process sends its in array of length len to the root
     * process (including the root itself). In the root process these
     * arrays are stored in rank order in the out array which must
     * have size len * number of processes.  @param[in] in The send
     * buffer with the data to send.  @param[out] out The buffer to
     * store the received data in. Might have length zero on non-root
     * tasks.  @param[in] len The number of elements to send on each
     * task.  @param[in] root The root task that gathers the data.
     */
    template<typename T>
    int gather (T* in, T* out, int len, int root) const // note out must have same size as in
    {
        // copy_n(in, len, out);
        for (int i=0; i<len; i++)
            out[i] = in[i];
        return 0;
    }

    /** @brief  Gather arrays of variable size on root task.
     *
     * Each process sends its in array of length sendlen to the root process
     * (including the root itself). In the root process these arrays are stored in rank
     * order in the out array.
     *
     * @param[in] in The send buffer with the data to be sent
     * @param[in] sendlen The number of elements to send on each task
     * @param[out] out The buffer to store the received data in. May have length zero on non-root
     *                 tasks.
     * @param[in] recvlen An array with size equal to the number of processes containing the number
     *                    of elements to receive from process i at position i, i.e. the number that
     *                    is passed as sendlen argument to this function in process i.
     *                    May have length zero on non-root tasks.
     * @param[out] displ An array with size equal to the number of processes. Data received from
     *                  process i will be written starting at out+displ[i] on the root process.
     *                  May have length zero on non-root tasks.
     * @param[in] root The root task that gathers the data.
     */
    template<typename T>
    int gatherv (T* in, int sendlen, T* out, int* recvlen, int* displ, int root) const
    {
        for (int i=*displ; i<sendlen; i++)
            out[i] = in[i];
        return 0;
    }

    /** @brief Scatter array from a root to all other task.
     *
     * The root process sends the elements with index from k*len to
     * (k+1)*len-1 in its array to task k, which stores it at index 0
     * to len-1.
     *
     * @param[in] send The array to scatter. Might have length zero on non-root
     *                  tasks.
     * @param[out] recv The buffer to store the received data in. Upon completion of the
     *                 method each task will have same data stored there as the one in
     *                 send buffer of the root task before.
     * @param[in] len The number of elements in the recv buffer.
     * @param[in] root The root task that gathers the data.
     */
    template<typename T>
    int scatter (T* send, T* recv, int len, int root) const // note out must have same size as in
    {
        for (int i=0; i<len; i++)
            recv[i] = send[i];
        return 0;
    }

    /** @brief Scatter arrays of variable length from a root to all other tasks.
     *
     * The root process sends the elements with index from
     * send+displ[k] to send+displ[k]-1 in * its array to task k,
     * which stores it at index 0 to recvlen-1.
     *
     * @param[in] send The array to scatter. May have length zero on non-root
     *                  tasks.
     * @param[in] sendlen An array with size equal to the number of processes containing the number
     *                    of elements to scatter to process i at position i, i.e. the number that
     *                    is passed as recvlen argument to this function in process i.
     * @param[in] displ An array with size equal to the number of processes. Data scattered to
     *                  process i will be read starting at send+displ[i] on root the process.
     * @param[out] recv The buffer to store the received data in. Upon completion of the
     *                  method each task will have the same data stored there as the one in
     *                  send buffer of the root task before.
     * @param[in] recvlen The number of elements in the recv buffer.
     * @param[in] root The root task that gathers the data.
     */
    template<typename T>
    int scatterv (T* send, int* sendlen, int* displ, T* recv, int recvlen, int root) const
    {
        for (int i=*displ; i<*sendlen; i++)
            recv[i] = send[i];
        return 0;
    }

    /**
     * @brief Gathers data from all tasks and distribute it to all.
     *
     * The block of data sent from the  jth  process  is  received  by  every
     *  process and placed in the jth block of the buffer recvbuf.
     *
     * @param[in] sbuf The buffer with the data to send. Has to be the same for
     *                 each task.
     * @param[in] count The number of elements to send by any process.
     * @param[out] rbuf The receive buffer for the data. Has to be of size
     *  notasks*count, with notasks being the number of tasks in the communicator.
     */
    template<typename T>
    int allgather(T* sbuf, int count, T* rbuf) const
    {
        for(T* end=sbuf+count; sbuf < end; ++sbuf, ++rbuf)
            *rbuf=*sbuf;
        return 0;
    }

    /**
     * @brief Gathers data of variable length from all tasks and distribute it to all.
     *
     * The block of data sent from the jth process is received by every
     *  process and placed in the jth block of the buffer out.
     *
     * @param[in] in The send buffer with the data to send.
     * @param[in] sendlen The number of elements to send on each task.
     * @param[out] out The buffer to store the received data in.
     * @param[in] recvlen An array with size equal to the number of processes containing the number
     *                    of elements to recieve from process i at position i, i.e. the number that
     *                    is passed as sendlen argument to this function in process i.
     * @param[in] displ An array with size equal to the number of processes. Data recieved from
     *                  process i will be written starting at out+displ[i].
     */
    template<typename T>
    int allgatherv (T* in, int sendlen, T* out, int* recvlen, int* displ) const
    {
        for (int i=*displ; i<sendlen; i++)
            out[i] = in[i];
        return 0;
    }

    /**
     * @brief Compute something over all processes
     * for each component of an array and return the result
     * in every process.
     *
     * The template parameter BinaryFunction is the type of
     * the binary function to use for the computation
     *
     * @param inout The array to compute on.
     * @param len The number of components in the array
     */
    template<typename BinaryFunction, typename Type>
    int allreduce(Type* inout, int len) const
    {
        return 0;
    }

    /**
     * @brief Compute something over all processes
     * for each component of an array and return the result
     * in every process.
     *
     * The template parameter BinaryFunction is the type of
     * the binary function to use for the computation
     *
     * @param in The array to compute on.
     * @param out The array to store the results in.
     * @param len The number of components in the array
     */
    template<typename BinaryFunction, typename Type>
    void allreduce(Type* in, Type* out, int len) const
    {
        std::copy(in, in+len, out);
        return;
    }
};

#ifdef GISMO_WITH_MPI

/**
 * @brief A parallel communicator class based on MPI
 *
 * @ingroup Mpi
 *
 */
class gsMpiComm
{
public:
    
    friend GISMO_EXPORT gsMpiComm & gsMpiCommSingleton(const int& argc, char** argv);
    
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
     * @brief Returns the singleton instance of gsMpiComm
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
    static gsMpiComm& instance(const int& argc = 0, char** argv = NULL)
    {
        return gsMpiCommSingleton(argc,argv);
    }
    
    /**
     * @brief Returns the rank of process
     */
    int rank () const { return rank_; }
    /**
     * @brief Returns the number of processes
     */
    int size () const { return size_; }

private:
    int rank_;
    int size_;

private:
    gsMpiComm();
    
    /// \brief calls MPI_Init with argc and argv as parameters
    gsMpiComm(const int& argc, char** argv)
    {
        init(const_cast<int*>(&argc), argv);
    }

    void init(int * argc = NULL, char** argv = NULL)
    {
        int initialized = -1;
        MPI_Initialized( &initialized );
        if( 0 == initialized )
        {
            rank_ = -1;
            size_ = -1;
            //Note: valgrind false positive here, see
            // https://www.open-mpi.org/faq/?category=debugging#valgrind_clean
            initialized = MPI_Init(argc, &argv);
            GISMO_ENSURE(MPI_SUCCESS==initialized, "MPI failed to initialize");
        }

        MPI_Comm_rank(MPI_COMM_WORLD,&rank_);
        MPI_Comm_size(MPI_COMM_WORLD,&size_);

        GISMO_ENSURE( rank_ >= 0, "Invalid processor rank");
        GISMO_ENSURE( size_ >= 1, "Invalid size");

        //gsDebug << "Called  MPI_Init on p=" << rank_ << "!" << std::endl;
    }
    
    /// \brief calls MPI_Finalize
    ~gsMpiComm()
    {
        int wasFinalized = -1;
        MPI_Finalized( &wasFinalized );
        if( 0 == wasFinalized)
        {
            MPI_Finalize();
        }
        //gsDebug << "Called MPI_Finalize on p=" << rank_ << "!" <<std::endl;
    }
    
    gsMpiComm(const gsMpiComm&);
    gsMpiComm& operator=(const gsMpiComm&);

public:

    /// @copydoc gsSerialComm::sum
    template<typename T>
    T sum (T& in) const
    {
        T out;
        allreduce<std::plus<T> >(&in,&out,1);
        return out;
    }

    /// @copydoc gsSerialComm::sum
    template<typename T>
    int sum (T* inout, int len) const
    {
        return allreduce<std::plus<T> >(inout,len);
    }

    template<typename T>
    int sum (T* inout, int len, int root) const
    {
        return reduce<std::plus<T> >(inout,len,root);
    }

    template<typename T>
    int sum (T* in, T* out, int len, int root) const
    {
        return reduce<std::plus<T> >(in,out,len,root);
    }

    template<typename T>
    int isum (T* in,T* out, int len, int root, MPI_Request* req) const
    {
        return iallreduce<std::plus<T> >(in, out,len,req);
    }

    template<typename T>
    int isum (T* inout, int len, int root, MPI_Request* req) const
    {
        return ireduce<std::plus<T> >(inout,len,root,req);
    }

    template<typename T>
    int isum (T* inout, int len, MPI_Request* req) const
    {
        return iallreduce<std::plus<T> >(inout,len,req);
    }


    /// @copydoc gsSerialComm::prod
    template<typename T>
    T prod (T& in) const
    {
        T out;
        allreduce<std::multiplies<T> >(&in,&out,1);
        return out;
    }

    /// @copydoc gsSerialComm::prod
    template<typename T>
    int prod (T* inout, int len) const
    {
        return allreduce<std::multiplies<T> >(inout,len);
    }

    /// @copydoc gsSerialComm::min
    template<typename T>
    T min (T& in) const
    {
        T out;
        allreduce<Min<T> >(&in,&out,1);
        return out;
    }

    /// @copydoc gsSerialComm::min
    template<typename T>
    int min (T* inout, int len) const
    {
        return allreduce<Min<T> >(inout,len);
    }


    /// @copydoc gsSerialComm::max
    template<typename T>
    T max (T& in) const
    {
        T out;
        allreduce<Max<T> >(&in,&out,1);
        return out;
    }

    /// @copydoc gsSerialComm::max
    template<typename T>
    int max (T* inout, int len) const
    {
        return allreduce<Max<T> >(inout,len);
    }

    /// @copydoc gsSerialComm::barrier
    int barrier () const
    {
        return MPI_Barrier(worldComm());
    }

    /// @copydoc gsSerialComm::broadcast
    template<typename T>
    int broadcast (T* inout, int len, int root) const
    {
        return MPI_Bcast(inout,len,MPITraits<T>::getType(),root,worldComm());
    }

    /// @copydoc gsSerialComm::gather()
    /// @note out must have space for P*len elements
    template<typename T>
    int gather (T* in, T* out, int len, int root) const
    {
        return MPI_Gather(in,len,MPITraits<T>::getType(),
                          out,len,MPITraits<T>::getType(),
                          root,worldComm());
    }

    /// @copydoc gsSerialComm::gatherv()
    template<typename T>
    int gatherv (T* in, int sendlen, T* out, int* recvlen, int* displ, int root) const
    {
        return MPI_Gatherv(in,sendlen,MPITraits<T>::getType(),
                           out,recvlen,displ,MPITraits<T>::getType(),
                           root,worldComm());
    }

    /// @copydoc gsSerialComm::scatter()
    /// @note out must have space for P*len elements
    template<typename T>
    int scatter (T* send, T* recv, int len, int root) const
    {
        return MPI_Scatter(send,len,MPITraits<T>::getType(),
                           recv,len,MPITraits<T>::getType(),
                           root,worldComm());
    }

    /// @copydoc gsSerialComm::scatterv()
    template<typename T>
    int scatterv (T* send, int* sendlen, int* displ, T* recv, int recvlen, int root) const
    {
        return MPI_Scatterv(send,sendlen,displ,MPITraits<T>::getType(),
                            recv,recvlen,MPITraits<T>::getType(),
                            root,worldComm());
    }


    operator MPI_Comm () const
    {
        return worldComm();
    }

    /// @copydoc gsSerialComm::allgather()
    template<typename T, typename T1>
    int allgather(T* sbuf, int count, T1* rbuf) const
    {
        return MPI_Allgather(sbuf, count, MPITraits<T>::getType(),
                             rbuf, count, MPITraits<T1>::getType(),
                             worldComm());
    }

    /// @copydoc gsSerialComm::allgatherv()
    template<typename T>
    int allgatherv (T* in, int sendlen, T* out, int* recvlen, int* displ) const
    {
        return MPI_Allgatherv(in,sendlen,MPITraits<T>::getType(),
                              out,recvlen,displ,MPITraits<T>::getType(),
                              worldComm());
    }

    /// @copydoc gsSerialComm::allreduce(Type* inout,int len) const
    template<typename BinaryFunction, typename Type>
    int allreduce(Type* inout, int len) const
    {
        /*
          Type* out = new Type[len];
          int ret = allreduce<BinaryFunction>(inout,out,len);
          std::copy(out, out+len, inout);
          delete[] out;
          return ret;
        */
        return MPI_Allreduce(MPI_IN_PLACE, inout, len, MPITraits<Type>::getType(),
                             (Generic_MPI_Op<Type, BinaryFunction>::get()),worldComm());
    }


    /// @copydoc gsSerialComm::allreduce(Type* in,Type* out,int len) const
    template<typename BinaryFunction, typename Type>
    int allreduce(Type* in, Type* out, int len) const
    {
        return MPI_Allreduce(in, out, len, MPITraits<Type>::getType(),
                             (Generic_MPI_Op<Type, BinaryFunction>::get()),worldComm());
    }



    /// @copydoc gsSerialComm::allreduce(Type* in,Type* out,int len) const
    template<typename BinaryFunction, typename Type>
    int iallreduce(Type* in, Type* out, int len,MPI_Request* req) const
    {
        return MPI_Iallreduce(in, out, len, MPITraits<Type>::getType(),
                              (Generic_MPI_Op<Type, BinaryFunction>::get()),worldComm(),req);
    }

    /// @copydoc gsSerialComm::allreduce(Type* inout,int len) const
    template<typename BinaryFunction, typename Type>
    int iallreduce(Type* inout, int len,MPI_Request* req) const
    {
        return MPI_Iallreduce(MPI_IN_PLACE, inout, len, MPITraits<Type>::getType(),
                              (Generic_MPI_Op<Type, BinaryFunction>::get()),worldComm(),req);
    }

    template<typename BinaryFunction, typename Type>
    int reduce(Type* inout, int len,int root) const
    {
        int ret;
        if(root == rank())
            ret = MPI_Reduce(MPI_IN_PLACE, inout, len, MPITraits<Type>::getType(),
                             (Generic_MPI_Op<Type, BinaryFunction>::get()),root,worldComm());
        else
            ret = MPI_Reduce(inout, NULL, len, MPITraits<Type>::getType(),
                             (Generic_MPI_Op<Type, BinaryFunction>::get()),root,worldComm());
        return ret;
    }

    template<typename BinaryFunction, typename Type>
    int reduce(Type* in,Type* out, int len,int root) const
    {
        int ret;
        ret = MPI_Reduce(in, out, len, MPITraits<Type>::getType(),
                         (Generic_MPI_Op<Type, BinaryFunction>::get()),root,worldComm());
        return ret;
    }

    template<typename BinaryFunction, typename Type>
    int ireduce(Type* inout, int len,int root, MPI_Request* req) const
    {
        int ret;
        if(root == rank())
            ret = MPI_Ireduce(MPI_IN_PLACE, inout, len, MPITraits<Type>::getType(),
                              (Generic_MPI_Op<Type, BinaryFunction>::get()),root,worldComm(),req);
        else
            ret = MPI_Ireduce(inout, inout, len, MPITraits<Type>::getType(),
                              (Generic_MPI_Op<Type, BinaryFunction>::get()),root,worldComm(),req);
        return ret;
    }

    template<typename BinaryFunction, typename Type>
    int ireduce(Type* in, Type* out, int len, int root, MPI_Request* req) const
    {
        return MPI_Ireduce(in, out, len, MPITraits<Type>::getType(),
                           (Generic_MPI_Op<Type, BinaryFunction>::get()),root,worldComm(),req);
    }

};

#endif

}

