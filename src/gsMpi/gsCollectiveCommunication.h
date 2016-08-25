/** @file gsCollectiveCommunication.h

    \brief Implements an utility class that provides
    collective communication methods for sequential programs.


    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): C.Hofer, taken from DUNE
    Created on: 2016-03-23
*/


#pragma once


namespace gismo
{

  /* define some type that definitely differs from MPI_Comm */
  struct No_Comm {};


  /*! @brief Collective communication interface and sequential default implementation

     gsCollectiveCommunication offers an abstraction to the basic methods
     of parallel communication, following the message-passing
     paradigm. It allows one to switch parallel features on and off, without
     changing the code. Currently only MPI and sequential code are
     supported.

     A gsCollectiveCommunication object is returned by all grids (also
     the sequential ones) in order to allow code to be written in
     a transparent way for sequential and parallel grids.

     This class provides a default implementation for sequential grids.
     The number of processes involved is 1, any sum, maximum, etc. returns
     just its input argument and so on.

     In specializations one can implement the real thing using appropriate
     communication functions, e.g. there exists an implementation using
     the Message Passing %Interface (MPI), see gismo::gsCollectiveCommunication<MPI_Comm>.

     \ingroup ParallelCommunication
   */
  template<typename C>
  class gsCollectiveCommunication
  {
  public:
    //! Construct default object
    gsCollectiveCommunication()
    {}
    gsCollectiveCommunication (const C&)
    {}

    //! Return rank, is between 0 and size()-1
    int rank () const
    {
      return 0;
    }

    //! Number of processes in set, is greater than 0
    int size () const
    {
      return 1;
    }

    /** @brief  Compute the sum of the argument over all processes and
            return the result in every process. Assumes that T has an operator+
     */
    template<typename T>
    T sum (T& in) const     // MPI does not know about const :-(
    {
      return in;
    }

    /** @brief Compute the sum over all processes for each component of an array and return the result
            in every process. Assumes that T has an operator+
     */
    template<typename T>
    int sum (T* inout, int len) const
    {
      return 0;
    }

    /** @brief  Compute the product of the argument over all processes and
            return the result in every process. Assumes that T has an operator*
     */
    template<typename T>
    T prod (T& in) const     // MPI does not know about const :-(
    {
      return in;
    }

    /** @brief Compute the product over all processes
            for each component of an array and return the result
            in every process. Assumes that T has an operator*
     */
    template<typename T>
    int prod (T* inout, int len) const
    {
      return 0;
    }

    /** @brief  Compute the minimum of the argument over all processes and
            return the result in every process. Assumes that T has an operator<
     */
    template<typename T>
    T min (T& in) const     // MPI does not know about const :-(
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

    /** @brief  Compute the maximum of the argument over all processes and
            return the result in every process. Assumes that T has an operator<
     */
    template<typename T>
    T max (T& in) const     // MPI does not know about const :-(
    {
      return in;
    }

    /** @brief Compute the maximum over all processes
            for each component of an array and return the result
            in every process. Assumes that T has an operator<
     */
    template<typename T>
    int max (T* inout, int len) const
    {
      return 0;
    }

    /** @brief Wait until all processes have arrived at this point in the program.
     */
    int barrier () const
    {
      return 0;
    }

    /** @brief Distribute an array from the process with rank root to all other processes
     */
    template<typename T>
    int broadcast (T* inout, int len, int root) const
    {
      return 0;
    }

    /** @brief  Gather arrays on root task.
     *
     * Each process sends its in array of length len to the root process
     * (including the root itself). In the root process these arrays are stored in rank
     * order in the out array which must have size len * number of processes.
     * @param[in] in The send buffer with the data to send.
     * @param[out] out The buffer to store the received data in. Might have length zero on non-root
     *                  tasks.
     * @param[in] len The number of elements to send on each task.
     * @param[in] root The root task that gathers the data.
     */
    template<typename T>
    int gather (T* in, T* out, int len, int root) const     // note out must have same size as in
    {
      for (int i=0; i<len; i++)
        out[i] = in[i];
      return 0;
    }

    /** @brief  Gather arrays of variable size on root task.
     *
     * Each process sends its in array of length sendlen to the root process
     * (including the root itself). In the root process these arrays are stored in rank
     * order in the out array.
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
     * The root process sends the elements with index from k*len to (k+1)*len-1 in its array to
     * task k, which stores it at index 0 to len-1.
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
     * The root process sends the elements with index from send+displ[k] to send+displ[k]-1 in
     * its array to task k, which stores it at index 0 to recvlen-1.
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

#ifdef HAVE_MPI

  //=======================================================
  // use singleton pattern and template specialization to
  // generate MPI operations
  //=======================================================


  /*! \brief Specialization of gsCollectiveCommunication for MPI
        \ingroup ParallelCommunication
   */
  template<>
  class gsCollectiveCommunication<MPI_Comm>
  {
  public:
    //! Instantiation using a MPI communicator
    gsCollectiveCommunication (const MPI_Comm& c = MPI_COMM_WORLD)
      : communicator(c)
    {
      if(communicator!=MPI_COMM_NULL) {
        int initialized = 0;
        MPI_Initialized(&initialized);
        if (!initialized)
          GISMO_ERROR("You must call gsMPIHelper::instance(argc,argv) in your main() function before using the MPI gsCollectiveCommunication!");
        MPI_Comm_rank(communicator,&me);
        MPI_Comm_size(communicator,&procs);
      }else{
        procs=0;
        me=-1;
      }
    }

    //! @copydoc gsCollectiveCommunication::rank
    int rank () const
    {
      return me;
    }

    //! @copydoc gsCollectiveCommunication::size
    int size () const
    {
      return procs;
    }

    //! @copydoc gsCollectiveCommunication::sum
    template<typename T>
    T sum (T& in) const     // MPI does not know about const :-(
    {
      T out;
      allreduce<std::plus<T> >(&in,&out,1);
      return out;
    }

    //! @copydoc gsCollectiveCommunication::sum
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


    //! @copydoc gsCollectiveCommunication::prod
    template<typename T>
    T prod (T& in) const     // MPI does not know about const :-(
    {
      T out;
      allreduce<std::multiplies<T> >(&in,&out,1);
      return out;
    }

    //! @copydoc gsCollectiveCommunication::prod
    template<typename T>
    int prod (T* inout, int len) const
    {
      return allreduce<std::multiplies<T> >(inout,len);
    }

    //! @copydoc gsCollectiveCommunication::min
    template<typename T>
    T min (T& in) const     // MPI does not know about const :-(
    {
      T out;
      allreduce<Min<T> >(&in,&out,1);
      return out;
    }

    //! @copydoc gsCollectiveCommunication::min
    template<typename T>
    int min (T* inout, int len) const
    {
      return allreduce<Min<T> >(inout,len);
    }


    //! @copydoc gsCollectiveCommunication::max
    template<typename T>
    T max (T& in) const     // MPI does not know about const :-(
    {
      T out;
      allreduce<Max<T> >(&in,&out,1);
      return out;
    }

    //! @copydoc gsCollectiveCommunication::max
    template<typename T>
    int max (T* inout, int len) const
    {
      return allreduce<Max<T> >(inout,len);
    }

    //! @copydoc gsCollectiveCommunication::barrier
    int barrier () const
    {
      return MPI_Barrier(communicator);
    }

    //! @copydoc gsCollectiveCommunication::broadcast
    template<typename T>
    int broadcast (T* inout, int len, int root) const
    {
      return MPI_Bcast(inout,len,MPITraits<T>::getType(),root,communicator);
    }

    //! @copydoc gsCollectiveCommunication::gather()
    //! @note out must have space for P*len elements
    template<typename T>
    int gather (T* in, T* out, int len, int root) const
    {
      return MPI_Gather(in,len,MPITraits<T>::getType(),
                        out,len,MPITraits<T>::getType(),
                        root,communicator);
    }

    //! @copydoc gsCollectiveCommunication::gatherv()
    template<typename T>
    int gatherv (T* in, int sendlen, T* out, int* recvlen, int* displ, int root) const
    {
      return MPI_Gatherv(in,sendlen,MPITraits<T>::getType(),
                         out,recvlen,displ,MPITraits<T>::getType(),
                         root,communicator);
    }

    //! @copydoc gsCollectiveCommunication::scatter()
    //! @note out must have space for P*len elements
    template<typename T>
    int scatter (T* send, T* recv, int len, int root) const
    {
      return MPI_Scatter(send,len,MPITraits<T>::getType(),
                         recv,len,MPITraits<T>::getType(),
                         root,communicator);
    }

    //! @copydoc gsCollectiveCommunication::scatterv()
    template<typename T>
    int scatterv (T* send, int* sendlen, int* displ, T* recv, int recvlen, int root) const
    {
      return MPI_Scatterv(send,sendlen,displ,MPITraits<T>::getType(),
                          recv,recvlen,MPITraits<T>::getType(),
                          root,communicator);
    }


    operator MPI_Comm () const
    {
      return communicator;
    }

    //! @copydoc gsCollectiveCommunication::allgather()
    template<typename T, typename T1>
    int allgather(T* sbuf, int count, T1* rbuf) const
    {
      return MPI_Allgather(sbuf, count, MPITraits<T>::getType(),
                           rbuf, count, MPITraits<T1>::getType(),
                           communicator);
    }

    //! @copydoc gsCollectiveCommunication::allgatherv()
    template<typename T>
    int allgatherv (T* in, int sendlen, T* out, int* recvlen, int* displ) const
    {
      return MPI_Allgatherv(in,sendlen,MPITraits<T>::getType(),
                            out,recvlen,displ,MPITraits<T>::getType(),
                            communicator);
    }

    //! @copydoc gsCollectiveCommunication::allreduce(Type* inout,int len) const
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
                                   (Generic_MPI_Op<Type, BinaryFunction>::get()),communicator);
    }


    //! @copydoc gsCollectiveCommunication::allreduce(Type* in,Type* out,int len) const
    template<typename BinaryFunction, typename Type>
    int allreduce(Type* in, Type* out, int len) const
    {
      return MPI_Allreduce(in, out, len, MPITraits<Type>::getType(),
                           (Generic_MPI_Op<Type, BinaryFunction>::get()),communicator);
    }



    //! @copydoc gsCollectiveCommunication::allreduce(Type* in,Type* out,int len) const
    template<typename BinaryFunction, typename Type>
    int iallreduce(Type* in, Type* out, int len,MPI_Request* req) const
    {
      return MPI_Iallreduce(in, out, len, MPITraits<Type>::getType(),
                           (Generic_MPI_Op<Type, BinaryFunction>::get()),communicator,req);
    }

    //! @copydoc gsCollectiveCommunication::allreduce(Type* inout,int len) const
    template<typename BinaryFunction, typename Type>
    int iallreduce(Type* inout, int len,MPI_Request* req) const
    {
      return MPI_Iallreduce(MPI_IN_PLACE, inout, len, MPITraits<Type>::getType(),
                                   (Generic_MPI_Op<Type, BinaryFunction>::get()),communicator,req);
    }

    template<typename BinaryFunction, typename Type>
    int reduce(Type* inout, int len,int root) const
    {
      int ret;
      if(root == rank())
        ret = MPI_Reduce(MPI_IN_PLACE, inout, len, MPITraits<Type>::getType(),
                       (Generic_MPI_Op<Type, BinaryFunction>::get()),root,communicator);
      else
        ret = MPI_Reduce(inout, NULL, len, MPITraits<Type>::getType(),
                         (Generic_MPI_Op<Type, BinaryFunction>::get()),root,communicator);
      return ret;
    }

    template<typename BinaryFunction, typename Type>
    int reduce(Type* in,Type* out, int len,int root) const
    {
      int ret;
      ret = MPI_Reduce(in, out, len, MPITraits<Type>::getType(),
                       (Generic_MPI_Op<Type, BinaryFunction>::get()),root,communicator);
      return ret;
    }

    template<typename BinaryFunction, typename Type>
    int ireduce(Type* inout, int len,int root, MPI_Request* req) const
    {
      int ret;
      if(root == rank())
        ret = MPI_Ireduce(MPI_IN_PLACE, inout, len, MPITraits<Type>::getType(),
                       (Generic_MPI_Op<Type, BinaryFunction>::get()),root,communicator,req);
      else
        ret = MPI_Ireduce(inout, inout, len, MPITraits<Type>::getType(),
                         (Generic_MPI_Op<Type, BinaryFunction>::get()),root,communicator,req);
      return ret;
    }

    template<typename BinaryFunction, typename Type>
    int ireduce(Type* in, Type* out, int len, int root, MPI_Request* req) const
    {
      return MPI_Ireduce(in, out, len, MPITraits<Type>::getType(),
                           (Generic_MPI_Op<Type, BinaryFunction>::get()),root,communicator,req);
    }

  private:
    MPI_Comm communicator;
    int me;
    int procs;
  };

#endif

}
