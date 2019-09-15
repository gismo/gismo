/**
 * @file gsMpiDynamicTraits.h
 *
 * @brief Specialized sender and receiver classes for dynamic sized types.
 */

#pragma once

#include <utility>
#include <cstddef>
#include "gsMpiComm.h"

namespace gismo
{
    struct gsActionableMpiRequest {
        gsMpiRequest req;
        std::function<int(MPI_Status*)> action;
    };


    // template<typename T, int _Rows, int _Options>
    // struct MPI_Dynamic_req {
    //     bool completed;
    //     gsVector<T, _Rows, _Options>* buffer;
    //     int* length;
    //     int source;
    //     MPI_Request req;
    // };

  /**
   * @brief Specialized sender and receiver classes for dynamic sized types.
   *
   * Specializations should provide static methods
   * \code
   *  static int send (T* in, int len, int dest, int tag = 0) const
   * \endcode
   */
    template<typename T, int _Rows, int _Options>
    static int dsend(gsMpiComm &comm, gsVector<T, _Rows, _Options>* in, int len, int dest, int tag = 0)
    {
        int length = in->size();
        comm.send(&length, 1, dest, tag);
        return comm.send(in->begin(), length, dest, tag);
    }

    /**
     * @brief Specialized dynamic non-blocking send for gsVector of any size.
     *
     * reqs should be an array that can hold at least two MPI_Requests.
     */
    template<typename T, int _Rows, int _Options>
    static int disend(gsMpiComm &comm, gsVector<T, _Rows, _Options>* in, int len, int dest, MPI_Request* reqs, int tag = 0)
    {
        int length = in->size();
        comm.isend(&length, 1, dest, &reqs[0], tag);

        return comm.isend(in->begin(), length, dest, &reqs[1], tag);
    }

    /**
     * @brief Specialized dynamic receive for gsVector of any size.
     *
     * statuses should be an array that can hold at least two MPI_Requests.
     */
    template<typename T, int _Rows, int _Options>
    static int drecv(gsMpiComm &comm, gsVector<T, _Rows, _Options>* out, int len, int source, int tag = 0, MPI_Status* statuses = NULL)
    {
        int length;
        comm.recv(&length, 1, source, tag, (statuses == NULL ? NULL : &statuses[0]));

        if(out->size() < length) {
            out->conservativeResize(length);
        }

        return comm.recv(out->begin(), length, source, tag, (statuses == NULL ? NULL : &statuses[1]));
    }

    template<typename T, int _Rows, int _Options>
    static gsActionableMpiRequest direcv(gsMpiComm &comm, gsVector<T, _Rows, _Options>* out, int len, int source, int tag = 0)
    {
        int* length = (int*) malloc(sizeof(int));
        gsActionableMpiRequest req;
        comm.irecv(length, 1, source, &req.req, tag);
        req.action = [comm, out, length, source, tag](MPI_Status* status) {
            printf("Received vector of length: %d. Current length: %d\n", *length, out->size());
            if(out->size() < *length) {
                printf("Resizing...\n");
                out->conservativeResize(*length);
            }

            return comm.recv(out->begin(), *length, source, tag, status);
        };

        return req;
    }

    static gsMpiStatus dwaitAny(int numberRequests, gsActionableMpiRequest reqs[], int* outIndex)
    {
        gsMpiRequest mpiReqs[numberRequests];
        for(int i = 0; i < numberRequests; i++) {
            mpiReqs[i] = reqs[i].req;
        }

        gsMpiRequest::waitAny(numberRequests, mpiReqs, outIndex);

        gsActionableMpiRequest req = reqs[*outIndex];
        gsMpiStatus status;
        req.action(&status);
        return status;
    }

    // template<typename T, int _Rows, int _Options>
    // static MPI_Dynamic_req<T, _Rows, _Options> direcv(gsMpiComm &comm, gsVector<T, _Rows, _Options>* out, int len, int source, int tag = 0)
    // {
    //     MPI_Request req;
    //     int* length = (int*) malloc(sizeof(int));
    //     comm.irecv(length, 1, source, mpi_req, tag);

    //     MPI_Dynamic_req req;
    //     req.buffer = out;
    //     req.length = length;
    //     req.source = source;
    //     req.tag = tag;
    //     req.req = mpi_req;

    //     return req;
    // }

    // template<typename T, int _Rows, int _Options>
    // static int dwaitAny(gsMpiComm &comm, MPI_Dynamic_req<T, _Rows, _Options> reqs[], int length, MPI_Status *status)
    // {
    //     int outIndex;
    //     MPI_Request mpi_reqs[length];
    //     for(int i = 0; i < length; i++) {
    //         mpi_reqs[i] = reqs[i].req;
    //     }
    //     comm.waitAny(length, mpi_reqs, &outIndex);
    //     MPI_Dynamic_req req = reqs[outIndex];
    //     if(req.buffer.size() < *req.length) {
    //         req.buffer.conservativeResize(*req.length);
    //     }

    //     return comm.recv(req.buffer, reqs.count, req.source, req.tag, status);
    // }
}
