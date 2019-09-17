/**
 * @file gsMpiDynamicTraits.h
 *
 * @brief Sender and receiver functions for dynamic sized types.
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

    /**
     * @brief Dynamic blocking send for gsVector of any size.
     */
    template<typename T, int _Rows, int _Options>
    static int dsend(gsMpiComm &comm, gsVector<T, _Rows, _Options>* in, int len, int dest, int tag = 0)
    {
        int length = in->size();
        comm.send(&length, 1, dest, tag);
        return comm.send(in->begin(), length, dest, tag);
    }

    /**
     * @brief Dynamic non-blocking send for gsVector of any size.
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
     * @brief Specialized dynamic blocking receive for gsVector of any size.
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

    /**
     * @brief Specialized dynamic non-blocking receive for gsVector of any size.
     *
     * Returns a gsActionableMpiRequest that should be waited on with a dynamic wait function (dwait*** family).
     */
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

    /**
     * @brief Dynamic waiting for any request in the list. Applies/completes the actionable request and returns.
     *
     * Returns the mpi status of the action, and stores the index of the request that has completed in outIndex.
     */
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
}
