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
    static int dsend(gsMpiComm &comm, gsVector<T, _Rows, _Options>& in, int dest, int tag = 0)
    {
        int length = in.size();
        comm.send(&length, 1, dest, tag);
        return comm.send(in.begin(), length, dest, tag);
    }

    /**
     * @brief Dynamic blocking send for gsSparseMatrix of any size.
     */
    template<typename T, int _Options, typename _Index>
    static int dsend(gsMpiComm &comm, gsSparseMatrix<T, _Options, _Index>& in, int dest, int tag = 0)
    {
        int dimensions[2];
        dimensions[0] = in.rows();
        dimensions[1] = in.cols();

        using Data = std::pair<std::pair<_Index, _Index>, T>;

        int size = in.nonZeros();

        Data* data = (Data*) malloc(sizeof(Data) * size);
        int i = 0;
        for(int outer = 0; outer < in.cols(); outer++) {
            auto it = in.begin(outer);
            while(it) {
                std::pair<_Index, _Index> coeff = std::make_pair(it.row(), it.col());
                data[i] = std::make_pair(coeff, it.value());
                ++i;
                ++it;
            }
        }

        comm.send(&dimensions, 2, dest, tag);
        int result = comm.send(data, size, dest, tag);
        free(data);
        return result;
    }

    /**
     * @brief Dynamic non-blocking send for gsVector of any size.
     *
     * reqs should be an array that can hold at least two MPI_Requests.
     */
    template<typename T, int _Rows, int _Options>
    static int disend(gsMpiComm &comm, gsVector<T, _Rows, _Options>& in, int dest, MPI_Request* reqs, int tag = 0)
    {
        int length = in.size();
        comm.isend(&length, 1, dest, &reqs[0], tag);

        return comm.isend(in.begin(), length, dest, &reqs[1], tag);
    }


   /**
     * @brief Dynamic non-blocking send for gsSparseMatrix of any size.
     *
     * reqs should be an array that can hold at least two MPI_Requests.
     */
    template<typename T, int _Options, typename _Index>
    static int disend(gsMpiComm &comm, gsSparseMatrix<T, _Options, _Index>& in, int dest,  MPI_Request* reqs, int tag = 0)
    {
        int dimensions[2];
        dimensions[0] = in.rows();
        dimensions[1] = in.cols();

        using Data = std::pair<std::pair<_Index, _Index>, T>;

        int size = in.nonZeros();

        Data* data = (Data*) malloc(sizeof(Data) * size);
        int i = 0;
        for(int outer = 0; outer < in.cols(); outer++) {
            auto it = in.begin(outer);
            while(it) {
                std::pair<_Index, _Index> coeff = std::make_pair(it.row(), it.col());
                data[i] = std::make_pair(coeff, it.value());
                ++i;
                ++it;
            }
        }

        comm.isend(&dimensions, 2, dest, &reqs[0], tag);
        int result = comm.isend(data, size, dest, &reqs[1], tag);
        free(data);
        return result;
    }

    /**
     * @brief Specialized dynamic blocking receive for gsVector of any size.
     *
     * statuses should be an array that can hold at least two MPI_Requests.
     */
    template<typename T, int _Rows, int _Options>
    static int drecv(gsMpiComm &comm, gsVector<T, _Rows, _Options>& out, int source, int tag = 0, MPI_Status* statuses = NULL)
    {
        int length;
        comm.recv(&length, 1, source, tag, (statuses == NULL ? NULL : &statuses[0]));

        if(out.size() < length) {
            out.resize(length);
        }

        return comm.recv(out.begin(), length, source, tag, (statuses == NULL ? NULL : &statuses[1]));
    }

    /**
     * @brief Specialized dynamic blocking receive for gsSparseMatrix of any size.
     *
     * statuses should be an array that can hold at least two MPI_Requests.
     */
    template<typename T, int _Options, typename _Index>
    static int drecv(gsMpiComm &comm, gsSparseMatrix<T, _Options, _Index>& out, int source, int tag = 0, MPI_Status* statuses = NULL)
    {
        int dimensions[2];
        comm.recv(&dimensions, 2, source, tag, (statuses == NULL ? NULL : &statuses[0]));

        if(dimensions[0] != out.rows() || dimensions[1] != out.cols()) {
            out.resize(dimensions[0], dimensions[1]);
        }

        using Data = std::pair<std::pair<_Index, _Index>, T>;

        gsMpiStatus status = comm.probe(source, tag);
        int size = status.size<Data>();

        Data* data = (Data*) malloc(sizeof(Data) * size);
        int result = comm.recv(data, size, source, tag);

        for(int i = 0; i < size; i++) {
            auto indices = std::get<0>(data[i]);
            T datum = std::get<1>(data[i]);
            out(std::get<0>(indices), std::get<1>(indices)) = datum;
        }

        free(data);

        return result;
    }

    /**
     * @brief Specialized dynamic non-blocking receive for gsVector of any size.
     *
     * Returns a gsActionableMpiRequest that should be waited on with a dynamic wait function (dwait*** family).
     */
    template<typename T, int _Rows, int _Options>
    static gsActionableMpiRequest direcv(gsMpiComm &comm, gsVector<T, _Rows, _Options>& out, int source, int tag = 0)
    {
        int* length = (int*) malloc(sizeof(int));
        gsActionableMpiRequest req;
        comm.irecv(length, 1, source, &req.req, tag);
        req.action = [&comm, &out, length, source, tag](MPI_Status* status) {
            if(out.size() != *length) {
                out.resize(*length);
            }

            int result = comm.recv(out.begin(), *length, source, tag, status);
            free(length);
            return result;
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
