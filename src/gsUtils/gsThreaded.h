/** @file gsThreaded.h

    @brief Wrapper for thread-local data members

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris
*/

#pragma once

#include <vector>
#include <algorithm>

#include <gsCore/gsDebug.h>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace gismo
{

namespace util
{

// Usage:
// gsThreaded<C> a;
// a.mine();
template<class C, class Allocator = std::allocator<C> >
class gsThreaded
{
#ifdef _OPENMP
    std::vector<C,Allocator> m_array;
    #else
    C m_c;
#endif

public:

#ifdef _OPENMP
    // Sized to max(omp_get_max_threads(), omp_get_num_procs()): neither
    // quantity alone is a hard bound. omp_get_num_procs() is the number of
    // processors available to the program, not a cap on the team size a
    // later omp_set_num_threads() may request; omp_get_max_threads() only
    // reflects the setting in force right now, at construction time, and
    // can itself be raised afterwards (a legal, runtime-visible OpenMP API
    // call this class cannot intercept). Taking the max of the two is a
    // generous estimate, not a guarantee, so _slot() still bounds-checks
    // every access; that check is what stands between a count raised past
    // the estimate and a worker thread indexing past the end of m_array.
    gsThreaded()
    : m_array(std::max(omp_get_max_threads(), omp_get_num_procs())) { }

    /// Casting to the local data
    operator C&()             { return m_array[_slot()]; }
    operator const C&() const { return m_array[_slot()]; }

    /// Returning the local data
    C&       mine()           { return m_array[_slot()]; }
    const C& mine() const     { return m_array[_slot()]; }

    /// Assigning to the local data
    C& operator = (C other)   { return m_array[_slot()] = give(other); }

private:
    /// Thread slot of the caller, bounds-checked. OpenMP has no hard upper
    /// bound on the team size (omp_set_num_threads accepts any positive n), so
    /// the constructor's sizing is a generous estimate, not a guarantee; this
    /// check is what stands between a count raised past it and an out-of-bounds
    /// write into the heap. GISMO_ENSURE, not GISMO_ASSERT: it must stay live
    /// under NDEBUG.
    inline size_t _slot() const
    {
        const int t = omp_get_thread_num();
        GISMO_ENSURE(static_cast<size_t>(t) < m_array.size(),
            "gsThreaded: thread id "<<t<<" exceeds the array size "<<m_array.size()
            <<" (omp_get_num_procs()="<<omp_get_num_procs()
            <<", omp_get_max_threads()="<<omp_get_max_threads()<<")");
        return static_cast<size_t>(t);
    }

#else
    /// Casting to the local data
    operator C&()             { return m_c; }
    operator const C&() const { return m_c; }

    /// Returning the local data
    C&       mine() { return m_c; }
    const C& mine() const { return m_c; }

    /// Assigning to the local data
    C& operator = (C other) { return m_c = give(other); }
#endif

};//gsThreaded

}//util

}//gismo
