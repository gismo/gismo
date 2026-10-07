/** @file gsThreaded.h

    @brief Wrapper for thread-local data members

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris
*/

#pragma once

#include <gsCore/gsDebug.h>

#ifdef _OPENMP
#include <omp.h>
#include <algorithm>
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
    /// One slot per thread, also for regions that run after the thread count
    /// was raised (up to the number of cores) after construction.
    gsThreaded() : m_array(std::max(omp_get_max_threads(), omp_get_num_procs())) { }

    /// Casting to the local data
    operator C&()             { return mine(); }
    operator const C&() const { return mine(); }

    /// Returning the local data
    C&       mine()       { return m_array[threadSlot()]; }
    const C& mine() const { return m_array[threadSlot()]; }

    /// Assigning to the local data
    C& operator = (C other) { return mine() = give(other); }

    /// Number of thread slots
    size_t size() const { return m_array.size(); }

private:
    size_t threadSlot() const
    {
        const size_t t = static_cast<size_t>(omp_get_thread_num());
        GISMO_ASSERT(t < m_array.size(), "gsThreaded: thread id " << t
                     << " exceeds the " << m_array.size() << " slots made at construction");
        return t;
    }

public:
#else
    /// Casting to the local data
    operator C&()             { return m_c; }
    operator const C&() const { return m_c; }
    
    /// Returning the local data
    C&       mine() { return m_c; }
    const C& mine() const { return m_c; }
    
    /// Assigning to the local data
    C& operator = (C other) { return m_c = give(other); }

    /// Number of thread slots
    size_t size() const { return 1; }
#endif
    
};//gsThreaded

}//util

}//gismo
