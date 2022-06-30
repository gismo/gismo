/** @file gsThreaded.h

    @brief Wrapper for thread-local data members

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris
*/

#pragma once

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
template<class C>
class gsThreaded
{
#ifdef _OPENMP
    std::vector<C> m_array;
    #else
    C m_c;
#endif

public:

#ifdef _OPENMP
    gsThreaded() : m_array(omp_get_max_threads()) { }

    /// Casting to the local data
    operator C&()             { return m_array[omp_get_thread_num()]; }
    operator const C&() const { return m_array[omp_get_thread_num()]; }

    /// Returning the local data
    C&       mine() { return m_array[omp_get_thread_num()]; }
    const C& mine() const { return m_array[omp_get_thread_num()]; }

    /// Assigning to the local data
    C& operator = (C other) { return m_array[omp_get_thread_num()] = give(other); }
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
