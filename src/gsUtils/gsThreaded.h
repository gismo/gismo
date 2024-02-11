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
template<class C, class Allocator = std::allocator<C> >
class gsThreaded
{
#ifdef _OPENMP
    std::vector< std::vector<C,Allocator> > m_array;
    #else
    C m_c;
#endif

public:

#ifdef _OPENMP
    gsThreaded()
    : m_array(omp_get_max_active_levels(), std::vector<C>(omp_get_max_threads()))
    { }

    /// Casting to the local data
    inline operator C&()             { return mine(); }
    inline operator const C&() const { return mine(); }

    /// Returning the local data
    inline const C& mine() const
    {
        //omp_get_ancestor_thread_num(omp_get_active_level())
        return m_array[omp_get_active_level()][omp_get_thread_num()];
    }

    inline C&       mine()
    {return const_cast<C &>(static_cast<const gsThreaded &>(*this).mine());}

    /// Assigning to the local data
    C& operator = (C other) { return mine() = give(other); }
#else
    /// Casting to the local data
    inline operator C&()             { return m_c; }
    inline operator const C&() const { return m_c; }

    /// Returning the local data
    inline C&       mine() { return m_c; }
    inline const C& mine() const { return m_c; }

    /// Assigning to the local data
    C& operator = (C other) { return m_c = give(other); }
#endif
    
};//gsThreaded

}//util

}//gismo
