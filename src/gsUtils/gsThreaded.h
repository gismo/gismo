/** @file gsThreaded.h

    @brief Wrapper for thread-local data members

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris
*/

namespace gismo
{

namespace util
{

template<class C>
class gsThreaded
{
    std::vector<C> m_array;
public:
    gsThreaded() : m_array(omp_get_max_threads()) { }

    // copying?

    // Usage:
    // gsThreaded<C> a;
    // a().function();
    C&       operator ()() { return m_array[omp_get_thread_num()]; }
    const C& operator ()() const { return m_array[omp_get_thread_num()]; }
    //operator C&() const { return m_array[omp_get_thread_num()]; }

    C&       get() { return m_array[omp_get_thread_num()]; }
    const C& get() const { return m_array[omp_get_thread_num()]; }

    C& operator = (C other) { return m_array[omp_get_thread_num()] = give(other); }

};//gsThreaded



}//util

}//gismo
