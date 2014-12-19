/** @file gsMemory.h

    @brief Provides utility function related to memory management.

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): C. Hofreither, A. Mantzaflaris
*/

# pragma once

#include <memory>

#ifdef TR1_SHARED_PTR_USE_TR1_MEMORY
#include <tr1/memory>
#else
#ifdef BOOST_SHARED_PTR_FOUND
#include <boost/shared_ptr.hpp>
#endif
#endif

namespace gismo {


namespace memory {

// Use the correct shared_ptr
#ifdef STD_SHARED_PTR_FOUND
    using std::shared_ptr;
#else 
#ifdef TR1_SHARED_PTR_FOUND
    using std::tr1::shared_ptr;
#else 
#ifdef BOOST_SHARED_PTR_FOUND
    using boost::shared_ptr;
#else 
    using NOT_FOUND::shared_ptr;
#endif
#endif
#endif

    using std::auto_ptr;

}

/// Takes a T* and wraps it in an auto_ptr. Useful for one-off
/// function return values to avoid memory leaks.
template <typename T>
inline std::auto_ptr<T> safe(T *x)
{
    return std::auto_ptr<T>(x);
}

/** Wrapper class for reference types which are allowed to be "moved
 *  from", i.e., which the caller does not need anymore. This is very
 *  similar to C++11 rvalue references, but without the language
 *  support.
 */
template <typename T>
class gsMovable
{
public:
    T& ref() const
    { return m_value; }

private:
    template <typename U>
    friend gsMovable<U> give(U& x);

    explicit gsMovable(T & x)
        : m_value(x)
    { }

    // disable default constructor
    gsMovable();

    // disable assignment operator
    gsMovable& operator= (const gsMovable& other);

    T & m_value;
};

/** Wrap a T& in a gsMovable<T> to indicate to the callee that the
 *  value is not needed anymore.
 */
template <typename T> 
inline gsMovable<T> give(T & x)
{ return gsMovable<T>(x); }


// Small, dynamically sized arrays on the stack.
// Only use this if the size is guaranteed not to be more than a few
// hundred bytes!
#if defined(__GNUC__)
#define STACK_ARRAY( T, name, sz )    T name[sz];
#else
#define STACK_ARRAY( T, name, sz )    T * name = (T*) alloca ( (sz) * sizeof(T) );
#endif


template <typename It, typename ItOut>
void cloneAll(It start, It end, ItOut out)
{
    for (It i = start; i != end; ++i)
        *out++ = (*i)->clone();
}

template <typename It>
void freeAll(It begin, It end)
{
    for (It it = begin; it != end; ++it)
        if ( *it )
            delete (*it);
}

template <typename Cont>
void freeAll(Cont& cont)
{
    freeAll( cont.begin(), cont.end() );
    cont.clear();
}

}; // namespace gismo
