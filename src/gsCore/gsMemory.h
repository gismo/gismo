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

/** @namespace gismo::memory

    @brief
    This namespace contains functions related to memory management.
    
    \ingroup Core
*/
namespace memory 
{

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
inline memory::auto_ptr<T> safe(T *x)
{
    return memory::auto_ptr<T>(x);
}

/// Takes a T* and wraps it in a shared_ptr. Useful for avoiding
/// memory leaks.
template <typename T>
inline memory::shared_ptr<T> shared(T *x)
{
    return memory::shared_ptr<T>(x);
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

/** Wrap a T& in a gsMovable<T> to indicate to the call that the
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


/// \brief Clones all pointers in the range [\a start \a end) and stores new
/// pointers in iterator \a out
template <typename It, typename ItOut>
void cloneAll(It start, It end, ItOut out)
{
    for (It i = start; i != end; ++i)
        *out++ = (*i)->clone();
}

/// \brief Frees all pointers in the range [\a begin \a end)
template <typename It>
void freeAll(It begin, It end)
{
    for (It it = begin; it != end; ++it)
    {
        delete (*it);
        *it = NULL;
    }
}

/// \brief Frees all pointers in the container \a Cont
template <typename Cont>
void freeAll(Cont& cont)
{
    freeAll( cont.begin(), cont.end() );
    cont.clear();
}

/// \brief Casts a vector of pointers 
template <typename Base, typename Derived>
std::vector<Base*> castVectorPtr(std::vector<Derived*> pVec)
{
    std::vector<Base*> result(pVec.size());
    std::copy(pVec.begin(), pVec.end(), result.begin() );
    return result;
}

/// \brief Small wrapper for std::copy, copies \a n positions starting
/// from \a begin into \a result. The latter is expected to have been
/// allocated in advance
template <class T>
inline void copyRange(const T * begin, T * result, const int n)
{
#   ifdef _MSC_VER
    // Take care of C4996 warning
    std::copy(begin, begin+n,
              stdext::checked_array_iterator<T*>(result,n));
#   else
    std::copy(begin, begin+n, result);
#   endif
}

/// \brief Deleter function that does not delete an object pointer
template <typename Obj>
void null_deleter(Obj *) {}

}; // namespace gismo
