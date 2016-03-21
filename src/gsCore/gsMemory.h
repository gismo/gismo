/** @file gsMemory.h

    @brief Provides utility function related to memory management.

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): C. Hofreither, A. Mantzaflaris
*/

#pragma once

#include <memory>

#ifdef TR1_SHARED_PTR_USE_TR1_MEMORY
#include <tr1/memory>
#elif defined(BOOST_SHARED_PTR_FOUND)
#include <boost/shared_ptr.hpp>
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
#elif defined(TR1_SHARED_PTR_FOUND)
    using std::tr1::shared_ptr;
#elif defined(BOOST_SHARED_PTR_FOUND)
    using boost::shared_ptr;
#else 
    using NOT_FOUND::shared_ptr;
#endif

using std::auto_ptr;

// typename unique<C>::ptr
template <typename C>
struct unique
{
#ifdef GISMO_BUILD_CPP11
    typedef std::unique_ptr<C> ptr;
#else
    typedef std::auto_ptr<C>   ptr;
#endif
};

template <typename C>
typename unique<C>::ptr make_unique(C * x)
{ return typename unique<C>::ptr(x); }

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

/**
   Wrapper for a reference that can be swapped with another object.
   Used by the give(.) function to implement argument passing
*/
template <typename T>
class gsMovable
{
public:

    // Moves resources to \a x
    inline void moveTo(T & x) { m_ref.swap(x); m_ref.clear();}

private:
    template<typename U>
    friend gsMovable<U> give(U & x);

    // Only give(.) can create the wrapper
    explicit gsMovable(T & x) : m_ref(x) { }
    
    // disable default constructor
    gsMovable();

    // disable assignment operator
    gsMovable& operator= (const gsMovable& other);

private:
    T & m_ref;
};

/** 
    Helper function for reference types which are allowed to be "moved
    from", i.e., which the caller does not need anymore. This is very
    similar to the use of std::move (C++11) for passing arguments by
    rvalue reference.
 */
template<typename T> inline
gsMovable<T> give(T & x) { return gsMovable<T>(x); }


// Small, dynamically sized arrays on the stack.
// Only use this if the size is guaranteed not to be more than a few
// hundred bytes! Be warned: overflow occurs without any warning
//#if defined(__GNUC__)
//
//#define STACK_ARRAY( T, name, sz )    T name[sz]; // Note: buggy on some compilers/versions
//#else
#define STACK_ARRAY( T, name, sz )    T * name = (T*) alloca ( (sz) * sizeof(T) );
//#endif


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
    for (typename Cont::iterator it = cont.begin(); 
         it != cont.end(); ++it)
        delete (*it);
    cont.clear();
}

/// \brief Constructs a vector of pointers from a vector of objects
template<typename obj> inline 
std::vector<obj*> asVectorPtr(const std::vector<obj> & matv)
{
    std::vector<obj*> result;
    const size_t d = matv.size();
    result.reserve(d);
    for ( size_t i = 0; i!=d; ++i)
        result.push_back( const_cast<obj*>(&matv[i]) );
    return result;
}

/// \brief Casts a vector of pointers 
template <typename Base, typename Derived>
std::vector<Base*> castVectorPtr(std::vector<Derived*> pVec)
{
    std::vector<Base*> result(pVec.size());
    std::copy(pVec.begin(), pVec.end(), result.begin() );
    return result;
}

/// \brief Returns true if all instances of \a Base cast to \a Derived
template <typename Derived, typename Base>
bool checkVectorPtrCast(std::vector<Base*> pVec)
{
    for (typename std::vector<Base*>::iterator it = pVec.begin();
         it != pVec.end(); ++it)
        if ( ! dynamic_cast<Derived*>(*it) )
            return false;
    return true;
}

/**
   \brief Small wrapper for std::copy mimicking memcpy (or
   std::copy_n), copies \a n positions starting from \a begin into
   \a result. The latter is expected to have been allocated in advance
*/
template <class T, class U>
inline void copy_n(const T * begin, const size_t n, U * result)
{
    std::copy(begin, begin+n,
#   ifdef _MSC_VER
              // Take care of C4996 warning
              //stdext::checked_array_iterator<U*>(result,n));
              stdext::unchecked_array_iterator<U*>(result));
#   else
    result);
// Note: in C++11 there is:
// std::copy_n(begin, n, result);
#   endif
}

/// \brief Deleter function that does not delete an object pointer
template <typename Obj>
void null_deleter(Obj *) {}

} // namespace gismo
