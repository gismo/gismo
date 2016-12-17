/** @file gsMemory.h

    @brief Provides utility function related to memory management.

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): C. Hofreither, A. Mantzaflaris
*/

#pragma once

#ifdef __MINGW32__
//#include <malloc/malloc.h> //xcode
#include <malloc.h>
#endif

#if defined(_LIBCPP_VERSION) || (_MSC_VER >= 1700) // libstdc++ ?
#include <memory>
#else // stdc++ (__GLIBCXX__)
//Assume GCC or a GCC-compliant compiler with TR1
#include <tr1/memory>
#endif

namespace gismo {

/** @namespace gismo::memory

    @brief
    This namespace contains functions related to memory management.
    
    \ingroup Core
*/
namespace memory 
{

/* \brief Adaptor for a shared pointer

usage:
\code
memory::shared_ptr<int> B;
\endcode
*/

#if defined(_LIBCPP_VERSION) || (_MSC_VER >= 1700)
using std::shared_ptr;
using std::weak_ptr;
#else
using std::tr1::shared_ptr;
using std::tr1::weak_ptr; 
#endif 

/* \brief Adaptor for a unique pointer

usage:
\code
memory::unique_ptr<int> B;
\endcode
*/
#if __cplusplus >= 201103
using std::unique_ptr;
#else
template <typename T>
class unique_ptr : public std::auto_ptr<T>
{
    typedef std::auto_ptr<T> Base;
    typedef std::auto_ptr_ref<T> unique_ptr_ref;
public :
    typedef T X;

    explicit unique_ptr(X* p = 0):Base(p) { }

    unique_ptr(unique_ptr& r ):Base(r) { }

    unique_ptr(unique_ptr_ref m):Base(m) { }

    using Base::operator=;

};
#endif

/// \brief Deleter function that does not delete an object pointer
template <typename T>
void null_deleter(T *) {}


/// Takes a T* and wraps it in a shared_ptr. Useful for avoiding
/// memory leaks.
template <typename T>
inline shared_ptr<T> make_shared(T *x)
{ return shared_ptr<T>(x); }

/// \brief Creates a shared pointer which does not eventually delete
/// the underlying raw pointer. Usefull to refer to objects which
/// should not be destroyed
template <typename T>
inline shared_ptr<T> make_shared_not_owned(const T *x)
{
    return shared_ptr<T>(const_cast<T*>(x), null_deleter<T>);
}

/// Takes a T* and wraps it in an unique_ptr. Useful for one-off
/// function return values to avoid memory leaks.
template <typename T>
inline unique_ptr<T> make_unique(T * x)
{ return unique_ptr<T>(x); }


} // namespace memory

/// Takes a T* and wraps it in an unique_ptr. Useful for one-off
/// function return values to avoid memory leaks.
/// Does the same as memory::make_unique.
template <typename T>
inline memory::unique_ptr<T> safe(T *x)
{ return memory::unique_ptr<T>(x); }

/**
   Wrapper for a reference that can be swapped with another object.
   Used by the give(.) function to implement argument passing
*/
template <typename T>
class gsMovable
{
public:

    /// Moves resource to \a x
    inline void moveTo(T & x) { m_ref.swap(x); m_ref.clear();}

    /// Read-only access resource
    const T & get() {return m_ref;}

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
template <typename T> inline
gsMovable<T> give(T & x) { return gsMovable<T>(x); }


// Small, dynamically sized arrays on the stack, for POD types.
// Only use this if the size is guaranteed not to be more than a few
// hundred bytes! Be warned: overflow occurs without any warning
#if defined(GISMO_WITH_MPQ) || defined(GISMO_WITH_MPFR)
 #define STACK_ARRAY( T, name, sz )    T name[sz];
#else
// Note: VLAs(following line) can be buggy on some compilers/versions,
// also not nececarily on the stack
// #define STACK_ARRAY( T, name, sz )    T name[sz];
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
   std::copy_n) for a raw pointer destination, copies \a n positions
   starting from \a begin into \a result. The latter is expected to
   have been allocated in advance
*/
template <class T, class U>
inline void copy_n(T begin, const size_t n, U* result)
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

namespace util
{
/**
   \brief Small wrapper for std::copy mimicking std::copy for a raw
   pointer destination, copies \a n positions starting from \a begin
   into \a result. The latter is expected to have been allocated in
   advance
*/
template <class T, class U>
inline void copy(T begin, T end, U* result)
{
    std::copy(begin, end,
#   ifdef _MSC_VER
              // Take care of C4996 warning
              //stdext::checked_array_iterator<U*>(result,n));
              stdext::unchecked_array_iterator<U*>(result));
#   else
    result);
#   endif
}

}

} // namespace gismo
