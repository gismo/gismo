/** @file gsMemory.h

    @brief Provides utility function related to memory management.

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): C. Hofreither, A. Mantzaflaris, J. Vogl
*/

#pragma once

#include <gsCore/gsTemplateTools.h>

#ifdef __MINGW32__
//#include <malloc/malloc.h> //xcode
#include <malloc.h>
#endif

#if  __cplusplus < 201103 && defined( __GLIBCXX__ )
#  if defined(__INTEL_COMPILER)
#    include <boost/shared_ptr.hpp>
#    include <boost/weak_ptr.hpp>
#  else
#    include <tr1/memory>
#  endif
#else // libc++ or other
#  include <memory>
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

#if   __cplusplus < 201103 && defined( __GLIBCXX__ )
#  if defined(__INTEL_COMPILER)
using boost::shared_ptr;
using boost::weak_ptr;
#  else
using std::tr1::shared_ptr;
using std::tr1::weak_ptr; 
#  endif
#else // libc++ or other
using std::shared_ptr;
using std::weak_ptr;
#endif 

/* \brief Adaptor for a unique pointer

usage:
\code
memory::unique_ptr<int> B;
\endcode
*/
#if __cplusplus >= 201103 || _MSC_VER >= 1600
using std::unique_ptr;
using std::nullptr_t;
#else

template <typename T>
class unique_ptr : public std::auto_ptr<T>
{
    typedef std::auto_ptr<T> Base;
    typedef std::auto_ptr_ref<T> unique_ptr_ref;

    //struct Cannot_Convert_Pointer;
    
public :
    explicit unique_ptr(T* p = 0)  throw() : Base(p) { }
    
    unique_ptr(const unique_ptr& r) : Base( const_cast<unique_ptr&>(r) ) { }
    
    unique_ptr(unique_ptr_ref m)  throw() : Base(m) { }
    
    template<typename U>
    unique_ptr(const unique_ptr<U> & r
               // unique_ptr<typename conditional<is_base_of<U,T>::value, U,
               //                      Cannot_Convert_Pointer >::type>
        )  throw()
    : Base( const_cast<unique_ptr<U>&>(r) ) { }
    
    unique_ptr & operator=(const unique_ptr& other)  throw()
    {
        Base::operator=(const_cast<unique_ptr&>(other));
        return *this;
    }
    
    template<class U>
    unique_ptr & operator=(const unique_ptr<U> & other) throw()
    {
        Base::operator=(const_cast<unique_ptr<U>&>(other));
        return *this;
    }
    
    //operator shared_ptr<T> () { return shared_ptr<T>(this->release()); }
    
    template<class U> operator shared_ptr<U>()
    // shared_ptr<typename conditional<is_base_of<U,T>::value, U,
    //                      Cannot_Convert_Pointer >::type> ()
    { return shared_ptr<U>(Base::release()); }
    
    bool operator!() const { return Base::get() == NULL; }

private:

    struct SafeBool
    { SafeBool(int) {}
       void dummy() {} };

    typedef void (SafeBool::*bool_cast_type)();

public:

    operator bool_cast_type() const
    { return !Base::get() ? 0 : &SafeBool::dummy; }
};

template<class T>
bool operator==(const unique_ptr<T> & p1, const unique_ptr<T> & p2)
{ return p1.get()==p2.get(); }
template<class T>
bool operator!=(const unique_ptr<T> & p1, const unique_ptr<T> & p2)
{ return p1.get()!=p2.get(); }
template<class T>
bool operator<(const unique_ptr<T> & p1, const unique_ptr<T> & p2)
{ return p1.get()<p2.get(); }
template<class T>
bool operator>(const unique_ptr<T> & p1, const unique_ptr<T> & p2)
{ return p1.get()>p2.get(); }
template<class T>
bool operator<=(const unique_ptr<T> & p1, const unique_ptr<T> & p2)
{ return p1.get()<=p2.get(); }
template<class T>
bool  operator>=(const unique_ptr<T> & p1, const unique_ptr<T> & p2)
{ return p1.get()>=p2.get(); }

class nullptr_t
{
public:
    /* Return 0 for any class pointer */
    template<typename T> operator T*() const {return 0;}
    /* Return 0 for any member pointer */
    template<typename T, typename U> operator T U::*() const {return 0;}
    /* Safe boolean conversion */
    operator void*() const  {return 0;}
private:
    /* Not allowed to get the address */
    void operator&() const;
};
#endif

/// \brief Deleter function that does not delete an object pointer
template <typename T> void null_deleter(T *) {}


/// Takes a T* and wraps it in a shared_ptr. Useful for avoiding
/// memory leaks.
///
/// This has a move semantics: the shared_ptr object takes
/// the ownership.
template <typename T>
inline shared_ptr<T> make_shared(T *x) { return shared_ptr<T>(x); }

/// \brief Creates a shared pointer which does not eventually delete
/// the underlying raw pointer. Usefull to refer to objects which
/// should not be destroyed.
///
/// The caller keeps the ownership.
template <typename T>
inline shared_ptr<T> make_shared_not_owned(const T *x)
{ return shared_ptr<T>(const_cast<T*>(x), null_deleter<T>); }

/// Takes a T* and wraps it in an unique_ptr. Useful for one-off
/// function return values to avoid memory leaks.
///
/// This has a move semantics: the unique_ptr object takes
/// the ownership.
template <typename T>
inline unique_ptr<T> make_unique(T * x) { return unique_ptr<T>(x); }

/// \brief Converts an uPtr \a p to an uPtr
/// of class \a toC and gives it back as return value.
template<class toC, typename from>
inline unique_ptr<toC> convert_ptr(from p)
{ return unique_ptr<toC>( dynamic_cast<toC*>(p.release()) ); }

/// Takes a vector of smart pointers and returns the corresponding raw pointers.
template <typename T>
inline std::vector<T*> get_raw(const std::vector< unique_ptr<T> >& cont)
{
    std::vector<T*> result;
    for (typename std::vector< unique_ptr<T> >::const_iterator it = cont.begin(); it != cont.end(); ++it)
        result.push_back(const_cast<T*>( (*it).get() ));
    return result;
}

/// Takes a vector of smart pointers and returns the corresponding raw pointers.
template <typename T>
inline std::vector<T*> get_raw(const std::vector< shared_ptr<T> >& cont)
{
    std::vector<T*> result;
    for (typename std::vector< shared_ptr<T> >::const_iterator it = cont.begin(); it != cont.end(); ++it)
        result.push_back(const_cast<T*>( (*it).get() ));
    return result;
}

/// Takes a vector of smart pointers, releases them and returns the corresponding raw pointers.
template <typename T>
inline std::vector<T*> release(std::vector< unique_ptr<T> >& cont)
{
    std::vector<T*> result;
    for (typename std::vector< unique_ptr<T> >::iterator it = cont.begin(); it != cont.end(); ++it)
        result.push_back( (*it).release() );
    cont.clear();
    return result;
}

} // namespace memory

#if __cplusplus >= 201103 || _MSC_VER >= 1900
// fix MSVC 2013- (_MSC_VER < 1900)
// MSVC < 1900 do not work probably. give makes a deep copy for return value,
// losses left value. But the alternative code results in segmentation vaults
// because a swap/give loop leads to a stack overflow.
// From the adresses, it seams that Eigen do not support rvalue with MSVC < 1900
// Therefore disabled EIGEN_HAS_RVALUE_REFERENCES for MSVC < 1900 and use
// alternative code.

/** 
    Alias for std::move, to be used instead of writing std::move for
    keeping backward c++98 compatibility
*/

template <class T> inline
auto give(T&& t) -> decltype(std::move(std::forward<T>(t)))
{
#if defined(GISMO_EXTRA_DEBUG) && ! defined(_MSC_VER)
    // TODO: is there way that also MS can check this?
    static_assert( util::has_move_constructor<typename std::remove_reference<T>::type>::value, "There is no move constructor. Copy would be created." );
#endif
    return std::move(std::forward<T>(t));
}

#else
/**
    Alias for std::move, to be used instead of std::move for backward
    c++98 compatibility and MSVC before 2015

    Based on swapping and copy elision.
*/
template <typename S> inline S give(S & x)
{ S t; t.swap(x); return t; }

template <typename T> inline
memory::unique_ptr<T> give(memory::unique_ptr<T> & x)
{ return memory::unique_ptr<T>(x.release()); }

template <typename T> inline
memory::shared_ptr<T> give(memory::shared_ptr<T> & x)
{ memory::shared_ptr<T> result = x; x.reset(); return result; }

#endif

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
/// raw pointers in iterator \a out.
template <typename It, typename ItOut>
void cloneAll(It start, It end, ItOut out)
{
    for (It i = start; i != end; ++i)
        *out++ = dynamic_cast<typename std::iterator_traits<ItOut>::value_type>((*i)->clone().release());
}

/// \brief Clones all pointers in the container \a in and stores them as raw
/// pointers in container \a out
template <typename ContIn, typename ContOut>
void cloneAll(const ContIn& in, ContOut& out)
{
    out.resize(in.size());
    cloneAll(in.begin(), in.end(), out.begin());
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
    for (typename Cont::iterator it = cont.begin(); it != cont.end(); ++it)
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
    for (typename std::vector<Base*>::iterator it = pVec.begin(); it != pVec.end(); ++it)
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

#if __cplusplus < 201103L && _MSC_VER < 1600 && !defined(nullptr)
// Define nullptr for compatibility with newer C++
static const gismo::memory::nullptr_t nullptr ={};
#endif
