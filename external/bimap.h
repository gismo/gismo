/* STL-like bidirectional map.
 * http://www.codeproject.com/Articles/3016/An-STL-like-bidirectional-map
 *
 * (C) 2002-2006 Joaquin M Lopez Munoz (joaquin@tid.es). All rights reserved.
 *
 * Permission is granted to use, distribute and modify this code provided that:
 *   - this copyright notice remain unchanged,
 *   - you submit all changes to the copyright holder and properly mark the
 *     changes so they can be told from the original code,
 *   - credits are given to the copyright holder in the documentation of any
 *     software using this code with the following line:
 *       "Portions copyright 2002-2006 Joaquin M Lopez Munoz (joaquin@tid.es)"
 *
 * The author welcomes any suggestions on the code or reportings of actual
 * use of the code. Please send your comments to joaquin@tid.es.
 *
 * The author makes NO WARRANTY or representation, either express or implied,
 * with respect to this code, its quality, accuracy, merchantability, or
 * fitness for a particular purpose.  This software is provided "AS IS", and
 * you, its user, assume the entire risk as to its quality and accuracy.
 *
 * Changes in v1.1:
 *
 *   - bimap::erase(to::iterator,to::iterator) incorrectly returned an
 *     iterator. Documentation was also erroneous about this point.
 *   - Incorrect use of allocator::allocate and allocator::deallocate
 *     was causing much more memory to be used than necessary.
 *   - Improved language conformance with respect to missing typename
 *     and template keywords, faulty friend declarations and broken
 *     implementation of some features of <iterator> in MSVC++.
 *   - allocator::rebind used if compiler supports it.
 *   - Fixed some non-conformances about construction of allocator
 *     and comparison objects in copy contructors.
 *   - The allocator object used to be protected for no good reason:
 *     changed to private as the rest of internal objects.
 *   - Some tweaks to make the thing compile under GNU GCC and Metrowerks
 *     CodeWarrior.
 *   - GCC didn't like a template parameter and a defined type to have
 *     the same name: I don't know if this is actually standard conformant.
 *
 * Changes in v1.2:
 *
 *   - Fixed the code to make it work under MSVC 7.1. Contributed by
 *     Steve Robb.
 *
 * Changes in v1.3:
 *
 *   - Fixed some incorrect (standardwise) friend declarations.
 *   - Sprinkled this-> throughout to cope with two-phase name lookup.
 *
 * Last modified: October 26th, 2006
 */

#ifndef BIMAP_H_9B698EF9_C6E9_4BC4_A7D2_5B4D71155761
#define BIMAP_H_9B698EF9_C6E9_4BC4_A7D2_5B4D71155761
#define VERSION_BIMAP_9B698EF9_C6E9_4BC4_A7D2_5B4D71155761 0x00010003
        
#ifdef _MSC_VER
#pragma warning(disable:4786)
#endif

#include <algorithm>
#include <functional>
#include <iterator>
#include <set>
#include <stddef.h>
#include <stdexcept>
#include <utility>

/* offsetof cannot be used on non-POD types (standard 18.1.5) altough bimap
 * does it in a safe manner. Starting with GCC 3.1, an annoying warning is
 * issued in this situation. Workarounded it thanks to a tip by Andrew Pollard.
 */

#if defined(__GNUC__)&&(__GNUC__>3||(__GNUC__==3&&__GNUC_MINOR__>= 1))
#define BIMAP_OFFSETOF_9B698EF9_C6E9_4BC4_A7D2_5B4D71155761(type,member) \
(__extension__                         \
  (                                    \
    {                                  \
      type* t=0;                       \
      reinterpret_cast<size_t>(        \
        reinterpret_cast<const char*>( \
          &(t->member)                 \
        )                              \
      );                               \
    }                                  \
  )                                    \
)
#else
#define BIMAP_OFFSETOF_9B698EF9_C6E9_4BC4_A7D2_5B4D71155761(type,member) offsetof(type,member)
#endif

/* MSVC++ 6.0 do not support allocator::rebind; in these cases, the only
 * option is use the original allocator_type unrebound, which VC++ 6.0
 * accepts merrily nevertheless.
 */

#if defined(_MSC_VER)&&_MSC_VER==1200 /* MSVC++ 6.0 */
#define BIMAP_REBIND_9B698EF9_C6E9_4BC4_A7D2_5B4D71155761(type1,type2) type1
#else
#define BIMAP_REBIND_9B698EF9_C6E9_4BC4_A7D2_5B4D71155761(type1,type2) \
typename type1::template rebind<type2>::other
#endif

namespace gismo{ // Angelos M.: was codeproject

namespace bimap_detail{

/* Template helper to check for type equality (save possibly constness.)
 * Refs:
 * - Alexandrescu, A.: Modern C++ Design: Generic Programming and Design
 *   Patterns Applied, ch. 2, Addison-Wesley, February 2001.
 * - Boost, type_traits library, boost::is_same<T,U> template,
 *   April 2001, http://boost.org/libs/type_traits/index.htm.
 */

typedef char         equal_types_yes;
struct               equal_types_no {char m[2];};
equal_types_no       equal_types_helper(...);
template<typename T>
equal_types_yes      equal_types_helper(T,T);
template <typename T>
struct equal_types_ptr
{
  static const T* make();
};

template<typename T,typename U> struct equal_types
{
  enum{
    value=
      (sizeof(equal_types_helper(equal_types_ptr<T>::make(),equal_types_ptr<U>::make()))==
        sizeof(equal_types_yes))};
};

/* Template stuff to select one or other type based on a compile-time
 * condition. This can be used to simulate PTS through derivation.
 * Refs:
 * - Marcus, M., Jones, J.: "Simulated partial Specialization for C++",
 *   September 2000, http://opensource.adobe.com/project4/project.shtml
 * - Czarnecki, K., Eisenecker, U.: Generative Programming - Methods, Tools,
 *   and Applications, Addison-Wesley, June 2000.
 */

struct select_then
{
  template<class then,class els> struct result
  {
    typedef then type;
  };
};

struct select_else
{
  template<class then,class els> struct result
  {
    typedef els type;
  };
};

template<bool test> struct selector_switch
{
  typedef select_then result;
};

template<> struct selector_switch<false>
{
  typedef select_else result;
};

template<bool test,typename then,typename els>
struct select
{
  typedef typename selector_switch<test>::result        sel;
  typedef typename sel::template result<then,els>::type result;
};

} /* namespace bimap_detail */

/* inv_pair provides the symmetrical counterpart to std::pair necessary for
 * the to memberspace of bimap. Its layout matches that of std::pair in the
 * sense that an inv_pair<T,U> can be reinterpret_cast'ed to be an std::pair<U,T>.
 * To preserve symmetry, std::pair's should provide the corresponding casting
 * operator to inv_pair. As this is not feasible (predefined classes cannot be
 * injected this type of operators), a derivation of std::pair called
 * direct_pair is defined that plays the role of std::pair.
 */

#if defined(_MSC_VER)&&_MSC_VER>=1300 /* MSVC++ 7.0 */
#pragma warning(push)
#pragma warning(disable:4512)
/* see http://support.microsoft.com/default.aspx?scid=kb;EN-US;Q87638 */
#endif
 
template<typename first_type,typename second_type>
struct direct_pair;

template<typename first_type_,typename second_type_>
struct inv_pair
{
  typedef first_type_  first_type;
  typedef second_type_ second_type;

  second_type second;
  first_type  first;

  inv_pair():second(second_type()),first(first_type()){}
  inv_pair(const first_type& first,const second_type& second):second(second),first(first){}
  inv_pair(const std::pair<first_type,second_type>& r):second(r.second),first(r.first){}
  template<typename F,typename S>
  inv_pair(const inv_pair<F,S>& r):second(r.second),first(r.first){}

  operator direct_pair<second_type,first_type>&()
  {
    return *reinterpret_cast<direct_pair<second_type,first_type> *>(this);
  }

  operator const direct_pair<second_type,first_type>&()const
  {
    return *reinterpret_cast<const direct_pair<second_type,first_type> *>(this);
  }
};

template<typename first_type,typename second_type>
struct direct_pair:public std::pair<first_type,second_type>
{
private:
  typedef std::pair<first_type,second_type> super;

public:
  direct_pair():super(first_type(),second_type()){}
  direct_pair(const first_type& first,const second_type& second):super(first,second){}
  direct_pair(const inv_pair<first_type,second_type>& r):super(r.first,r.second){}
  template<typename F,typename S>
    direct_pair(const std::pair<F,S>& r):super(r.first,r.second){}

  operator inv_pair<second_type,first_type>&()
  {
    return *reinterpret_cast<inv_pair<second_type,first_type> *>(this);
  }

  operator const inv_pair<second_type,first_type>&()const
  {
    return *reinterpret_cast<const inv_pair<second_type,first_type> *>(this);
  }
};

#if defined(_MSC_VER)&&_MSC_VER>=1300 /* MSVC++ 7.0 */
#pragma warning(pop)
#endif

template<typename first_type,typename second_type>
inv_pair<first_type,second_type> make_inv_pair(
    const first_type& first,
    const second_type& second)
{
  return inv_pair<first_type,second_type>(first,second);
}

/* comparison operators for inv_pair */

/* == */

template<typename first_type,typename second_type>
bool operator==(
  const inv_pair<first_type,second_type>& x,
  const inv_pair<first_type,second_type>& y)
{
  return x.first==y.first&&x.second==y.second;
}

template<typename first_type,typename second_type>
bool operator==(
  const inv_pair<first_type,second_type>& x,
  const std::pair<first_type,second_type>& y)
{
  return x.first==y.first&&x.second==y.second;
}

template<typename first_type,typename second_type>
bool operator==(
  const std::pair<first_type,second_type>& x,
  const inv_pair<first_type,second_type>& y)
{
  return x.first==y.first&&x.second==y.second;
}

/* != */

template<typename first_type,typename second_type>
bool operator!=(
  const inv_pair<first_type,second_type>& x,
  const inv_pair<first_type,second_type>& y)
{
  return !(x==y);
}

template<typename first_type,typename second_type>
bool operator!=(
  const inv_pair<first_type,second_type>& x,
  const std::pair<first_type,second_type>& y)
{
  return !(x==y);
}

template<typename first_type,typename second_type>
bool operator!=(
  const std::pair<first_type,second_type>& x,
  const inv_pair<first_type,second_type>& y)
{
  return !(x==y);
}

/* < */

template<typename first_type,typename second_type>
bool operator<(
  const inv_pair<first_type,second_type>& x,
  const inv_pair<first_type,second_type>& y)
{
    return x.first<y.first|| (x.first==y.first&&x.second<y.second);
}

template<typename first_type,typename second_type>
bool operator<(
  const inv_pair<first_type,second_type>& x,
  const std::pair<first_type,second_type>& y)
{
    return x.first<y.first||(x.first==y.first&&x.second<y.second);
}

template<typename first_type,typename second_type>
bool operator<(
  const std::pair<first_type,second_type>& x,
  const inv_pair<first_type,second_type>& y)
{
    return x.first<y.first||(x.first==y.first&&x.second<y.second);
}

/* > */

template<typename first_type,typename second_type>
bool operator>(
  const inv_pair<first_type,second_type>& x,
  const inv_pair<first_type,second_type>& y)
{
  return y<x;
}

template<typename first_type,typename second_type>
bool operator>(
  const inv_pair<first_type,second_type>& x,
  const std::pair<first_type,second_type>& y)
{
  return y<x;
}

template<typename first_type,typename second_type>
bool operator>(
  const std::pair<first_type,second_type>& x,
  const inv_pair<first_type,second_type>& y)
{
  return y<x;
}

/* <= */

template<typename first_type,typename second_type>
bool operator<=(
  const inv_pair<first_type,second_type>& x,
  const inv_pair<first_type,second_type>& y)
{
  return !(y<x);
}

template<typename first_type,typename second_type>
bool operator<=(
  const inv_pair<first_type,second_type>& x,
  const std::pair<first_type,second_type>& y)
{
  return !(y<x);
}

template<typename first_type,typename second_type>
bool operator<=(
  const std::pair<first_type,second_type>& x,
  const inv_pair<first_type,second_type>& y)
{
  return !(y<x);
}

/* >= */

template<typename first_type,typename second_type>
bool operator>=(
  const inv_pair<first_type,second_type>& x,
  const inv_pair<first_type,second_type>& y)
{
  return !(x<y);
}

template<typename first_type,typename second_type>
bool operator>=(
  const inv_pair<first_type,second_type>& x,
  const std::pair<first_type,second_type>& y)
{
  return !(x<y);
}

template<typename first_type,typename second_type>
bool operator>=(
  const std::pair<first_type,second_type>& x,
  const inv_pair<first_type,second_type>& y)
{
  return !(x<y);
}

/* bimap_base is the common base for all bimaps */

struct bimap_base
{
  class value_not_found:public std::logic_error
  {
  public:
    value_not_found():logic_error("value not found"){}
    value_not_found(const std::string& str):logic_error(str){}
  };

  class duplicate_value:public std::logic_error
  {
  public:
    duplicate_value():logic_error("duplicate value"){}
    duplicate_value(const std::string& str):logic_error(str){}
  };

protected:
  /* from_binding and company implement operator[]s. They're
   * defined outside from and to memberspaces because
   * otherwise VC++ chokes and produces an internal compiler
   * error under some circumstances.
   */

  template<typename bimap_type>
  class from_binding
  {
  public:
    from_binding(bimap_type& bm,const typename bimap_type::from_type& f):bm(bm),f(f){}

    const typename bimap_type::to_type& get()const
    {
      typename bimap_type::from_set::const_iterator it=bm.fset.find(&f);
      if(it==bm.fset.end())throw value_not_found();
      return bimap_type::element_by_from(*it)->second;
    }

    operator const typename bimap_type::to_type&()const
    {
      return get();
    }

    from_binding& operator=(const typename bimap_type::to_type& t)
    {
      /* MSVC++ chokes if a ctor for typename bimap_type::element
       * is called directly.
       */

      typedef typename bimap_type::element bimap_element;

      typename bimap_type::fset_iterator fit=bm.fset.find(&f);
      typename bimap_type::tset_iterator tit=bm.tset.find(&t);
      if(tit!=bm.tset.end()){ /* v.second shouldn't be already in */
        /* small chance the pair (f,t) is already inserted */
        if(fit!=bm.fset.end()&&
           bimap_type::element_by_from(*fit)==bimap_type::element_by_to(*tit)){
          return *this;
        }
        else throw duplicate_value();
      }

      bimap_element *                    pne=0;
      typename bimap_type::tset_iterator tnit=bm.tset.end();
      try{
        pne=bm.new_element(bimap_element(f,t));
        tnit=bm.tset.insert(&pne->second).first;
        if(fit==bm.fset.end()){
          bm.fset.insert(&pne->first);
          return *this;
        }
        else{ // rebound fit
          bimap_element * pe=bimap_type::element_by_from(*fit);
          bm.tset.erase(bm.tset.find(&pe->second));
          bm.delete_element(pe);
          const_cast<typename bimap_type::from_type *&>(*fit)=
            &const_cast<typename bimap_type::from_type &>(pne->first);
          return *this;
        }
      }catch(...){
        if(tnit!=bm.tset.end())bm.tset.erase(tnit);
        if(pne)bm.delete_element(pne);
        throw;
      }         
    }

  private:
    from_binding& operator=(const from_binding&);

    bimap_type&                          bm;
    const typename bimap_type::from_type f;
  };

  template<typename bimap_type>
  class const_from_binding
  {
  public:
    const_from_binding(const bimap_type& bm,const typename bimap_type::from_type& f):
      bm(bm),f(f)
    {}

    const typename bimap_type::to_type& get()const
    {
      typename bimap_type::from_set::const_iterator it=bm.fset.find(&f);
      if(it==bm.fset.end())throw value_not_found();
      return bimap_type::element_by_from(*it)->second;
    }

    operator const typename bimap_type::to_type&()const
    {
      return get();
    }

  private:
    const_from_binding& operator=(const const_from_binding&);

    const bimap_type&                    bm;
    const typename bimap_type::from_type f;
  };

  template<typename bimap_type>
  class to_binding
  {
  public:
    to_binding(bimap_type& bm,const typename bimap_type::to_type& t):bm(bm),t(t){}

    const typename bimap_type::from_type& get()const
    {
      typename bimap_type::to_set::const_iterator it=bm.tset.find(&t);
      if(it==bm.tset.end())throw value_not_found();
      return bimap_type::element_by_to(*it)->first;
    }

    operator const typename bimap_type::from_type&()const
    {
      return get();
    }

    to_binding& operator=(const typename bimap_type::from_type& f)
    {
      /* MSVC++ chokes if a ctor for typename bimap_type::element
       * is called directly.
       */

      typedef typename bimap_type::element bimap_element;

      typename bimap_type::tset_iterator tit=bm.tset.find(&t);
      typename bimap_type::fset_iterator fit=bm.fset.find(&f);
      if(fit!=bm.fset.end()){ /* v.second shouldn't be already in */
        /* small chance the pair (f,t) is already inserted */
        if(tit!=bm.tset.end()&&
           bimap_type::element_by_to(*tit)==bimap_type::element_by_from(*fit)){
          return *this;
        }
        else throw duplicate_value();
      }

      bimap_element *                    pne=0;
      typename bimap_type::fset_iterator fnit=bm.fset.end();
      try{
        pne=bm.new_element(bimap_element(f,t));
        fnit=bm.fset.insert(&pne->first).first;
        if(tit==bm.tset.end()){
          bm.tset.insert(&pne->second);
          return *this;
        }
        else{ // rebound tit
          bimap_element * pe=bimap_type::element_by_to(*tit);
          bm.fset.erase(bm.fset.find(&pe->first));
          bm.delete_element(pe);
          const_cast<typename bimap_type::to_type *&>(*tit)=
            &const_cast<typename bimap_type::to_type &>(pne->second);
          return *this;
        }
      }catch(...){
        if(fnit!=bm.fset.end())bm.fset.erase(fnit);
        if(pne)bm.delete_element(pne);
        throw;
      }         
    }

  private:
    to_binding& operator=(const to_binding&);

    bimap_type&                        bm;
    const typename bimap_type::to_type t;
  };

  template<typename bimap_type>
  class const_to_binding
  {
  public:
    const_to_binding(const bimap_type& bm,const typename bimap_type::to_type& t):
      bm(bm),t(t)
    {}

    const typename bimap_type::from_type& get()const
    {
      typename bimap_type::to_set::const_iterator it=bm.tset.find(&t);
      if(it==bm.tset.end())throw value_not_found();
      return bimap_type::element_by_to(*it)->first;
    }

    operator const typename bimap_type::from_type&()const
    {
      return get();
    }

  private:
    const_to_binding& operator=(const const_to_binding&);

    const bimap_type&                  bm;
    const typename bimap_type::to_type t;
  };
};

/* prebimap holds the entire code for bimap except for the 
 * global memberspace, which is constructed via simulated PTS
 * in bimap (derived from prebimap).
 */

template<
  typename from_type_,typename to_type_,
  typename from_compare,typename to_compare,
  typename allocator_type_>
class prebimap:public bimap_base
{
private:
  /* Data structure. The bidirectional map is implemented with two sets
   * fset and tset indexing the corresponding members of elements of type
   * direct_pair<from_type,to_type>.
   */

  typedef from_type_ from_type;
  typedef to_type_   to_type;

  template<typename type,typename compare>
  struct p_compare /* pointer compare based on value compare */
  {
    p_compare(const compare& c=compare()):c(c){}
    bool operator()(const type* p1,const type* p2)const{return c(*p1,*p2);}
    compare get_compare()const{return c;}
  private:
    compare c;
  };

  typedef std::set<
    const from_type*,
    p_compare<
      from_type,
      from_compare>,
    BIMAP_REBIND_9B698EF9_C6E9_4BC4_A7D2_5B4D71155761(allocator_type_,const from_type*)>
                                            from_set;
  typedef std::set<
    const to_type*,
    p_compare<
      to_type,
      to_compare>,
    BIMAP_REBIND_9B698EF9_C6E9_4BC4_A7D2_5B4D71155761(allocator_type_,const to_type*)>
                                            to_set;
  typedef typename from_set::allocator_type fset_allocator_type;
  typedef typename to_set::allocator_type   tset_allocator_type;
  typedef typename from_set::iterator       fset_iterator;
  typedef typename from_set::const_iterator const_fset_iterator;
  typedef typename to_set::iterator         tset_iterator;
  typedef typename to_set::const_iterator   const_tset_iterator;
  typedef direct_pair<
    const from_type,
    const to_type>                          element;
  allocator_type_                           allocator;
  from_set                                  fset;
  to_set                                    tset;

  /* Basic data management */

  /* new_element and delete_element deal only with allocation/
   * deallocation of elements, i.e. they do not update fset and tset.
   */

  element * new_element(const element& e)
  {
    element * pe=allocator.allocate(1,0);
    try{
      allocator.construct(pe,e);
    }catch(...){
      allocator.destroy(pe);
      throw;
    }
    return pe;
  }

  void delete_element(element *pe)
  {
    allocator.destroy(pe);
    allocator.deallocate(pe,1);
  }

  static element * element_by_from(const from_type* pf)
  {
    return
      reinterpret_cast<element*>(
        reinterpret_cast<char*>(
          const_cast<from_type *>(pf))-
            BIMAP_OFFSETOF_9B698EF9_C6E9_4BC4_A7D2_5B4D71155761(element,first));
  }

  static element * element_by_to(const to_type* pt)
  {
    return
      reinterpret_cast<element *>(
        reinterpret_cast<char*>(
          const_cast<to_type *>(pt))-
            BIMAP_OFFSETOF_9B698EF9_C6E9_4BC4_A7D2_5B4D71155761(element,second));
  }

public:
  /* memberspace from */

  class from_impl
  {
  public:

    /* assigning a from_impl object is the same as assigning its owner */

    from_impl& operator=(const from_impl& r)
    {
      prebimap tmp(r.owner());
      swap(tmp.from);
      return *this;
    }

    /* Comparison */

    bool operator==(const from_impl& r)const
    {
      return size()==r.size()&&std::equal(begin(),end(),r.begin());
    }

    bool operator!=(const from_impl& r)const
    {
      return !(*this==r);
    }

    bool operator<(const from_impl& r)const
    {
      return
        std::lexicographical_compare(
          begin(),end(),r.begin(),r.end());
    }

    bool operator>(const from_impl& r)const
    {
      return r<*this;
    }

    bool operator<=(const from_impl& r)const
    {
      return !(r<*this);
    }

    bool operator>=(const from_impl& r)const
    {
      return !(*this<r);
    }

    /* Standard member types */

    typedef from_type_      key_type;
    typedef to_type_        mapped_type;
    typedef to_type_        referent_type; /* prestandard synonim */
    typedef to_type_        data_type;     /* prestandard synonim */
    typedef from_compare    key_compare;
    typedef allocator_type_ allocator_type;
    typedef
      direct_pair<
        const from_type_,
        const to_type_>     value_type;

    /* value_compare lexicographically orders value_type's. This is
     * compatible with the weaker value_compare implemented by maps.
     */

    class value_compare:public std::binary_function<value_type,value_type,bool>
    {
    public:
      bool operator()(const value_type& x,const value_type& y)
      {
        if(kcomp(x.first,y.first))return true;
        if(kcomp(y.first,x.first))return false;
        return tcomp(x.second,y.second);
      }
    protected:
      value_compare(key_compare kcomp,to_compare tcomp):kcomp(kcomp),tcomp(tcomp){}
      key_compare kcomp;
      to_compare  tcomp;
    };

    typedef typename allocator_type_::size_type       size_type;
    typedef typename allocator_type_::difference_type difference_type;
    typedef value_type *                              pointer;
    typedef const value_type *                        const_pointer;
    typedef value_type&                               reference;
    typedef const value_type&                         const_reference;

    /* Iterators */

    class const_iterator;
    class iterator:public std::iterator<std::bidirectional_iterator_tag,const value_type>
    {
      friend class from_impl;
      friend class const_iterator;
      friend class prebimap<from_type,to_type,from_compare,to_compare,allocator_type>;

#if defined(_MSC_VER)&&_MSC_VER==1200 /* MSVC++ 6.0 */
      /* MSVC++ 6.0 fails to define reference in std::iterator */

      typedef const value_type& reference;
#endif

    public:
      iterator(){}

#ifdef __MWERKS__
      /* strange bug */

      typename reference operator*()const{return *element_by_from(*fit);}
      typename value_type * operator->()const{return &operator*();}
#else
      typename iterator::reference operator*()const{return *element_by_from(*fit);}
      typename iterator::value_type * operator->()const{return &operator*();}
#endif

      iterator& operator++(){++fit;return *this;}
      const iterator operator++(int){const iterator tmp=*this;++*this;return tmp;}
      iterator& operator--(){--fit;return *this;}
      const iterator operator--(int){const iterator tmp=*this;--*this;return tmp;}
      bool operator==(const iterator& it)const{return fit==it.fit;}
      bool operator!=(const iterator& it)const{return !(*this==it);}
    private:
      iterator(const fset_iterator& fit):fit(fit){}
      fset_iterator fit;
    };

    class const_iterator:public std::iterator<std::bidirectional_iterator_tag,const value_type>
    {
      friend class from_impl;
      friend class prebimap<from_type,to_type,from_compare,to_compare,allocator_type>;

#if defined(_MSC_VER)&&_MSC_VER==1200 /* MSVC++ 6.0 */
      /* MSVC++ 6.0 fails to define reference in std::iterator */

      typedef const value_type& reference;
#endif

    public:
      const_iterator(){}
      const_iterator(const typename prebimap::from_impl::iterator& it):fit(it.fit){} // Angelos: bug fix for MSVC

#ifdef __MWERKS__
      /* strange bug */

      typename reference operator*()const{return *element_by_from(*fit);}
      const typename value_type * operator->()const{return &operator*();}
#else
      typename const_iterator::reference operator*()const{return *element_by_from(*fit);}
      const typename const_iterator::value_type * operator->()const{return &operator*();}
#endif

      const_iterator& operator++(){++fit;return *this;}
      const const_iterator operator++(int){const const_iterator tmp=*this;++*this;return tmp;}
      const_iterator& operator--(){--fit;return *this;}
      const const_iterator operator--(int){const const_iterator tmp=*this;--*this;return tmp;}
      bool operator==(const const_iterator& it)const{return fit==it.fit;}
      bool operator!=(const const_iterator& it)const{return !(*this==it);}
    private:
      const_iterator(const const_fset_iterator& fit):fit(fit){}
      const_fset_iterator fit;
    };

#if defined(_MSC_VER)&&_MSC_VER==1200 /* MSVC++ 6.0 */
    typedef
      std::reverse_bidirectional_iterator<
        const_iterator,
        const value_type,
        const_reference,
        const value_type *,
        difference_type>                          reverse_iterator;
#else
    typedef std::reverse_iterator<const_iterator> reverse_iterator;
#endif

    typedef reverse_iterator                      const_reverse_iterator;

    /* Iterator retrieval methods */

    iterator begin(){return iterator(owner().fset.begin());}
    const_iterator begin()const{return const_iterator(owner().fset.begin());}
    iterator end(){return iterator(owner().fset.end());}
    const_iterator end()const{return const_iterator(owner().fset.end());}

    reverse_iterator rbegin(){return reverse_iterator(end());}
    const_reverse_iterator rbegin()const{return const_reverse_iterator(end());}
    reverse_iterator rend(){return reverse_iterator(begin());}
    const_reverse_iterator rend()const{return const_reverse_iterator(begin());}
    
    /* Utility standard methods */

    size_type size()const{return owner().fset.size();}
    size_type max_size()const{return owner().fset.max_size();}
    bool empty()const{return owner().fset.empty();}
    allocator_type get_allocator()const{return owner().allocator;}

    /* operator []. Uses wrapper classes from_binding and const_from_binding. */

    from_binding<prebimap<from_type_,to_type_,from_compare,to_compare,allocator_type_> >
    operator[](const from_type_& f)
    {
      return
        from_binding<prebimap<from_type_,to_type_,from_compare,to_compare,allocator_type_> >
          (owner(),f);
    }

    const_from_binding<prebimap<from_type_,to_type_,from_compare,to_compare,allocator_type_> >
    operator[](const from_type_& f)const
    {
      return 
        const_from_binding<prebimap<from_type_,to_type_,from_compare,to_compare,allocator_type_> >
          (owner(),f);
    }

    /* Insertion and erasing */

    std::pair<iterator,bool> insert(const value_type& x)
    {
      fset_iterator fit=owner().fset.find(&x.first);
      if(fit!=owner().fset.end()){
        return std::make_pair(iterator(fit),false);
      }
      tset_iterator tit=owner().tset.find(&x.second);
      if(tit!=owner().tset.end())throw duplicate_value();

      element * pe=0;
      tset_iterator tnit=owner().tset.end();
      fset_iterator fnit=owner().fset.end();
      try{
        pe=owner().new_element(x);
        tnit=owner().tset.insert(&pe->second).first;
        fnit=owner().fset.insert(&pe->first).first;
      }catch(...){
        if(tnit!=owner().tset.end())owner().tset.erase(tnit);
        if(pe)owner().delete_element(pe);
        throw;
      }
      return std::make_pair(iterator(fnit),true);
    }

    iterator insert(iterator it,const value_type& x)
    {
      if(!adjacent(it.fit,x.first))return insert(x).first;

      tset_iterator tit=owner().tset.find(&x.second);
      if(tit!=owner().tset.end())throw duplicate_value();

      element * pe=0;
      tset_iterator tnit=owner().tset.end();
      fset_iterator fnit=owner().fset.end();
      try{
        pe=owner().new_element(x);
        tnit=owner().tset.insert(&pe->second).first;
        fnit=owner().fset.insert(it.fit,&pe->first);
      }catch(...){
        if(tnit!=owner().tset.end())owner().tset.erase(tnit);
        if(pe)owner().delete_element(pe);
        throw;
      }
      return iterator(fnit);
    }

    template<typename it_type>
    void insert(it_type first,it_type last)
    {
      while(first!=last){
        insert(*first);
        ++first;
      }
    }


#ifdef _MSC_VER
    /* The standard says the return type for iterator-based erases
     * in associative containers is void. Strangely enough, VC++
     * implementation returns an iterator. We keep the iterator return
     * for the first version of erase.
     */

    iterator erase(iterator it)
    {
      fset_iterator& fit=it.fit;
      element *      pe=element_by_from(*fit);
      tset_iterator  tit=owner().tset.find(&pe->second);
      owner().delete_element(pe);
      owner().tset.erase(tit);
      return(iterator(owner().fset.erase(fit)));
    }
#else
    void erase(iterator it)
    {
      fset_iterator& fit=it.fit;
      element *      pe=element_by_from(*fit);
      tset_iterator  tit=owner().tset.find(&pe->second);
      owner().delete_element(pe);
      owner().tset.erase(tit);
      owner().fset.erase(fit);
    }
#endif

    void erase(iterator first,iterator last)
    {
      while(first!=last)erase(first++);
    }

    size_type erase(const key_type& key)
    {
      fset_iterator fit=owner().fset.find(&key);
      if(fit==owner().fset.end())return 0;
      element * pe=element_by_from(*fit);
      owner().tset.erase(owner().tset.find(&pe->second));
      owner().fset.erase(fit);
      owner().delete_element(pe);
      return 1;
    }
    
    void clear()
    {
      erase(begin(),end());
    }

    void swap(from_impl& x)
    {
      /* assumes allocator equivalence */

      owner().fset.swap(x.owner().fset);
      owner().tset.swap(x.owner().tset);
    }

    /* Search methods */

    key_compare key_comp()const
    {
      return owner().fset.key_comp().get_compare();
    }

    value_compare value_comp()const
    {
      return 
        value_compare(
          owner().fset.key_comp().get_compare(),
          owner().tset.key_comp().get_compare());
    }

    iterator find(const key_type& key)
    {
      return iterator(owner().fset.find(&key));
    }

    const_iterator find(const key_type& key)const
    {
      return const_iterator(owner().fset.find(&key));
    }

    size_type count(const key_type& key)const
    {
      return owner().fset.count(&key);
    }

    iterator lower_bound(const key_type& key)
    {
      return iterator(owner().fset.lower_bound(&key));
    }

    const_iterator lower_bound(const key_type& key)const
    {
      return const_iterator(owner().fset.lower_bound(&key));
    }

    iterator upper_bound(const key_type& key)
    {
      return iterator(owner().fset.upper_bound(&key));
    }

    const_iterator upper_bound(const key_type& key)const
    {
      return const_iterator(owner().fset.upper_bound(&key));
    }

    std::pair<iterator,iterator> equal_range(const key_type& key)
    {
      return std::make_pair(lower_bound(key),upper_bound(key));
    }

    std::pair<const_iterator,const_iterator> equal_range(const key_type& key)const
    {
      return std::make_pair(lower_bound(key),upper_bound(key));
    }

  protected:
    from_impl(){}
    from_impl(const from_impl&);

    prebimap& owner()
    {
      return *reinterpret_cast<prebimap*>(
        reinterpret_cast<char*>(this)-
          BIMAP_OFFSETOF_9B698EF9_C6E9_4BC4_A7D2_5B4D71155761(prebimap,from));
    };
    const prebimap& owner()const
    {
      return *reinterpret_cast<const prebimap*>(
        reinterpret_cast<const char*>(this)-
          BIMAP_OFFSETOF_9B698EF9_C6E9_4BC4_A7D2_5B4D71155761(prebimap,from));
    };

    bool adjacent(const_fset_iterator fit,const from_type_& f)const
    {
      if(fit==owner().fset.end()){
        if(owner().fset.size()==0)return true;
        const_fset_iterator fit2=fit;
        --fit2;
        return owner().fset.key_comp()(*fit2,&f);
      }
      else if(owner().fset.key_comp()(&f,*fit)){
        if(fit==owner().fset.begin())return true;
        const_fset_iterator fit2=fit;
        --fit2;
        return owner().fset.key_comp()(*fit2,&f);

      }
      else if(owner().fset.key_comp()(*fit,&f)){
        const_fset_iterator fit2=fit;
        ++fit2;
        if(fit2==owner().fset.end())return true;
        return owner().fset.key_comp()(&f,*fit2);
      }
      else return false;
    }
  };
  
  friend class from_impl;

  class from:public from_impl
  {
    friend class prebimap<from_type_,to_type_,from_compare,to_compare,allocator_type_>;

  public:
    friend void swap(from& x,from& y)
    {
      x.swap(y);
    }
  }from; /* from memberspace */

#ifdef __MWERKS__
  /* strange bug */

  friend class from_impl::iterator;
  friend class from_impl::const_iterator;
#elif defined(_MSC_VER) 
  friend typename from::iterator;
  friend typename from::const_iterator;
#endif

  friend class from_binding<prebimap>;
  friend class const_from_binding<prebimap>;

  /* memberspace to, symetrical to from */

  class to_impl
  {
  public:

    to_impl& operator=(const to_impl& r)
    {
      prebimap tmp(r.owner());
      swap(tmp.to);
      return *this;
    }

    /* Comparison */

    bool operator==(const to_impl& r)const
    {
      return size()==r.size()&&std::equal(begin(),end(),r.begin());
    }

    bool operator!=(const to_impl& r)const
    {
      return !(*this==r);
    }

    bool operator<(const to_impl& r)const
    {
      return
        std::lexicographical_compare(
          begin(),end(),r.begin(),r.end());
    }

    bool operator>(const to_impl& r)const
    {
      return r<*this;
    }

    bool operator<=(const to_impl& r)const
    {
      return !(r<*this);
    }

    bool operator>=(const to_impl& r)const
    {
      return !(*this<r);
    }
    
    /* Standard member types */

    typedef to_type_        key_type;
    typedef from_type_      mapped_type;
    typedef from_type_      referent_type; /* prestandard synonim */
    typedef from_type_      data_type;     /* prestandard synonim */
    typedef to_compare      key_compare;
    typedef allocator_type_ allocator_type;
    typedef
      inv_pair<
        const to_type_,
        const from_type_>   value_type;

    class value_compare:public std::binary_function<value_type,value_type,bool>
    {
    public:
      bool operator()(const value_type& x,const value_type& y)
      {
        if(kcomp(x.first,y.first))return true;
        if(kcomp(y.first,x.first))return false;
        return fcomp(x.second,y.second);
      }
    protected:
      value_compare(key_compare kcomp,from_compare fcomp):kcomp(kcomp),fcomp(fcomp){}
      key_compare  kcomp;
      from_compare fcomp;
    };

    typedef typename allocator_type_::size_type       size_type;
    typedef typename allocator_type_::difference_type difference_type;
    typedef value_type *                              pointer;
    typedef const value_type *                        const_pointer;
    typedef value_type&                               reference;
    typedef const value_type&                         const_reference;

    /* Iterators */

    class const_iterator;
    class iterator:public std::iterator<std::bidirectional_iterator_tag,const value_type>
    {
      friend class to_impl;
      friend class const_iterator;
      friend class prebimap<from_type,to_type,from_compare,to_compare,allocator_type>;

#if defined(_MSC_VER)&&_MSC_VER==1200 /* MSVC++ 6.0 */
      /* MSVC++ 6.0 fails to define reference in std::iterator */

      typedef const value_type& reference;
#endif
    public:
      iterator(){}

#ifdef __MWERKS__
      /* strange bug */

      typename reference operator*()const{return *element_by_to(*tit);}
      typename value_type * operator->()const{return &operator*();}
#else
      typename iterator::reference operator*()const{return *element_by_to(*tit);}
      typename iterator::value_type * operator->()const{return &operator*();}
#endif

      iterator& operator++(){++tit;return *this;}
      const iterator operator++(int){const iterator tmp=*this;++*this;return tmp;}
      iterator& operator--(){--tit;return *this;}
      const iterator operator--(int){const iterator tmp=*this;--*this;return tmp;}
      bool operator==(const iterator& it)const{return tit==it.tit;}
      bool operator!=(const iterator& it)const{return !(*this==it);}
    private:
      iterator(const tset_iterator& tit):tit(tit){}
      tset_iterator tit;  
    };

    class const_iterator:public std::iterator<std::bidirectional_iterator_tag,const value_type>
    {
      friend class to_impl;
      friend class prebimap<from_type,to_type,from_compare,to_compare,allocator_type>;

#if defined(_MSC_VER)&&_MSC_VER==1200 /* MSVC++ 6.0 */
      /* MSVC++ 6.0 fails to define reference in std::iterator */

      typedef const value_type& reference;
#endif

    public:
      const_iterator(){}
	const_iterator(const typename prebimap::to_impl::iterator& it):tit(it.tit){} // Angelos, MSVC fix

#ifdef __MWERKS__
      /* strange bug */

      typename reference operator*()const{return *element_by_to(*tit);}
      const typename value_type * operator->()const{return &operator*();}
#else
      typename const_iterator::reference operator*()const{return *element_by_to(*tit);}
      const typename const_iterator::value_type * operator->()const{return &operator*();}
#endif

      const_iterator& operator++(){++tit;return *this;}
      const const_iterator operator++(int){const const_iterator tmp=*this;++*this;return tmp;}
      const_iterator& operator--(){--tit;return *this;}
      const const_iterator operator--(int){const const_iterator tmp=*this;--*this;return tmp;}
      bool operator==(const const_iterator& it)const{return tit==it.tit;}
      bool operator!=(const const_iterator& it)const{return !(*this==it);}
    private:
      const_iterator(const const_tset_iterator& tit):tit(tit){}
      const_tset_iterator tit;
    };

#if defined(_MSC_VER)&&_MSC_VER==1200 /* MSVC++ 6.0 */
    typedef
      std::reverse_bidirectional_iterator<
        const_iterator,
        const value_type,
        const_reference,
        const value_type *,
        difference_type>                          reverse_iterator;
#else
    typedef std::reverse_iterator<const_iterator> reverse_iterator;
#endif

    typedef reverse_iterator                      const_reverse_iterator;

    /* Iterator retrieval methods */

    iterator begin(){return iterator(owner().tset.begin());}
    const_iterator begin()const{return const_iterator(owner().tset.begin());}
    iterator end(){return iterator(owner().tset.end());}
    const_iterator end()const{return const_iterator(owner().tset.end());}

    reverse_iterator rbegin(){return reverse_iterator(end());}
    const_reverse_iterator rbegin()const{return const_reverse_iterator(end());}
    reverse_iterator rend(){return reverse_iterator(begin());}
    const_reverse_iterator rend()const{return const_reverse_iterator(begin());}
    
    /* Utility standard methods */

    size_type size()const{return owner().tset.size();}
    size_type max_size()const{return owner().tset.max_size();}
    bool empty()const{return owner().tset.empty();}
    allocator_type get_allocator()const{return owner().allocator;}

    /* operator []. Uses wrapper classes to_binding and const_to_binding. */

    to_binding<prebimap<from_type_,to_type_,from_compare,to_compare,allocator_type_> >
    operator[](const to_type_& t)
    {
      return
        to_binding<prebimap<from_type_,to_type_,from_compare,to_compare,allocator_type_> >
          (owner(),t);
    }

    const_to_binding<prebimap<from_type_,to_type_,from_compare,to_compare,allocator_type_> >
    operator[](const to_type_& t)const
    {
      return 
        const_to_binding<prebimap<from_type_,to_type_,from_compare,to_compare,allocator_type_> >
          (owner(),t);
    }

    /* Insertion and erasing */

    std::pair<iterator,bool> insert(const value_type& x)
    {
      tset_iterator tit=owner().tset.find(&x.first);
      if(tit!=owner().tset.end()){
        return std::make_pair(iterator(tit),false);
      }
      fset_iterator fit=owner().fset.find(&x.second);
      if(fit!=owner().fset.end())throw duplicate_value();

      element * pe=0;
      fset_iterator fnit=owner().fset.end();
      tset_iterator tnit=owner().tset.end();
      try{
        pe=owner().new_element(x);
        fnit=owner().fset.insert(&pe->first).first;
        tnit=owner().tset.insert(&pe->second).first;
      }catch(...){
        if(fnit!=owner().fset.end())owner().fset.erase(fnit);
        if(pe)owner().delete_element(pe);
        throw;
      }
      return std::make_pair(iterator(tnit),true);
    }

    iterator insert(iterator it,const value_type& x)
    {
      if(!adjacent(it.tit,x.first))return insert(x).first;

      fset_iterator fit=owner().fset.find(&x.second);
      if(fit!=owner().fset.end())throw duplicate_value();

      element * pe=0;
      fset_iterator fnit=owner().fset.end();
      tset_iterator tnit=owner().tset.end();
      try{
        pe=owner().new_element(x);
        fnit=owner().fset.insert(&pe->first).first;
        tnit=owner().tset.insert(it.tit,&pe->second);
      }catch(...){
        if(fnit!=owner().fset.end())owner().fset.erase(fnit);
        if(pe)owner().delete_element(pe);
        throw;
      }
      return iterator(tnit);
    }

    template<typename it_type>
    void insert(it_type first,it_type last)
    {
      while(first!=last){
        insert(*first);
        ++first;
      }
    }


#ifdef _MSC_VER
    /* see note in from::erase */

    iterator erase(iterator it)
    {
      tset_iterator& tit=it.tit;
      element *      pe=element_by_to(*tit);
      fset_iterator  fit=owner().fset.find(&pe->first);
      owner().delete_element(pe);
      owner().fset.erase(fit);
      return(iterator(owner().tset.erase(tit)));
    }
#else
    void erase(iterator it)
    {
      tset_iterator& tit=it.tit;
      element *      pe=element_by_to(*tit);
      fset_iterator  fit=owner().fset.find(&pe->first);
      owner().delete_element(pe);
      owner().fset.erase(fit);
      owner().tset.erase(tit);
    }
#endif

    void erase(iterator first,iterator last)
    {
      while(first!=last)erase(first++);
    }

    size_type erase(const key_type& key)
    {
      tset_iterator tit=owner().tset.find(&key);
      if(tit==owner().tset.end())return 0;
      element * pe=element_by_to(*tit);
      owner().fset.erase(owner().fset.find(&pe->first));
      owner().tset.erase(tit);
      owner().delete_element(pe);
      return 1;
    }
    
    void clear()
    {
      erase(begin(),end());
    }

    void swap(to_impl& x)
    {
      /* assumes allocator equivalence */

      owner().tset.swap(x.owner().tset);
      owner().fset.swap(x.owner().fset);
    }

    /* Search methods */

    key_compare key_comp()const
    {
      return owner().tset.key_comp().get_compare();
    }

    value_compare value_comp()const
    {
      return 
        value_compare(
          owner().tset.key_comp().get_compare(),
          owner().fset.key_comp().get_compare());
    }

    iterator find(const key_type& key)
    {
      return iterator(owner().tset.find(&key));
    }

    const_iterator find(const key_type& key)const
    {
      return const_iterator(owner().tset.find(&key));
    }

    size_type count(const key_type& key)const
    {
      return owner().tset.count(&key);
    }

    iterator lower_bound(const key_type& key)
    {
      return iterator(owner().tset.lower_bound(&key));
    }

    const_iterator lower_bound(const key_type& key)const
    {
      return const_iterator(owner().tset.lower_bound(&key));
    }

    iterator upper_bound(const key_type& key)
    {
      return iterator(owner().tset.upper_bound(&key));
    }

    const_iterator upper_bound(const key_type& key)const
    {
      return const_iterator(owner().tset.upper_bound(&key));
    }

    std::pair<iterator,iterator> equal_range(const key_type& key)
    {
      return std::make_pair(lower_bound(key),upper_bound(key));
    }

    std::pair<const_iterator,const_iterator> equal_range(const key_type& key)const
    {
      return std::make_pair(lower_bound(key),upper_bound(key));
    }

  protected:
    to_impl(){}
    to_impl(const to_impl&);

    prebimap& owner()
    {
      return *reinterpret_cast<prebimap*>(
        reinterpret_cast<char*>(this)-
          BIMAP_OFFSETOF_9B698EF9_C6E9_4BC4_A7D2_5B4D71155761(prebimap,to));
    };
    const prebimap& owner()const
    {
      return *reinterpret_cast<const prebimap*>(
        reinterpret_cast<const char*>(this)-
          BIMAP_OFFSETOF_9B698EF9_C6E9_4BC4_A7D2_5B4D71155761(prebimap,to));
    };

    bool adjacent(const_tset_iterator tit,const to_type_& t)const
    {
      if(tit==owner().tset.end()){
        if(owner().tset.size()==0)return true;
        const_tset_iterator tit2=tit;
        --tit2;
        return owner().tset.key_comp()(*tit2,&t);
      }
      else if(owner().tset.key_comp()(&t,*tit)){
        if(tit==owner().tset.begin())return true;
        const_tset_iterator tit2=tit;
        --tit2;
        return owner().tset.key_comp()(*tit2,&t);

      }
      else if(owner().tset.key_comp()(*tit,&t)){
        const_tset_iterator tit2=tit;
        ++tit2;
        if(tit2==owner().tset.end())return true;
        return owner().tset.key_comp()(&t,*tit2);
      }
      else return false;
    }
  };
  
  friend class to_impl;

  class to:public to_impl
  {
    friend class prebimap<from_type_,to_type_,from_compare,to_compare,allocator_type_>;

  public:
    friend void swap(to& x,to& y)
    {
      x.swap(y);
    }
  }to; /* to memberspace */

#ifdef __MWERKS__
  /* strange bug */

  friend class to_impl::iterator;
  friend class to_impl::const_iterator;
#elif defined(_MSC_VER) 
  friend typename to::iterator;
  friend typename to::const_iterator;
#endif

  friend class to_binding<prebimap>;
  friend class const_to_binding<prebimap>;

  /* Double-hint insertion. This does not naturally belong into any
   * memberspace.
   */

  std::pair<typename from_impl::iterator,typename to_impl::iterator>
  insert(
    typename from_impl::iterator fit,typename to_impl::iterator tit,
    const typename from_impl::value_type& x)
  {
    typedef typename from_impl::iterator from_iterator;
    typedef typename to_impl::iterator   to_iterator;

    if(!from.adjacent(fit.fit,x.first)){
      fit=fset.find(&x.first);
      if(fit!=fset.end()){ /* small chance x is already inserted */
        tset_iterator tnit=tset.find(&x.second);
        if(tnit!=tset.end()&&element_by_from(*(fit.fit))==element_by_to(*tnit)){
          return std::make_pair(fit,to_iterator(tnit));
        }
        else throw duplicate_value();
      }
    }

    if(!to.adjacent(tit.tit,x.second)){
      tit=tset.find(&x.second);
      if(tit!=tset.end()) throw duplicate_value();
      /* no need to check for x being inserted (already done above) */
    }

    element * pe=0;
    tset_iterator tnit=tset.begin();
    fset_iterator fnit=fset.begin();
    try{
      pe=new_element(x);
      tnit=tset.insert(tit.tit,&pe->second);
      fnit=fset.insert(fit.fit,&pe->first);
    }catch(...){
      if(tnit!=tset.end())tset.erase(tnit);
      if(pe)delete_element(pe);
      throw;
    }
    return std::make_pair(from_iterator(fnit),to_iterator(tnit));
  }

protected:
  prebimap(
      const from_compare& from_comp,
      const to_compare& to_comp,
      const allocator_type_& al):
    allocator(al),
    fset(p_compare<from_type,from_compare>(from_comp),fset_allocator_type(allocator)),
    tset(p_compare<to_type,to_compare>(to_comp),tset_allocator_type(allocator))
  {}

  prebimap(const prebimap& r):
    allocator(r.allocator),
    fset(r.fset.key_comp(),fset_allocator_type(allocator)),
    tset(r.tset.key_comp(),tset_allocator_type(allocator))
  {
    try{

#if defined(_MSC_VER)&&_MSC_VER==1200 /* MSVC++ 6.0 */
/* Amortized constant insertion in VC++ 6.0 happens if hint iterator is
 * right *after* the insertion point, in disagreement with the standard.
 */

      typename from::iterator fhint=from.end();
      typename to::iterator   thint=to.end();
      for(typename from::const_iterator it=r.from.begin();it!=r.from.end();++it){
        insert(fhint,thint,*it);
      }
#else
      typename from_impl::iterator fhint=from.end();
      typename to_impl::iterator   thint=to.end();
      for(typename from_impl::const_iterator it=r.from.begin();it!=r.from.end();++it){
        std::pair<typename from_impl::iterator,typename to_impl::iterator>
          p=insert(fhint,thint,*it);
        fhint=p.first;
        thint=p.second;
      }
#endif

    }catch(...){
      from.clear();
      throw;
    }
  }

  ~prebimap()
  {
    for(fset_iterator it=fset.begin();it!=fset.end();++it){
      delete_element(element_by_from(*it));
    }
  }
};

#ifdef __GNUC__
/* GCC chokes when trying to derive from a template class with memberspaces.
 * See http://gcc.gnu.org/cgi-bin/gnatsweb.pl?cmd=view%20audit-trail&database=gcc&pr=9159
 * for details. prebimap_identity is nothing but compiler sugar to workaround
 * the problem.
 */

template<
  typename from_type,typename to_type,
  typename from_compare,typename to_compare,
  typename allocator_type>
struct prebimap_identity
{
  typedef prebimap<from_type,to_type,from_compare,to_compare,allocator_type> type; 
};
#endif

/* When from_type and to_type are equal, we promote only the from memberspace
 * to the global memberspace ans members of the to space posing no ambiguity.
 */
template<
  typename from_type_,typename to_type_,
  typename from_compare,typename to_compare,
  typename allocator_type_>
class bimap_equal_types:

#ifdef __GNUC__
  public prebimap_identity<from_type_,to_type_,from_compare,to_compare,allocator_type_>::type
#else
  public prebimap<from_type_,to_type_,from_compare,to_compare,allocator_type_>
#endif

{
  typedef 
    prebimap<from_type_,to_type_,from_compare,to_compare,allocator_type_> super;
  typedef typename super::template from_binding<super>                    from_binding_t;
  typedef typename super::template const_from_binding<super>              const_from_binding_t;
  typedef typename super::to::value_type                                  to_value_type;
  typedef typename super::to::iterator                                    to_iterator;

public:
  typedef typename  super::from::key_type               key_type;
  typedef typename  super::from::mapped_type            mapped_type;
  typedef typename  super::from::referent_type          referent_type;
  typedef typename  super::from::data_type              data_type;
  typedef typename  super::from::key_compare            key_compare;
  typedef typename  super::from::allocator_type         allocator_type;
  typedef typename  super::from::value_type             value_type;
  typedef typename  super::from::value_compare          value_compare;
  typedef typename  super::from::size_type              size_type;
  typedef typename  super::from::difference_type        difference_type;
  typedef typename  super::from::pointer                pointer;
  typedef typename  super::from::const_pointer          const_pointer;
  typedef typename  super::from::reference              reference;
  typedef typename  super::from::const_reference        const_reference;
  typedef typename  super::from::iterator               iterator;
  typedef typename  super::from::const_iterator         const_iterator;
  typedef typename  super::from::reverse_iterator       reverse_iterator;
  typedef typename  super::from::const_reverse_iterator const_reverse_iterator;

  iterator                 begin(){return this->from.begin();}
  const_iterator           begin()const{return this->from.begin();}
  iterator                 end(){return this->from.end();}
  const_iterator           end()const{return this->from.end();}
  reverse_iterator         rbegin(){return this->from.rbegin();}
  const_reverse_iterator   rbegin()const{return this->from.rbegin();}
  reverse_iterator         rend(){return this->from.rend();}
  const_reverse_iterator   rend()const{return this->from.rend();}
  size_type                size()const{return this->from.size();}
  size_type                max_size()const{return this->from.max_size();}
  bool                     empty()const{return this->from.empty();}
  allocator_type           get_allocator()const{return this->from.get_allocator();}
  from_binding_t           operator[](const key_type& key){return this->from[key];}
  const_from_binding_t     operator[](const key_type& key)const{return this->from[key];}
  using                    super::insert;
  std::pair<iterator,bool> insert(const value_type& x){return this->from.insert(x);}
  iterator                 insert(iterator it,const value_type& x){return this->from.insert(it,x);}
  template<
    typename it_type>
  void                     insert(it_type first,it_type last){this->from.insert(first,last);}

#ifdef _MSC_VER /* see note in from::erase */
  iterator                 erase(iterator it){return this->from.erase(it);}
#else
  void                     erase(iterator it){this->from.erase(it);}
#endif

  void                     erase(iterator first,iterator last){this->from.erase(first,last);}
  size_type                erase(const key_type& key){return this->from.erase(key);}
  void                     clear(){this->from.clear();}
  void                     swap(bimap_equal_types& x){this->from.swap(x.from);}
  friend void              swap(bimap_equal_types& x,bimap_equal_types& y){x.swap(y);}
  key_compare              key_comp()const{return this->from.key_comp();}
  value_compare            value_comp()const{return this->from.value_comp();}
  iterator                 find(const key_type& key){return this->from.find(key);}
  const_iterator           find(const key_type& key)const{return this->from.find(key);}
  size_type                count(const key_type& key)const{return this->from.count(key);}
  iterator                 lower_bound(const key_type& key){return this->from.lower_bound(key);}
  const_iterator           lower_bound(const key_type& key)const{return this->from.lower_bound(key);}
  iterator                 upper_bound(const key_type& key){return this->from.upper_bound(key);}
  const_iterator           upper_bound(const key_type& key)const{return this->from.upper_bound(key);}
  std::pair<
    iterator,
    iterator>              equal_range(const key_type& key){return this->from.equal_range(key);}
  std::pair<
    const_iterator,
    const_iterator>        equal_range(const key_type& key)const{return this->from.equal_range(key);}

  /* Promotion of unambiguous parts of the to memberspace */

  std::pair<
    to_iterator,
    bool>                  insert(const to_value_type& x){return this->to.insert(x);}
  to_iterator              insert(to_iterator it,const to_value_type& x){return this->to.insert(it,x);}
  void                     insert(const to_value_type *first,const to_value_type *last)
                           {this->to.insert(first,last);}

#ifdef _MSC_VER /* see note in from::erase */
  to_iterator              erase(to_iterator it){return this->to.erase(it);}
#else
  void                     erase(to_iterator it){this->to.erase(it);}
#endif

  void                     erase(to_iterator first,to_iterator last){this->to.erase(first,last);}

protected:
  bimap_equal_types(
      const from_compare& from_comp,
      const to_compare& to_comp,
      const allocator_type_& al):
    super(from_comp,to_comp,al)
  {
  }

  /* default copy ctor serves well */
};

/* If from_type and to_type are distinct, some more members of the to
 * memberspace can be promoted to the global memberspace without
 * ambiguity (notably operator[]).
 */
template<
  typename from_type_,typename to_type_,
  typename from_compare,typename to_compare,
  typename allocator_type_>
class bimap_different_types:

#ifdef __GNUC__
  public prebimap_identity<from_type_,to_type_,from_compare,to_compare,allocator_type_>::type
#else
  public prebimap<from_type_,to_type_,from_compare,to_compare,allocator_type_>
#endif

{
  typedef 
    prebimap<from_type_,to_type_,from_compare,to_compare,allocator_type_> super;
  typedef typename super::template from_binding<super>                    from_binding_t;
  typedef typename super::template const_from_binding<super>              const_from_binding_t;
  typedef typename super::template to_binding<super>                      to_binding_t;
  typedef typename super::template const_to_binding<super>                const_to_binding_t;
  typedef typename super::to::key_type                                    to_key_type;
  typedef typename super::to::value_type                                  to_value_type;
  typedef typename super::to::iterator                                    to_iterator;
  typedef typename super::to::const_iterator                              const_to_iterator;

public:
  typedef typename super::from::key_type               key_type;
  typedef typename super::from::mapped_type            mapped_type;
  typedef typename super::from::referent_type          referent_type;
  typedef typename super::from::data_type              data_type;
  typedef typename super::from::key_compare            key_compare;
  typedef typename super::from::allocator_type         allocator_type;
  typedef typename super::from::value_type             value_type;
  typedef typename super::from::value_compare          value_compare;
  typedef typename super::from::size_type              size_type;
  typedef typename super::from::difference_type        difference_type;
  typedef typename super::from::pointer                pointer;
  typedef typename super::from::const_pointer          const_pointer;
  typedef typename super::from::reference              reference;
  typedef typename super::from::const_reference        const_reference;
  typedef typename super::from::iterator               iterator;
  typedef typename super::from::const_iterator         const_iterator;
  typedef typename super::from::reverse_iterator       reverse_iterator;
  typedef typename super::from::const_reverse_iterator const_reverse_iterator;

  iterator                 begin(){return this->from.begin();}
  const_iterator           begin()const{return this->from.begin();}
  iterator                 end(){return this->from.end();}
  const_iterator           end()const{return this->from.end();}
  reverse_iterator         rbegin(){return this->from.rbegin();}
  const_reverse_iterator   rbegin()const{return this->from.rbegin();}
  reverse_iterator         rend(){return this->from.rend();}
  const_reverse_iterator   rend()const{return this->from.rend();}
  size_type                size()const{return this->from.size();}
  size_type                max_size()const{return this->from.max_size();}
  bool                     empty()const{return this->from.empty();}
  allocator_type           get_allocator()const{return this->from.get_allocator();}
  from_binding_t           operator[](const key_type& key){return this->from[key];}
  const_from_binding_t     operator[](const key_type& key)const{return this->from[key];}
  using                    super::insert;
  std::pair<iterator,bool> insert(const value_type& x){return this->from.insert(x);}
  iterator                 insert(iterator it,const value_type& x){return this->from.insert(it,x);}
  template<
    typename it_type>
  void                     insert(it_type first,it_type last){this->from.insert(first,last);}

#ifdef _MSC_VER /* see note in from::erase */
  iterator                 erase(iterator it){return this->from.erase(it);}
#else
  void                     erase(iterator it){this->from.erase(it);}
#endif

  void                     erase(iterator first,iterator last){this->from.erase(first,last);}
  size_type                erase(const key_type& key){return this->from.erase(key);}
  void                     clear(){this->from.clear();}
  void                     swap(bimap_different_types& x){this->from.swap(x.from);}
  friend void              swap(bimap_different_types& x,bimap_different_types& y)
                           {x.swap(y);}
  key_compare              key_comp()const{return this->from.key_comp();}
  value_compare            value_comp()const{return this->from.value_comp();}
  iterator                 find(const key_type& key){return this->from.find(key);}
  const_iterator           find(const key_type& key)const{return this->from.find(key);}
  size_type                count(const key_type& key)const{return this->from.count(key);}
  iterator                 lower_bound(const key_type& key){return this->from.lower_bound(key);}
  const_iterator           lower_bound(const key_type& key)const{return this->from.lower_bound(key);}
  iterator                 upper_bound(const key_type& key){return this->from.upper_bound(key);}
  const_iterator           upper_bound(const key_type& key)const{return this->from.upper_bound(key);}
  std::pair<
    iterator,
    iterator>              equal_range(const key_type& key){return this->from.equal_range(key);}
  std::pair<
    const_iterator,
    const_iterator>        equal_range(const key_type& key)const{return this->from.equal_range(key);}

  /* Promotion of unambiguous parts of the to memberspace */

  to_binding_t             operator[](const to_key_type& key){return this->to[key];}
  const_to_binding_t       operator[](const to_key_type& key)const{return this->to[key];}
  std::pair<
    to_iterator,
    bool>                  insert(const to_value_type& x){return this->to.insert(x);}
  to_iterator              insert(to_iterator it,const to_value_type& x){return this->to.insert(it,x);}
  void                     insert(const to_value_type *first,const to_value_type *last)
                           {this->to.insert(first,last);}

#ifdef _MSC_VER /* see note in this->from.erase */
  to_iterator              erase(to_iterator it){return this->to.erase(it);}
#else
  void                     erase(to_iterator it){this->to.erase(it);}
#endif

  void                     erase(to_iterator first,to_iterator last){this->to.erase(first,last);}
  size_type                erase(const to_key_type& key){return this->to.erase(key);}
  to_iterator              find(const to_key_type& key){return this->to.find(key);}
  const_to_iterator        find(const to_key_type& key)const{return this->to.find(key);}
  size_type                count(const to_key_type& key)const{return this->to.count(key);}
  to_iterator              lower_bound(const to_key_type& key){return this->to.lower_bound(key);}
  const_to_iterator        lower_bound(const to_key_type& key)const{return this->to.lower_bound(key);}
  to_iterator              upper_bound(const to_key_type& key){return this->to.upper_bound(key);}
  const_to_iterator        upper_bound(const to_key_type& key)const{return this->to.upper_bound(key);}
  std::pair<
    to_iterator,
    to_iterator>           equal_range(const to_key_type& key){return this->to.equal_range(key);}
  std::pair<
    const_to_iterator,
    const_to_iterator>     equal_range(const to_key_type& key)const{return this->to.equal_range(key);}

protected:
  bimap_different_types(
      const from_compare& from_comp,
      const to_compare& to_comp,
      const allocator_type_& al):
    super(from_comp,to_comp,al)
  {
  }

  /* default copy ctor serves well */
};

/* bimap finally inherits from bimap_equal_types or bimap_different_types
 * depending on the equality of from_type and to_type, as a way to
 * simulate PTS.
 */
template<
  typename from_type_,typename to_type_,
  typename from_compare=std::less<from_type_>,
  typename to_compare=std::less<to_type_>,
  typename allocator_type=std::allocator<direct_pair<const from_type_,const to_type_> > >
class bimap:
  public 
    bimap_detail::select<
      bimap_detail::equal_types<from_type_,to_type_>::value,
      bimap_equal_types<
        from_type_,to_type_,
        from_compare,to_compare,allocator_type>,
      bimap_different_types<
        from_type_,to_type_,
        from_compare,to_compare,allocator_type>
    >::result
{
protected:
  typedef typename
    bimap_detail::select<
      bimap_detail::equal_types<from_type_,to_type_>::value,
      bimap_equal_types<
        from_type_,to_type_,
        from_compare,to_compare,allocator_type>,
      bimap_different_types<
        from_type_,to_type_,
        from_compare,to_compare,allocator_type>
    >::result
    super;

public:
  explicit bimap(
    const from_compare& from_comp=from_compare(),
    const to_compare& to_comp=to_compare(),
    const allocator_type& al=allocator_type()):
    super(from_comp,to_comp,al)
  {
  }

  /* default copy ctor serves well */

  bimap& operator=(const bimap& r)
  {
    bimap tmp(r);
    this->swap(tmp);
    return *this;
  }

  /* inverse copy ctor (from a bimap<to_type,from_type>) */

#if defined(_MSC_VER)&&_MSC_VER==1200 /* MSVC++ 6.0 */
  /* no allocator::rebind, assume allocator_type==std::allocator */

  typedef 
    bimap<
      to_type_,from_type_,
      to_compare,from_compare,
      std::allocator<direct_pair<const to_type_,const from_type_> > >
    inv_bimap;

  explicit bimap(const inv_bimap& r):
    super(r.to.key_comp(),r.from.key_comp(),allocator_type())
#else
  typedef 
    bimap<
      to_type_,from_type_,
      to_compare,from_compare,
      typename allocator_type::template rebind<
        direct_pair<const to_type_,const from_type_> >::other>
    inv_bimap;

  explicit bimap(const inv_bimap& r):
    super(r.to.key_comp(),r.from.key_comp(),r.get_allocator())
#endif

/* body of bimap(const inv_bimap& r) follows */

  {
    try{

#if defined(_MSC_VER)&&_MSC_VER==1200 /* MSVC++ 6.0 */
/* Amortized constant insertion in VC++ 6.0 happens if hint iterator is
 * right *after* the insertion point, in disagreement with the standard.
 */
      typename super::from::iterator fhint=from.end();
      typename super::to::iterator   thint=to.end();
      for(typename inv_bimap::to::const_iterator it=r.to.begin();it!=r.to.end();++it){
        insert(fhint,thint,*it);
      }
#else
      typename super::from::iterator fhint=this->from.end();
      typename super::to::iterator   thint=this->to.end();
      for(typename inv_bimap::to::const_iterator it=r.to.begin();it!=r.to.end();++it){
        std::pair<typename super::from::iterator,typename super::to::iterator>
          p=this->insert(fhint,thint,*it);
        fhint=p.first;
        thint=p.second;
      }
#endif

    }catch(...){
      this->clear();
      throw;
    }
  }

  template<typename it_type>
  bimap(
    it_type first,it_type last,
    const from_compare& from_comp=from_compare(),
    const to_compare& to_comp=to_compare(),
    const allocator_type& al=allocator_type()):
    super(from_comp,to_comp,al)
  {
    try{

#if defined(_MSC_VER)&&_MSC_VER==1200 /* MSVC++ 6.0 */
/* Amortized constant insertion in VC++ 6.0 happens if hint iterator is
 * right *after* the insertion point, in disagreement with the standard.
 */

      typename super::from::iterator fhint=from.end();
      typename super::to::iterator   thint=to.end();
      while(first!=last){
        insert(fhint,thint,first++);
      }
#else
      typename super::from::iterator fhint=this->from.end();
      typename super::to::iterator   thint=this->to.end();
      while(first!=last){
        std::pair<typename super::from::iterator,typename super::to::iterator>
          p=this->insert(fhint,thint,first++);
        fhint=p.first;
        thint=p.second;
      }
#endif

    }catch(...){
      this->clear();
      throw;
    }
  }

  /* Comparison: we simply forward to from */

  bool operator==(const bimap& r)const{return this->from==r.from;}
  bool operator!=(const bimap& r)const{return this->from!=r.from;}
  bool operator< (const bimap& r)const{return this->from<r.from;}
  bool operator> (const bimap& r)const{return this->from>r.from;}
  bool operator<=(const bimap& r)const{return this->from<=r.from;}
  bool operator>=(const bimap& r)const{return this->from>=r.from;}
};

} /* namespace gismo */

#elif VERSION_BIMAP_9B698EF9_C6E9_4BC4_A7D2_5B4D71155761!=0x00010003
#error You have included two BIMAP.H with different version numbers
#endif
