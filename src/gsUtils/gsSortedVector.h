/** @file gsSortedVector.h

    @brief An std::vector with sorting capabilities

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): A. Mantzaflaris, based on 
    https://www.thedigitalmachine.com/wiki/C%2B%2B_sorted_vector
*/

#pragma once

//#include <vector>
//#include <algorithm>			// For lower_bound


namespace std
{

template<typename T1, typename T2>
std::ostream& operator << ( std::ostream& os, 
    			const std::pair<T1,T2>& rhs )
{
    os << rhs.first << ", " << rhs.second;
    return os;
}

}//namespace std

namespace gismo {

/** \brief This class is derived from std::vector, and adds sort tracking.

    There are two basic ways to use a sorted vector:
    
    METHOD 1
    Always maintain sort order by inserting with push_sorted() -
    the location of new items is determined before inserting;
    since the vector remains sorted, this doesn't take too long
    (although for large batch insertions METHOD 2 is definitely
    faster);
    
    METHOD 2
    Allow batch insertion without sorting with push_unsorted(); then
    provide an additional call to sort the vector;  before
    searching for an item, the vector is always sorted if needed;

    Of course you need to provide an operator()< for the type of object
    you're sorting, if it doesn't have one.  Example:
    \code{.cpp}
    class MyClass
    {
    public:
    bool operator< (const MyClass left) const
    {
    if ( left.m_nMostImportant == m_nMostImportant )
    return left.m_nLeastImportant < m_nLeastImportant;
    
    return left.m_nMostImportant < m_nMostImportant;
    }
    
    int m_nMostImportant;
    int m_nLeastImportant;
    }
    \endcode
    
    NOTE: C++ doesn't let you use an operator()< for POINTERS.  This
    breaks down when creating the template code, as you end up with
    a ref to a ref which is not allowed.
    So if you have a vector of pointers, here's what you have to do:
    
    FOR NOW, with straight C++, create a less-than functor, then
    pass that in to the functor versions of the class methods below.
    Create a functor, aka function object, as follows:
    \code{.cpp}
    struct my_class_lessthan
    {
    bool operator()(const MyClass* left, const MyClass* right)
    {
    return left->get_timestamp() < right->get_timestamp();
    }
    };
    \endcode
    Usage example:
    \code{.cpp}
    gsSortedVector<MyClass*> svpMC;
    svpMC.push_unsorted(new MyClass(blah, blah);
    svpMC.push_unsorted(new MyClass(blah, blah);
    vpMC.sort( my_class_lessthan() );
    \endcode
    Once C++0x is available, I need to update this class to use
    a function object wrapper, and allow the user to set it
    in the constructor, then always use it automatically.
    http:en.wikipedia.org/wiki/C%2B%2B0x#Polymorphic_wrappers_for_function_objects
    
    \warning if you change the key value of any object in the vector,
    you have unsorted the array without marking it as such.  Make sure
    you call SetSorted(false) where appropriate.
    
    \ingroup Utils
*/

template<class T, class _A = std::allocator<T> >
class gsSortedVector : public std::vector<T, _A>
{
    typedef std::vector<T,_A> inherited;

public:
    typedef typename std::vector<T>::const_iterator const_iterator;
    typedef typename std::vector<T>::iterator       iterator;

public:

//
    gsSortedVector( )
        : std::vector<T, _A>( ), m_bSorted(true) {  }

    gsSortedVector( typename std::vector<T>::const_iterator start
                   , typename std::vector<T>::const_iterator end)
        : std::vector<T, _A>(start,end), m_bSorted(0) { };


//
    //-------------------------------------------------------------------//
    // SetSorted()																			//
    //-------------------------------------------------------------------//
    // I didn't want to override every constructor to set this
    // member variable, so this function is publicly accessible.
    // You should call SetSorted( true ) right after construction.
    // TO DO: if you feel like it... derive all constructors to avoid
    // the need for this.  There are 4 last time I checked.
    //-------------------------------------------------------------------//
    void SetSorted( bool bSorted = true ) { m_bSorted = bSorted; }

    //-------------------------------------------------------------------//
    // sort()                                                            //
    //-------------------------------------------------------------------//
    // This function sorts the data as needed.  Call it after repeated calls to
    // push_unsorted(), or just let other members call it for you on next access.
    // It calls std::sort(), which defaults to using operator<() for
    // comparisons.
    //-------------------------------------------------------------------//
    void sort()
    {
        if ( !m_bSorted )
        {
            std::sort( inherited::begin(), inherited::end() );
            SetSorted();
        }
    }

    bool bContains( const T& t ) const
    {
      if ( !m_bSorted )
          gsWarn<<"gsSortedVector is not sorted, bContains("<<t<<")"
                << "is not guaranteed to be correct.\n";

        return std::binary_search( inherited::begin(), inherited::end(), t );
    }

    typename std::vector<T>::iterator lower_bound_it( const T& key )
    {
        if ( !m_bSorted )
          sort();

        typename std::vector<T>::iterator it = std::lower_bound( inherited::begin(), inherited::end(), key );
        return it;
    }

    /*const*/ T* lower_bound_ptr( const T& key )
    {
        typename std::vector<T>::iterator it = lower_bound_it( key );

        if (it==inherited::end())
            return 0;

        /*const*/ T* t = &(*it);
        return t;
    }


    //-------------------------------------------------------------------//
    // find_it_or_fail()                                                 //
    //-------------------------------------------------------------------//
    // This function takes the given object and determines if there is
    // a match in the vector.  It returns an iterator to the actual
    // object in the vector, if found.  Otherwise returns std::vector::end().
    //
    // This is the function you want to use most of the time
    // (or the predicate version if you are using object pointers).
    //
    // USAGE: it makes most sense to use this function if you have
    // an object with a key, other member variables, and operator<()
    // that uses the key to test for equality.  You then set up a dummy
    // "search" object with the key set to the search value, call the
    // function, and use the result to extract the additional information
    // from the object.
    //-------------------------------------------------------------------//
     typename std::vector<T>::iterator find_it_or_fail( const T& key )
     {
         typename std::vector<T>::iterator it = lower_bound_it( key );

          if ( it != inherited::end() )

              // lower_bound() does not necessarily indicate a successful search.
              // The iterator points to the object where an insertion
              // should take place.  We check that result to see if we actually
              // had an exact match.

              // NOTE: This is how the STL determines equality using only operator()<.
              // Two comparisons, ugg, but it is a nice little trick.
              if( !((*it)<key) && !(key<(*it)) )

                  return it;

          return inherited::end();
     }

     typename std::vector<T>::const_iterator find_it_or_fail( const T& key ) const
     {
         typename std::vector<T>::const_iterator it = 
             std::lower_bound( inherited::begin(), inherited::end(), key );

          if ( it != inherited::end() )
              if( !((*it)<key) && !(key<(*it)) )
                  return it;

          return inherited::end();
     }

     //-------------------------------------------------------------------//
     // find_ptr_or_fail()                                                 //
     //-------------------------------------------------------------------//
     // A variation of find_it_or_fail() that provides a pointer to result.
     //-------------------------------------------------------------------//
      T* find_ptr_or_fail( const T& key )
      {
          typename std::vector<T>::iterator it = find_it_or_fail( key );
          if ( it != inherited::end() )
              return &(*it);

          return 0;
      }

    // Same as find_it_or_fail but return an index instead of an iterator
    size_t getIndex(const T& key ) const
    {
        return find_it_or_fail(key) - inherited::begin();
    }

    //-------------------------------------------------------------------//
    // push_sorted()																		//
    //-------------------------------------------------------------------//
    // This is used to insert into a vector that always remains sorted.
    // Because we have a sorted vector, finding the insertion location
    // with std::lower_bound() is relatively cheap.
    //
    // If you have multiple insertions, consider
    // using push_unsorted() for each, then calling sort().
    //-------------------------------------------------------------------//
    void push_sorted( const T& t )
    {
        if ( !m_bSorted )
        {
            sort();
        }

        // Insert at "lower_bound" (the proper sorted location).
        inherited::insert( std::lower_bound( inherited::begin(), inherited::end(), t ), t );
    }

    // Same as push_sorted but only inserts the element if it does not exist already
    void push_sorted_unique( const T& t )
    {
        if ( !m_bSorted )
        {
            sort();
        }
        
        // iterator itr// Assumes sorted in ascending order
        //     = std::lower_bound(vec.begin(), vec.end(), e, std::greater<C>() );
        iterator pos = std::lower_bound(inherited::begin(), inherited::end(), t );
        
        if ( pos == inherited::end() || *pos != t )// If not found
            inherited::insert(pos, t);

    }

    //-------------------------------------------------------------------//
    // push_unsorted()																	//
    //-------------------------------------------------------------------//
    // This is similar to push_back(), but in addition, it sets the
    // unsorted flag.
    //-------------------------------------------------------------------//
    void push_unsorted( const T& t )
    {
        SetSorted( false );
        this->push_back(t);
    }

    size_t uniqueSize() const
    {
        if ( inherited::begin() ==  inherited::end() ) 
            return 0;

        if ( !m_bSorted )
            gsWarn<<"gsSortedVector is not sorted, uniqueSize()"
                  << "is not guaranteed to be correct.\n";

        size_t cnt = 1;
        for (const_iterator it = inherited::begin()+1; it != inherited::end(); ++it)
            if ( *(it-1) != *(it) ) ++cnt;

        return cnt;
    }

    //-------------------------------------------------------------------//
    // operator=()																	//
    //-------------------------------------------------------------------//
    // This allows us to set the gsSortedVector from a std::vector.
    //-------------------------------------------------------------------//
    gsSortedVector<T>& operator=(std::vector<T> v)
    {
        v.swap(*this);
        sort();
        return *this;
    }

    // CALLS WHERE YOU PROVIDE THE FUNCTOR OR FUNCTION POINTER
    // If you need to use a predicate sort function, ALWAYS use these methods
    // instead of the non-functor versions.
    // NOTE: UPDATE THIS when C++0x polymorphic function wrappers are available.
   template<class _Pr> inline
    void sort( _Pr pr )
    {
        if ( !m_bSorted )
        {
            std::sort( inherited::begin(), inherited::end(), pr );
            SetSorted();
        }
    }
    template<class _Pr> inline
     typename std::vector<T>::iterator lower_bound_it( const T& key, _Pr pr )
     {
         if ( !m_bSorted )
         {
             std::sort( inherited::begin(), inherited::end(), pr );
             SetSorted();
         }
         typename std::vector<T>::iterator it = std::lower_bound( inherited::begin(), inherited::end(), key, pr );
         return it;
     }
    
    template<class _Pr> inline
    /*const*/ T* lower_bound_ptr( const T& key, _Pr pr )
    {
        typename std::vector<T>::iterator it = lower_bound_it( key, pr );

        if (it==inherited::end())
            return 0;

        /*const*/ T* t = &(*it);
        return t;
    }

     template<class _Pr> inline
     void push_sorted( const T& t, _Pr pr )
     {
         if ( !m_bSorted )
         {
             std::sort( inherited::begin(), inherited::end(), pr );
             SetSorted();
         }

         // Insert at "lower_bound" (the proper sorted location).
         insert( std::lower_bound( inherited::begin(), inherited::end(), t, pr ), t );
     }

     template<class _Pr> inline
      typename std::vector<T>::iterator find_it_or_fail( const T& key, _Pr pr )
      {
          typename std::vector<T>::iterator it = lower_bound_it( key, pr );

          if ( it != inherited::end() )
              // NOTE: We have to apply this using the predicate function, be careful...
              if (!(pr((*it), key)) && !(pr(key,(*it))))
                  return it;

          return inherited::end();
      }

      template<class _Pr> inline
      T* find_ptr_or_fail( const T& key, _Pr pr )
      {
          typename std::vector<T>::iterator it = find_it_or_fail( key, pr );
          if ( it != inherited::end() )
              return &(*it);

          return 0;
      }

protected:
    bool m_bSorted;

private:

    // to do
    //void push_back( const T& t);
    
};


} //namespace gismo

