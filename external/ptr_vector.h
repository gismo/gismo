// see
// http://www.codeproject.com/Articles/7351/ptr_vector-A-Container-For-Pointers

#ifndef PTR_VECTOR_H
#define PTR_VECTOR_H

/** @file
<pre>
  Copyright (c) 2004, 2005, 2006 Roland Pibinger
  All rights reserved.

  Redistribution and use in source and binary forms, with or without modification,
  are permitted provided that the following conditions are met:

  - Redistributions of source code must retain the above copyright notice, this
    list of conditions and the following disclaimer.
  - Neither the name of the copyright holders nor the names of contributors may
    be used to endorse or promote products derived from this software without
    specific prior written permission.

  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
  ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
  WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
  DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
  ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
  (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
  LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
  ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
  (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
  SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
</pre>
 *
 * @email rpbg123@yahoo.com
 *
 * $Log: $
 *
 * @version
 * @date     31-May-2004 - first version
 * @date     16-Dec-2004 - add: wrapper function templates for std::algorithms
 * @date     12-Mar-2005 - chg: iterators contain ptr_vec_iter_imp instead of deriving private(ly)
 * @date     21-Oct-2006 - chg: operator== (const ptr_vector<T>& x, const ptr_vector<T>& y) because of 'depricated' warings in VC++8.0
 * @date     25-Oct-2006 - chg: ptr_vec_iter_imp: assert in op_pp() and op_mm() corrected
 */

/*
 *
 *
 *  ptr_vector - A Container For Pointers
 *               Convenient STL-compliant vector for pointers
 * ptr_vector<T> is a wrapper for Standard vector<T*> that cuts one level of
 * indirection. In essence ptr_vector lets you treat a vector of pointers as if
 * it were a vector of values.
 *
 *     for a detailed description see
 *     http://www.codeproject.com/vcpp/stl/ptr_vecto.asp
 *
*/


#if defined (_MSC_VER) && _MSC_VER <= 1200  // MSVC++ 6.0
// disable warning C4786: identifier longer than 255 characters (identifier was truncated to '255' characters)
#   pragma warning(disable: 4786)
#endif


#include <assert.h>
#include <iterator>   //  std::back_inserter
#include <vector>
#include <functional>  // std::less
#include <algorithm>


// determine internal container to be used - std::vector<T*>  or  std::vector<void*>
#if defined (DEBUG) || defined (_DEBUG) || defined (PTV_USE_TYPED_VEC)
#  define PTV_INTERNAL_T T*     // adapt std::vector<T*>
#  define PTV_USE_TYPED_VEC
#else
#  define PTV_INTERNAL_T void*  // adapt std::vector<void*>
#endif

// workaround for MSVC++ 6.0 - keyword 'class' cannot always be used where it should
#if defined (_MSC_VER) && _MSC_VER <= 1200  // MSVC++ 6.0
#   define PTV_KW_CLASS
#else
#   define PTV_KW_CLASS class
#endif

/**
 * stdx: namespace for classes, functions, and templates related to the C++ Standard library
 */
namespace stdx {

#ifndef DOXYGEN_SHOULD_SKIP_THIS  // see Doxygen FAQ 4

// forward declaration of class templates used
template <typename T> class ptr_vector;
template <typename T> class ptr_vec_iter;
template <typename T> class const_ptr_vec_iter;
template <typename T> class ptr_iter_util;

// forward declared friend function templates (necessary according to the C++ Standard)
template <typename T> bool operator== (const ptr_vec_iter<T>& left, const ptr_vec_iter<T>& right);
template <typename T> bool operator!= (const ptr_vec_iter<T>& left, const ptr_vec_iter<T>& right);
template <typename T> bool operator<  (const ptr_vec_iter<T>& left, const ptr_vec_iter<T>& right);
template <typename T> bool operator<= (const ptr_vec_iter<T>& left, const ptr_vec_iter<T>& right);
template <typename T> bool operator>  (const ptr_vec_iter<T>& left, const ptr_vec_iter<T>& right);
template <typename T> bool operator>= (const ptr_vec_iter<T>& left, const ptr_vec_iter<T>& right);
template <typename T> bool operator== (const const_ptr_vec_iter<T>& left, const const_ptr_vec_iter<T>& right);
template <typename T> bool operator!= (const const_ptr_vec_iter<T>& left, const const_ptr_vec_iter<T>& right);
template <typename T> bool operator<  (const const_ptr_vec_iter<T>& left, const const_ptr_vec_iter<T>& right);
template <typename T> bool operator<= (const const_ptr_vec_iter<T>& left, const const_ptr_vec_iter<T>& right);
template <typename T> bool operator>  (const const_ptr_vec_iter<T>& left, const const_ptr_vec_iter<T>& right);
template <typename T> bool operator>= (const const_ptr_vec_iter<T>& left, const const_ptr_vec_iter<T>& right);


// ==================================================================
//
// ptr_iterator_tag
//
// ==================================================================

/**
 * <h2> ptr_iterator_tag </h2>
 *
 * ptr_iterator_tag is a marker interface to be inherited by pointer iterators
 * (see ptr_vec_iter and const_ptr_vec_iter). ptr_iterator_tag lets you distinguish
 * between Standard iterators and pointer iterators in cases where you need this
 * distinction.
 * The primary purpose of ptr_iterator_tag is to facilitate the implementation
 * of algorithms that take both Standard iterators and pointer iterators as arguments
 * but handle them differently. E.g. a sort algorithms can be created that
 * swap objects if Standard iterators are passed and pointers to objects
 * if pointer iterators are passed (see accompanying testcases for an example).
 *
*/
struct ptr_iterator_tag {
protected:
   ~ptr_iterator_tag() {}
};


// ==================================================================
//
// ptr_vec_iter_imp
//
// ==================================================================

/**
 * <h2> ptr_vec_iter_imp </h2>
 * contains the implentation common to both ptr_vec_iter and const_ptr_vec_iter. <br>
 * This class template is for internal use only i.e. it is not part of the
 * <code> ptr_vector</code> interface and may change in the future.
*/
template <typename PtrT>
class ptr_vec_iter_imp {
public:
   typedef std::vector<PtrT> base_vector;

   typedef std::random_access_iterator_tag       iterator_category;
   typedef typename base_vector::difference_type difference_type;
   typedef typename base_vector::size_type       size_type;

   ptr_vec_iter_imp (const base_vector* v, size_type idx): baseVec(v), index (idx) {}

   // ---------------------------------------------------------------------------
   static bool is_equal (const ptr_vec_iter_imp& left, const ptr_vec_iter_imp& right) {
      assert (check_not_null (left, right)); // 'singular' iterators are not comparable
      // assert (same_base_vector (left, right));   //#### ???
      return left.get_index() == right.get_index() && left.get_base_vec() == right.get_base_vec();
   }
   static bool is_less (const ptr_vec_iter_imp& left, const ptr_vec_iter_imp& right) {
      assert (check_not_null (left, right));
      assert (same_base_vector (left, right));
      return  left.get_index() < right.get_index();
   }

   // ---------------------------------------------------------------------------
   ptr_vec_iter_imp  op_plus  (difference_type n) const { ptr_vec_iter_imp tmp (*this); tmp.add (n); return tmp; }
   void              op_incr  (difference_type n) { add (n); }
   ptr_vec_iter_imp  op_minus (difference_type n) const { return op_plus (n * (-1)); }
   void              op_decr  (difference_type n) { op_incr (n * (-1)); }
   difference_type   op_minus (const ptr_vec_iter_imp& other) const {
      assert (check_not_null (*this, other));
      assert (this->get_base_vec() == other.get_base_vec());
      return this->get_index() - other.get_index();
   }
   void op_pp() { ++index; assert (check_bounds (*this, get_index())); }
   void op_mm() { --index; assert (check_bounds (*this, get_index())); }
   ptr_vec_iter_imp op_pp_post() { ptr_vec_iter_imp tmp (*this); op_pp(); return tmp; }
   ptr_vec_iter_imp op_mm_post() { ptr_vec_iter_imp tmp (*this); op_mm(); return tmp; }

   // ---------------------------------------------------------------------------
   typename base_vector::value_type get_ptr() const { // base_vector::value_type is a pointer!
      assert (check_bounds (*this, get_index()));
      return (*get_base_vec()) [get_index()];
   }

   const base_vector* get_base_vec() const { return baseVec; }
   size_type get_index() const { return index; }

   typename base_vector::iterator base_iter() const {
      assert (check_bounds (*this, get_index()));
      // we need to call the non-const begin() here
      base_vector* bv = const_cast<base_vector*> (get_base_vec());
      return bv->begin() + get_index();
   }

private:
   void add (difference_type n) {
      assert (check_bounds (*this, get_index() + n));
      index += n;
   }

   static bool check_bounds (const ptr_vec_iter_imp& position, size_type n) {
      bool ret = position.get_base_vec() != 0;
      if (ret) {
         ret = n <= position.get_base_vec()->size();
      }
      return ret;
   }

   static bool check_not_null (const ptr_vec_iter_imp& left, const ptr_vec_iter_imp& right) {
      return left.get_base_vec() && right.get_base_vec();
   }
   
   static bool same_base_vector (const ptr_vec_iter_imp& left, const ptr_vec_iter_imp& right) {
      return left.get_base_vec() == right.get_base_vec();
   }

   const base_vector* baseVec;
   size_type   index;
};


// ==================================================================
//
// ptr_vec_iter
//
// ==================================================================

/**
 * <h2> ptr_vec_iter </h2>
 * iterator class template for ptr_vector (typedef-ed as ptr_vector<T>::iterator)
 *  that fulfills the requirements of a Random Access Iterator;   <br>
 * iterates over pointed-to objects not pointers;  <br>
 * Similar to <code> std::vector<T>::iterator </code> this iterator represents a
 * pointer to a <b>position</b> within ptr_vector rather than a pointer to a
 * concrete object. In contrast to <code> std::vector<T>::iterator </code>
 * this iterator <b>remains valid</b> when ptr_vector expands
 * (e.g. in the progress of <code> push_back() </code>).
*/
template <typename T>
class ptr_vec_iter : public ptr_iterator_tag {
private:
   typedef stdx::ptr_vec_iter_imp<PTV_INTERNAL_T> imp;
   typedef typename imp::base_vector base_vector;

public:
   typedef T  value_type;
   typedef T& reference;
   typedef T* pointer;

   typedef typename imp::iterator_category iterator_category;
   typedef typename imp::difference_type   difference_type;
   typedef typename imp::size_type         size_type;

   ptr_vec_iter (): impl (0, size_type()) {}

private:
   typename base_vector::iterator base_iter() const { return impl.base_iter(); }
   typename base_vector::value_type get_ptr() const { return impl.get_ptr(); }
   const base_vector* get_base_vec() const { return impl.get_base_vec(); }
   size_type get_index() const { return impl.get_index(); }

   friend PTV_KW_CLASS stdx::ptr_vector<T>;
   friend PTV_KW_CLASS stdx::const_ptr_vec_iter<T>;
   friend PTV_KW_CLASS stdx::ptr_iter_util<ptr_vec_iter<T> >;

   ptr_vec_iter (const base_vector* v,   size_type idx ) : impl (v, idx)  {}
   ptr_vec_iter (const imp& pvi): impl (pvi) {}

public:
   pointer operator-> () const { return static_cast<T*> (impl.get_ptr()); }
   reference operator* () const { return * static_cast<T*> (impl.get_ptr()); }

   friend bool operator== (const ptr_vec_iter& left, const ptr_vec_iter& right) { return imp::is_equal   (left.impl, right.impl); }
   friend bool operator!= (const ptr_vec_iter& left, const ptr_vec_iter& right) { return ! imp::is_equal (left.impl, right.impl); }
   friend bool operator<  (const ptr_vec_iter& left, const ptr_vec_iter& right) { return imp::is_less    (left.impl, right.impl); }
   friend bool operator<= (const ptr_vec_iter& left, const ptr_vec_iter& right) { return ! imp::is_less  (right.impl, left.impl); }
   friend bool operator>  (const ptr_vec_iter& left, const ptr_vec_iter& right) { return imp::is_less    (right.impl, left.impl); }
   friend bool operator>= (const ptr_vec_iter& left, const ptr_vec_iter& right) { return ! imp::is_less  (left.impl, right.impl); }

   ptr_vec_iter    operator+  (difference_type n) const { return impl.op_plus (n); }
   ptr_vec_iter&   operator+= (difference_type n)  { impl.op_incr (n); return *this; }
   ptr_vec_iter    operator-  (difference_type n) const  { return impl.op_minus (n); }
   ptr_vec_iter&   operator-= (difference_type n)  { impl.op_decr (n); return *this; }
   difference_type operator-  (const ptr_vec_iter& other)  const { return impl.op_minus (other.impl); }
   reference       operator[] (difference_type n)  const{ ptr_vec_iter tmp = impl.op_plus (n); return *tmp; }

   ptr_vec_iter& operator++ ()    { impl.op_pp(); return *this; }
   ptr_vec_iter  operator++ (int) { return impl.op_pp_post(); }
   ptr_vec_iter& operator-- ()    { impl.op_mm(); return *this; }
   ptr_vec_iter  operator-- (int) { return impl.op_mm_post(); }
   
private:
   stdx::ptr_vec_iter_imp <PTV_INTERNAL_T> impl;
};


// ==================================================================
//
// const_ptr_vec_iter
//
// ==================================================================

/**
 * <h2> const_ptr_vec_iter </h2>
 * const_iterator class template for ptr_vector (ptr_vector<T>::const_iterator by typedef). <br>
 * @see ptr_vec_iter
*/
template <typename T>
class const_ptr_vec_iter : public ptr_iterator_tag {
private:
   typedef stdx::ptr_vec_iter_imp<PTV_INTERNAL_T> imp;
   typedef typename imp::base_vector base_vector;
   friend PTV_KW_CLASS stdx::ptr_vector<T>;
   
public:
   typedef T  value_type;
   typedef const T& reference;
   typedef const T* pointer;

   typedef typename imp::iterator_category iterator_category;
   typedef typename imp::difference_type   difference_type;
   typedef typename imp::size_type         size_type;

   const_ptr_vec_iter (): impl (0, size_type()) {}
   const_ptr_vec_iter (const stdx::ptr_vec_iter<T>& other) : impl (other.get_base_vec(), other.get_index()) {}

private:
   const_ptr_vec_iter (const base_vector* v,   size_type idx ) : impl (v, idx)  {}
   const_ptr_vec_iter (const imp& pvi): impl (pvi) {}

public:
   pointer operator-> () const { return static_cast<T*> (impl.get_ptr()); }
   reference operator* () const { return * static_cast<T*> (impl.get_ptr()); }

   friend bool operator== (const const_ptr_vec_iter& left, const const_ptr_vec_iter& right) { return imp::is_equal   (left.impl, right.impl); }
   friend bool operator!= (const const_ptr_vec_iter& left, const const_ptr_vec_iter& right) { return ! imp::is_equal (left.impl, right.impl); }
   friend bool operator<  (const const_ptr_vec_iter& left, const const_ptr_vec_iter& right) { return imp::is_less    (left.impl, right.impl); }
   friend bool operator<= (const const_ptr_vec_iter& left, const const_ptr_vec_iter& right) { return ! imp::is_less  (right.impl, left.impl); }
   friend bool operator>  (const const_ptr_vec_iter& left, const const_ptr_vec_iter& right) { return imp::is_less    (right.impl, left.impl); }
   friend bool operator>= (const const_ptr_vec_iter& left, const const_ptr_vec_iter& right) { return ! imp::is_less  (left.impl, right.impl); }

   const_ptr_vec_iter    operator+  (difference_type n) const { return impl.op_plus (n); }
   const_ptr_vec_iter&   operator+= (difference_type n)  { impl.op_incr (n); return *this; }
   const_ptr_vec_iter    operator-  (difference_type n) const  { return impl.op_minus (n); }
   const_ptr_vec_iter&   operator-= (difference_type n)  { impl.op_decr (n); return *this; }
   difference_type       operator-  (const const_ptr_vec_iter& other)  const { return impl.op_minus (other.impl); }
   reference             operator[] (difference_type n)  const{ const_ptr_vec_iter tmp = impl.op_plus (n); return *tmp; }

   const_ptr_vec_iter& operator++ ()    { impl.op_pp(); return *this; }
   const_ptr_vec_iter  operator++ (int) { return impl.op_pp_post(); }
   const_ptr_vec_iter& operator-- ()    { impl.op_mm(); return *this; }
   const_ptr_vec_iter  operator-- (int) { return impl.op_mm_post(); }

private:
   stdx::ptr_vec_iter_imp <PTV_INTERNAL_T> impl;
};

#endif /* DOXYGEN_SHOULD_SKIP_THIS */



// ==================================================================
//
// ptr_vector
//
// ==================================================================

/**
 * <h2> ptr_vector </h2>
 * <code> stdx::ptr_vector<T> </code> is a wrapper <code> for std::vector<T*> </code>
 * that cuts one level of indirection for iterators and member functions.
 * In essence ptr_vector lets you treat a vector of pointers as if it were a vector
 * of values.
 * <p>
 * @see
 *     http://www.codeproject.com/vcpp/stl/ptr_vecto.asp    <br>
 * for a detailed description
 */
template <typename T>
class ptr_vector {
   typedef std::vector<PTV_INTERNAL_T> base_vector;

public:
   // types:
   typedef typename base_vector::allocator_type allocator_type;
   typedef T  value_type;
   typedef T* pointer;
   typedef T& reference;
   typedef const T* const_pointer;
   typedef const T& const_reference;
   typedef typename base_vector::size_type size_type;
   typedef typename base_vector::difference_type difference_type;
   typedef stdx::ptr_vec_iter<T>       iterator;
   typedef stdx::const_ptr_vec_iter<T> const_iterator;

   // workaround for MSVC++ 6.0 non-Standard (const_)reverse_iterator definition
#if defined (_MSC_VER) && _MSC_VER <= 1200  // MSVC++ 6.0
private:
   // help VC++ 6.0 compiler to perform correct name lookup
   typedef stdx::ptr_vec_iter<T> vc60_iterator;

   typedef std::reverse_iterator
      <  vc60_iterator,
         vc60_iterator::value_type,
         vc60_iterator::reference,
         vc60_iterator::pointer,
         vc60_iterator::difference_type
      >
      vc60_reverse_iterator;
public:
   struct reverse_iterator : vc60_reverse_iterator {
      reverse_iterator() {}
      explicit reverse_iterator (const iterator& iter)     : vc60_reverse_iterator (iter) {}
      reverse_iterator (const vc60_reverse_iterator& iter) : vc60_reverse_iterator (iter) {}
      reverse_iterator (const reverse_iterator& iter)      : vc60_reverse_iterator (iter) {}
   };

private:
   typedef std::reverse_iterator
      <  const_iterator,
         const_iterator::value_type,
         const_iterator::reference,
         const_iterator::pointer,
         const_iterator::difference_type
      >
      vc60_const_reverse_iterator;
public:
   struct const_reverse_iterator : vc60_const_reverse_iterator {
      const_reverse_iterator() {}
      explicit const_reverse_iterator (const iterator& iter)        : vc60_const_reverse_iterator (iter) {}
      const_reverse_iterator (const vc60_reverse_iterator& riter)   : vc60_const_reverse_iterator (riter.base()) {}
      const_reverse_iterator (const const_reverse_iterator& criter) : vc60_const_reverse_iterator (criter.base()) {}
   };

#else
   typedef std::reverse_iterator<iterator> reverse_iterator;
   typedef std::reverse_iterator<const_iterator> const_reverse_iterator;
#endif


   // ---- construct/copy/destroy ----

   /**
   * constructs an empty ptr_vector
   */
   ptr_vector () {}

   /**
   * constructs a ptr_vector from a pair of ptr_vector<T>::iterators in the
   * range [first, last), e.g.
   * <pre>
   * ptr_vector one;
   * one.push_back (new T);
   * ptr_vector two (one.begin(), one.end());
   * </pre>
   */
   ptr_vector (iterator first, iterator last) {
      assert (same_base_vector (first, last));
      base_vector tmpVec (first.base_iter(), last.base_iter());
      baseVec.swap (tmpVec);
   }

   /**
   * constructs a ptr_vector from a pair of input iterators with value_type T*
   * in the range [first, last), e.g.
   * <pre>
   * std::list<T*> lst;
   * lst.push_back (new T);
   * ptr_vector ptv (lst.begin(), lst.end());
   * </pre>
   */
   template <typename StlPtrInputIterator>
#if defined (_MSC_VER) && _MSC_VER <= 1200  // MSVC++ 6.0
   ptr_vector (StlPtrInputIterator first, StlPtrInputIterator last, bool notUsedVc60Workaround = true) {
      base_vector tmpVec;
      for (StlPtrInputIterator start = first; start != last; ++start) {
         tmpVec.push_back (&(**start));
   }
#else
   ptr_vector (StlPtrInputIterator first, StlPtrInputIterator last) {
      base_vector tmpVec (first, last);
#endif
      assert (check_not_null (first, last));
      baseVec.swap (tmpVec);
   }

   /**
   * destroys pointers in ptr_vector; pointed-to objets are NOT destroyed!!
   */
   ~ptr_vector() {}

   // ---- iterators ----

   /** @return iterator to the first element in the ptr_vector */
   iterator begin() { return iterator (&baseVec, 0); }

   /** @return const_iterator to the first element in the ptr_vector */
   const_iterator begin() const { return const_iterator (&baseVec, 0); }

   /** @return iterator indicating one past the end */
   iterator end() { return iterator (&baseVec, size()); }

   /** @return const_iterator indicating one past the end */
   const_iterator end() const   { return const_iterator (&baseVec, size()); }

   /** @return reverse_iterator (end()) */
   reverse_iterator rbegin()             { return reverse_iterator (end()); }

   /** @return const_reverse_iterator (end()) */
   const_reverse_iterator rbegin() const { return const_reverse_iterator (end()); }

   /** @return reverse_iterator (begin()) */
   reverse_iterator rend()               { return reverse_iterator (begin()); }

   /** @return const_reverse_iterator (begin()) */
   const_reverse_iterator rend() const   { return const_reverse_iterator (begin()); }

   // ---- capacity ----

   /** @return current number of elements */
   size_type size() const      { return baseVec.size(); }

   /** @return maximal number of elements the ptr_vector can hold */
   size_type max_size() const  { return baseVec.max_size(); }

   /** @return is the ptr_vector empty */
   bool empty() const          { return baseVec.empty(); }

   /** 
     @return number of elements that the ptr_vector can hold without requiring reallocation;  <br/>
	 Note: reallocation neither copies pointed-to objects nor invalidates ptr_vector iterators!
   */
   size_type capacity() const  { return baseVec.capacity(); }

   /** ensures that capacity() is greater or equal to n */
   void reserve (size_type n)  { baseVec.reserve (n); }

   // ---- element access ----

   /** @return reference to element n */
   reference operator[] (size_type n)     { assert (check_bounds (n)); return * static_cast<pointer> (baseVec[n]); }

   /** @return const reference to element n*/
   const_reference operator[](size_type n) const  { assert (check_bounds (n)); return * static_cast<pointer> (baseVec[n]); }

   /** @return reference to element at position n */
   reference at(size_type n)              { return * static_cast<pointer> (baseVec.at(n)); }

   /** @return reference to element at position n */
   const_reference at(size_type n) const  { return * static_cast<pointer> (baseVec.at(n)); }

   /** @return reference to the first element */
   reference front() { assert (!empty()); return * static_cast<pointer> (baseVec.front()); }

   /** @return const reference to the first element */
   const_reference front() const { assert (!empty()); return * static_cast<pointer> (baseVec.front()); }

   /** @return reference to the last element */
   reference back() { assert (!empty()); return * static_cast<pointer> (baseVec.back()); }

   /** @return const reference to the last element */
   const_reference back() const  { assert (!empty()); return * static_cast<pointer> (baseVec.back()); }


   // ---- modifiers ----

   /**
   * appends element x to the ptr_vector; x must not be 0 or NULL
   */
   void push_back (pointer x)  { assert (x); baseVec.push_back (x); }

   /**
   * removes the last element from the ptr_vector;
   * @return pointer to the removed element
   */
   pointer pop_back()  {
      assert (!empty());
      typename base_vector::value_type p = baseVec.back();
      baseVec.pop_back();
      return static_cast<pointer> (p);
   }

   /**
   * inserts element x before position; x must not be 0 or NULL;
   * @return iterator to the inserted element x
   */
   iterator insert (iterator position, pointer x) {
      assert (same_base_vector (this->begin(), position));
      assert (check_iter_bounds (position));
      assert (x);
      typename base_vector::iterator iter = baseVec.insert (position.base_iter(), x);
      return iterator (&baseVec, iter - baseVec.begin());
   }

   /**
   * inserts a sequence delimited by [first, last) before position;
   */
   void insert (iterator position, iterator first, iterator last) {
      assert (same_base_vector (this->begin(), position));
      assert (same_base_vector (first, last));
      assert (check_iter_bounds (position));
      difference_type diff = last - first;
      assert (diff >= 0);
      reserve (size() + diff);
      baseVec.insert (position.base_iter(), first.base_iter(), last.base_iter());
   }

   /**
   * inserts a sequence delimited by input iterators with value_type T* in the
   * range [first, last) before position;
   */
   template <typename StlPtrInputIterator>
#if defined (_MSC_VER) && _MSC_VER <= 1200  // MSVC++ 6.0
   void insert(iterator position, StlPtrInputIterator first, StlPtrInputIterator last, bool notUsedVc60Workaround = true) {
      assert (same_base_vector (this->begin(), position));
      assert (check_iter_bounds (position));
      base_vector tmpVec;
      for (StlPtrInputIterator start = first; start != last; ++start) {
         tmpVec.push_back (&(**start));
   }
#else
   void insert(iterator position, StlPtrInputIterator first, StlPtrInputIterator last) {
      assert (same_base_vector (this->begin(), position));
      assert (check_iter_bounds (position));
      base_vector tmpVec (first, last);
#endif
      assert (check_not_null (first, last));
      reserve (size() + tmpVec.size());
      baseVec.insert (position.base_iter(), tmpVec.begin(), tmpVec.end());
   }

   /**
   * removes element denoted by position;
   * @return pointer to removed element
   */
   pointer detach (iterator position) {
      assert (same_base_vector (this->begin(), position));
      assert (check_iter_bounds (position));
      typename base_vector::value_type p = position.get_ptr();
      baseVec.erase (position.base_iter());
      return static_cast<pointer> (p);
   }

   /**
   * copies elements in the range [first, last) into the range [result, result + (last - first))
   * starting from first and proceeding to last; result must point to the begin of a range
   * large enough to hold the copied elements;
   * removes copied elements from ptr_vector  <br>
   * @see std::copy
   */
   template <typename StlPtrOutputIterator>
   void detach (iterator first, iterator last, StlPtrOutputIterator result) {
      assert (same_base_vector (this->begin(), first));
      assert (same_base_vector (this->begin(), last));
      assert (check_iter_bounds (first));
      assert (check_iter_bounds (last));
      for (iterator start = first; start != last; ++start, ++result) {
         *result = static_cast<pointer> (&(*start));
      }
      baseVec.erase(first.base_iter(), last.base_iter());
   }

   /** exchanges elements with other ptr_vector */
   void swap (ptr_vector& other)  { baseVec.swap (other.baseVec); }

   /** exchanges elements with std::vector<T*> */
   void swap (std::vector<T*>& other) {
      if (other.get_allocator() == baseVec.get_allocator()) {
         baseVec.swap (* static_cast<base_vector*>(static_cast<void*>(&other)));
      } else { // fall back
         base_vector tmpVec;
         tmpVec.reserve (other.size());
         other.reserve (baseVec.size());
         std::copy (other.begin(), other.end(), std::back_inserter (tmpVec));
         other.clear();
         baseVec.swap (tmpVec);
         for (typename base_vector::iterator iter = tmpVec.begin(); iter != tmpVec.end(); ++iter) {
            other.push_back (static_cast<pointer> (*iter));
         }
      }
      assert (check_not_null (baseVec.begin(), baseVec.end()));
   }

   /**
   * sorts ptr_vector according to operator<
   */
   void sort() {
      this->sort (std::less<T>());
   }

private:
   /**
    * helper class for comparing pointed-to objects providing pointers as input;
   */
   template <typename Predicate>
   struct deref_compare
   {
      deref_compare (const Predicate& predicate): pred (predicate) {}
      typedef typename std::vector<PTV_INTERNAL_T>::value_type base_ptr;

      bool operator()(const base_ptr left, const base_ptr right) const {
         return pred (* static_cast<const_pointer> (left), * static_cast<const_pointer> (right));
      }
      Predicate pred;
   };

public:
  /**
   * sorts ptr_vector according to Compare function (object) comp
   */
   template <typename Compare>
      void sort (Compare comp) {
      deref_compare<Compare> derefComp (comp);
      std::sort (baseVec.begin(), baseVec.end(), derefComp);

   }

private:
   bool check_bounds (difference_type idx) const {
      return idx >= 0 && size_type (idx) < size();
   }

   bool check_iter_bounds (iterator iter) const {
      return iter.get_index() <= size();
   }

   template <typename StlPtrInputIterator>
   static bool check_not_null (StlPtrInputIterator first, StlPtrInputIterator last) {
      bool found = false;
      for (StlPtrInputIterator start = first; start != last; ++ start) {
         if (*start == 0) {
            found = true;
            break;
         }
      }
      return !found;
   }
   
   static bool same_base_vector (iterator left, iterator right) {
      return left.get_base_vec() == right.get_base_vec();
   }


   // STL: ptr_vector is not a model of concept Assignable
   ptr_vector (const ptr_vector<T>&);
   ptr_vector& operator= (const ptr_vector<T>&);

   base_vector baseVec;

#if defined (_MSC_VER) && _MSC_VER <= 1200  // MSVC++ 6.0
   friend void swap (ptr_vector<T>& left,   ptr_vector<T>& right)   { left.swap (right); }
   friend void swap (std::vector<T*>& left, ptr_vector<T>& right)   { right.swap (left); }
   friend void swap (ptr_vector<T>& left,   std::vector<T*>& right) { left.swap (right); }
#else
   // see below
#endif
}; // ptr_vector


#if ! (defined (_MSC_VER) && _MSC_VER <= 1200) // not MSVC++ 6.0
// specialized algorithms:
template <typename T>
inline void swap (ptr_vector<T>& left,   ptr_vector<T>& right)   { left.swap (right); }
template <typename T>
inline void swap (std::vector<T*>& left, ptr_vector<T>& right)   { right.swap (left); }
template <typename T>
inline void swap (ptr_vector<T>& left,   std::vector<T*>& right) { left.swap (right); }
#endif


template <typename T>
inline bool operator== (const ptr_vector<T>& x, const ptr_vector<T>& y) {
  bool eq = x.size() == y.size();
  if (eq) {
#if defined (_MSC_VER) && _MSC_VER >= 1400  // VC++ 8.0
    // std::equal is depricated for VC++ 8.0; 
    eq = (!(x < y) && !(y < x));
#else
    eq = std::equal (x.begin(), x.end(), y.begin());
#endif
  }
  return eq;
}

/** compares left to right element by element */
template <typename T>
inline bool operator< (const ptr_vector<T>& left, const ptr_vector<T>& right) {
  return std::lexicographical_compare (left.begin(), left.end(), right.begin(), right.end());
}

template <typename T>
inline bool operator!=(const ptr_vector<T>& left, const ptr_vector<T>& right) { return !(left == right); }
template <typename T>
inline bool operator> (const ptr_vector<T>& left, const ptr_vector<T>& right) { return right < left; }
template <typename T>
inline bool operator>=(const ptr_vector<T>& left, const ptr_vector<T>& right) { return !(left < right); }
template <typename T>
inline bool operator<=(const ptr_vector<T>& left, const ptr_vector<T>& right) { return !(right < left); }

#ifndef DOXYGEN_SHOULD_SKIP_THIS  // see Doxygen FAQ 4

// ==================================================================
//
// wrappers for Standard algorithms
//
// ==================================================================

/**
 * <h2> ptr_iter_util </h2>
 * is a helper class template. The static member functions are
 * especially used to ease the definiton of wrappers for Standard algorithms
 * that take pointer iterators as arguments (instead of Standard iterators). <br>
 * This class template is for internal use only i.e. it is not part of the
 * <code> ptr_vector</code> interface and may change in the future.
*/

template <typename It>
class ptr_iter_util {
   typedef stdx::ptr_iter_util<It>      util;
public:
   typedef typename It::value_type      value_type;
   typedef typename It::size_type       size_type;
   typedef typename It::difference_type difference_type;

private:
   template <typename Predicate>
   struct deref_compare {
      deref_compare (Predicate predicate): pred (predicate) {}

      bool operator()(const void* left, const void* right) const {
         return pred (* static_cast<const value_type* > (left), * static_cast<const value_type* > (right));
      }
      Predicate pred;
   };

   template <typename Predicate>
   struct deref_pred {
      deref_pred (Predicate predicate): pred (predicate) {}

      bool operator()(const void* ptr) const {
         return pred (* static_cast<const value_type* > (ptr));
      }
      Predicate pred;
   };

   /** does not copy value_type objects */
   struct compare_equal {
      compare_equal (const value_type& r) : ref (r) {}

      bool operator()(const value_type& r) const {
         return r == ref;
      }
      const value_type& ref;
   };

   static void advance (It& iter, difference_type d) {
#if defined (_MSC_VER) && _MSC_VER <= 1200  // MSVC++ 6.0
      for (difference_type i = 0; i < d; ++i) {
         ++iter;
      }
#else
      std::advance (iter, d);
#endif
   }


public:
   static void ptr_swap (It first, It second) {
      std::swap (*first.base_iter(), *second.base_iter());
   }

   template <typename It2>
   static It2 swap_ranges (It first1, It last1, It2 first2) {
      difference_type d =
         std::distance (first1.base_iter(),
                        std::swap_ranges (first1.base_iter(), last1.base_iter(), first2.base_iter()));
      util::advance (first1, d);      
      
      return first1;
   }

   template <typename T>
   static It remove (It first, It last, const T& value) {
      /*util::*/compare_equal compEq (value);
      return util::remove_if (first, last, compEq);
   }

   template <typename Predicate>
   static It remove_if (It first, It last, Predicate pred) {
      It current = first,
         end     = last;
         
      while (current != end && ((current = std::find_if (current, end, pred)) != end)) {
         It next = current; ++next;
         if (next != end) {
            util::rotate (current, next, end);
         }
         --end;
      }

      return end;
   }

   static It unique (It first, It last) {
      return util::unique (first, last, std::equal_to<value_type>());
   }

   template <typename Predicate>
   static It unique (It first, It last, Predicate pred) {
      It current = first,
         end     = last;
         
      if (current != end) {
         ++current;
         while (current != end) {
            It prev = current; --prev;
            if (pred (*prev, *current)) {
               It next = current; ++next;
               if (next != end) {
                  util::rotate (current, next, end);
               }
               --end;
            } else {
               ++current;
            }
         }
      }

      return current;
   }

   static void reverse (It first, It last) {
      std::reverse (first.base_iter(), last.base_iter());
   }

   static void rotate (It first, It middle, It last) {
      std::rotate (first.base_iter(), middle.base_iter(), last.base_iter());
   }

   static void random_shuffle (It first, It last) {
      std::random_shuffle (first.base_iter(), last.base_iter());
   }

   template <typename RandomNumberGenerator>
   static void random_shuffle (It first, It last, RandomNumberGenerator& rand) {
      std::random_shuffle (first.base_iter(), last.base_iter(), rand);
   }

   template <typename Predicate>
   static It partition (It first, It last, Predicate pred) {
      /*util::*/deref_pred<Predicate> derefPred (pred);

      difference_type d =
         std::distance (first.base_iter(),
                        std::partition (first.base_iter(), last.base_iter(), derefPred));
      util::advance (first, d); 
      
      return first;
   }

   template <typename Predicate>
   static It stable_partition (It first, It last, Predicate pred) {
      /*util::*/deref_pred<Predicate> derefPred (pred);

      difference_type d =
         std::distance (first.base_iter(),
                        std::stable_partition (first.base_iter(), last.base_iter(), derefPred));
      util::advance (first, d); 

      return first;
   }

   static void sort (It first, It last) {
      util::sort (first, last, std::less<value_type>());
   }
   template <typename Compare>
   static void sort (It first, It last, Compare comp) {
      /*util::*/deref_compare<Compare> derefComp (comp);
      std::sort (first.base_iter(), last.base_iter(), derefComp);
   }

   static void stable_sort (It first, It last) {
      util::stable_sort (first, last, std::less<value_type>());
   }
   template <typename Compare>
   static void stable_sort (It first, It last, Compare comp) {
      /*util::*/deref_compare<Compare> derefComp (comp);
      std::stable_sort (first.base_iter(), last.base_iter(), derefComp);
   }

   static void partial_sort (It first, It middle, It last) {
      util::partial_sort (first, middle, last, std::less<value_type>());
   }
   template <typename Compare>
   static void partial_sort (It first, It middle, It last, Compare comp) {
       /*util::*/deref_compare<Compare> derefComp (comp);
       std::partial_sort (first.base_iter(), middle.base_iter(), last.base_iter(), derefComp);
   }

   static void nth_element (It first, It nth, It last) {
      util::nth_element (first, nth, last, std::less<value_type>());
   }
   template <typename Compare>
   static void nth_element (It first, It nth, It last, Compare comp) {
       /*util::*/deref_compare<Compare> derefComp (comp);
       std::nth_element (first.base_iter(), nth.base_iter(), last.base_iter(), derefComp);
   }

   static void inplace_merge (It first, It middle, It last) {
      util::inplace_merge (first, middle, last, std::less<value_type>());
   }
   template <typename Compare>
   static void inplace_merge (It first, It middle, It last, Compare comp) {
      /*util::*/deref_compare<Compare> derefComp (comp);
      std::inplace_merge (first.base_iter(), middle.base_iter(), last.base_iter(), derefComp);
   }

   static void push_heap (It first, It last) {
      util::push_heap (first, last, std::less<value_type>());
   }
   template <typename Compare>
   static void push_heap (It first, It last, Compare comp) {
      /*util::*/deref_compare<Compare> derefComp (comp);
      std::push_heap (first.base_iter(), last.base_iter(), derefComp);
   }

   static void pop_heap (It first, It last) {
      util::pop_heap (first, last, std::less<value_type>());
   }
   template <typename Compare>
   static void pop_heap (It first, It last, Compare comp) {
      /*util::*/deref_compare<Compare> derefComp (comp);
      std::pop_heap (first.base_iter(), last.base_iter(), derefComp);
   }

   static void make_heap (It first, It last) {
      util::make_heap (first, last, std::less<value_type>());
   }
   template <typename Compare>
   static void make_heap (It first, It last, Compare comp) {
      /*util::*/deref_compare<Compare> derefComp (comp);
      std::make_heap (first.base_iter(), last.base_iter(), derefComp);
   }

   static void sort_heap (It first, It last) {
      util::sort_heap (first, last, std::less<value_type>());
   }
   template <typename Compare>
   static void sort_heap (It first, It last, Compare comp) {
      /*util::*/deref_compare<Compare> derefComp (comp);
      std::sort_heap (first.base_iter(), last.base_iter(), derefComp);
   }

   static bool next_permutation (It first, It last) {
      return util::next_permutation (first, last, std::less<value_type>());
   }
   template <typename Compare>
   static bool next_permutation (It first, It last, Compare comp) {
      /*util::*/deref_compare<Compare> derefComp (comp);
      return std::next_permutation (first.base_iter(), last.base_iter(), derefComp);
   }

   static bool prev_permutation (It first, It last) {
      return util::prev_permutation (first, last, std::less<value_type>());
   }
   template <typename Compare>
   static bool prev_permutation (It first, It last, Compare comp) {
      /*util::*/deref_compare<Compare> derefComp (comp);
      return std::prev_permutation (first.base_iter(), last.base_iter(), derefComp);
   }
}; // ptr_iter_util
#endif /* DOXYGEN_SHOULD_SKIP_THIS */


/** converts a ptr_vector<T>::iterator to a pointer */
template <typename T>
inline T* iter2ptr (stdx::ptr_vec_iter<T> ptrIter) {
	return &(*ptrIter);
}

/** converts a ptr_vector<T>::const_iterator to a pointer to const */
template <typename T>
inline const T* iter2ptr (stdx::const_ptr_vec_iter<T> ptrIter) {
	return &(*ptrIter);
}

/** swaps underlying pointers in the container(!) of given ptr_vector<T>::iterators */
template <typename T>
inline void ptr_swap (stdx::ptr_vec_iter<T> first, stdx::ptr_vec_iter<T> second) {
   stdx::ptr_iter_util<stdx::ptr_vec_iter<T> >::ptr_swap (first, second);
}

/**
 * @see std::swap_ranges
 * @return first2 + (last1 - first1)
 */
template <typename PtrBidirectionalIterator1, typename PtrBidirectionalIterator2>
inline PtrBidirectionalIterator2 swap_ranges
   (PtrBidirectionalIterator1 first1, PtrBidirectionalIterator1 last1, PtrBidirectionalIterator2 first2) {
   return stdx::ptr_iter_util<PtrBidirectionalIterator1>::swap_ranges (first1, last1, first2);
}

/**
 * moves elements that are equal to <code>value</code> to the end of the
 * range [first, last); in contrast to std::remove all elements initially
 * within the range [first, last) remain valid after the call to stdx::remove
 * (elements are not duplicated and not destroyed). <br>
 * Complexity differs from std::remove
 *
 * @return end of resulting range
 */
template <typename PtrBidirectionalIterator, typename T>
inline PtrBidirectionalIterator remove (PtrBidirectionalIterator first, PtrBidirectionalIterator last, const T& value) {
   return stdx::ptr_iter_util<PtrBidirectionalIterator>::remove (first, last, value);
}

/**
 * moves elements for which <code>pred</code> is true to the end of the
 * range [first, last); in contrast to std::remove_if all elements initially
 * within the range [first, last) remain valid after the call to stdx::remove_if
 * (elements are not duplicated and not destroyed). <br>
 * Complexity differs from std::remove_if
 *
 * @return end of resulting range
 */
template <typename PtrBidirectionalIterator, typename Predicate>
inline PtrBidirectionalIterator remove_if (PtrBidirectionalIterator first, PtrBidirectionalIterator last, Predicate pred) {
   return stdx::ptr_iter_util<PtrBidirectionalIterator>::remove_if (first, last, pred);
}

/**
 * moves duplicates in the sorted range [first, last) to the end of the range;
 * in contrast to std::unique all elements initially
 * within the range [first, last) remain valid after the call to stdx::unique
 * (elements are not duplicated and not destroyed). <br>
 * Complexity differs from std::unique
 *
 * @return end of resulting range
 */
template <typename PtrBidirectionalIterator>
inline PtrBidirectionalIterator unique (PtrBidirectionalIterator first, PtrBidirectionalIterator last) {
   return stdx::ptr_iter_util<PtrBidirectionalIterator>::unique (first, last);
}

/**
 * moves duplicates in the range [first, last) sorted according to Predicate
 * <code>pred</code> to the end of the range;
 * in contrast to std::unique all elements initially
 * within the range [first, last) remain valid after the call to stdx::unique
 * (elements are not duplicated and not destroyed). <br>
 * Complexity differs from std::unique
 *
 * @return end of resulting range
 */
template <typename PtrBidirectionalIterator, typename Predicate>
inline PtrBidirectionalIterator unique (PtrBidirectionalIterator first, PtrBidirectionalIterator last, Predicate pred) {
   return stdx::ptr_iter_util<PtrBidirectionalIterator>::unique (first, last, pred);
}

/** @see std::reverse */
template <typename PtrBidirectionalIterator>
inline void reverse (PtrBidirectionalIterator first, PtrBidirectionalIterator last) {
   stdx::ptr_iter_util<PtrBidirectionalIterator>::reverse (first, last);
}

/** @see std::rotate */
template <typename PtrForwardIterator>
inline void rotate (PtrForwardIterator first, PtrForwardIterator middle, PtrForwardIterator last) {
   stdx::ptr_iter_util<PtrForwardIterator>::rotate (first, middle, last);
}

/** @see std::random_shuffle */
template <typename PtrRandomAccessIterator>
inline void random_shuffle (PtrRandomAccessIterator first, PtrRandomAccessIterator last) {
   stdx::ptr_iter_util<PtrRandomAccessIterator>::random_shuffle (first, last);
}

/** @see std::random_shuffle */
template <typename PtrRandomAccessIterator, typename RandomNumberGenerator>
inline void random_shuffle (PtrRandomAccessIterator first, PtrRandomAccessIterator last, RandomNumberGenerator& rand) {
   stdx::ptr_iter_util<PtrRandomAccessIterator>::random_shuffle (first, last, rand);
}

/** @see std::partition */
template <typename PtrBidirectionalIterator, typename Predicate>
inline PtrBidirectionalIterator partition (PtrBidirectionalIterator first, PtrBidirectionalIterator last, Predicate pred) {
   return stdx::ptr_iter_util<PtrBidirectionalIterator>::partition (first, last, pred);
}

/** @see std::stable_partition */
template <typename PtrBidirectionalIterator, typename Predicate>
inline PtrBidirectionalIterator stable_partition (PtrBidirectionalIterator first, PtrBidirectionalIterator last, Predicate pred) {
   return stdx::ptr_iter_util<PtrBidirectionalIterator>::stable_partition (first, last, pred);
}

/** @see std::sort */
template <typename PtrRandomAccessIterator>
inline void sort (PtrRandomAccessIterator first, PtrRandomAccessIterator last) {
   stdx::ptr_iter_util<PtrRandomAccessIterator>::sort (first, last);
}
/** @see std::sort */
template <typename PtrRandomAccessIterator, typename Compare>
inline void sort (PtrRandomAccessIterator first, PtrRandomAccessIterator last, Compare comp) {
   stdx::ptr_iter_util<PtrRandomAccessIterator>::sort (first, last, comp);
}

/** @see std::stable_sort */
template <typename PtrRandomAccessIterator>
inline void stable_sort (PtrRandomAccessIterator first, PtrRandomAccessIterator last) {
   stdx::ptr_iter_util<PtrRandomAccessIterator>::stable_sort (first, last);
}
/** @see std::stable_sort */
template <typename PtrRandomAccessIterator, typename Compare>
inline void stable_sort (PtrRandomAccessIterator first, PtrRandomAccessIterator last, Compare comp) {
   stdx::ptr_iter_util<PtrRandomAccessIterator>::stable_sort (first, last, comp);
}

/** @see std::partial_sort */
template <typename PtrRandomAccessIterator>
inline void partial_sort (PtrRandomAccessIterator first, PtrRandomAccessIterator middle, PtrRandomAccessIterator last) {
   stdx::ptr_iter_util<PtrRandomAccessIterator>::partial_sort (first, middle, last);
}

/** @see std::partial_sort */
#ifndef __GLIBCPP__
template <typename PtrRandomAccessIterator, typename Compare>
inline void partial_sort (PtrRandomAccessIterator first, PtrRandomAccessIterator middle, PtrRandomAccessIterator last, Compare comp) {
   stdx::ptr_iter_util<PtrRandomAccessIterator>::partial_sort (first, middle, last, comp);
}
#else
// lookup bug in g++ std library (overloaded function is ambiguous)
// see also: http://www.open-std.org/jtc1/sc22/wg21/docs/lwg-defects.html#225
// allegedly fixed in g++ 3.4.0
template <typename T, typename Compare>
inline void partial_sort (stdx::ptr_vec_iter<T> first, stdx::ptr_vec_iter<T> middle, stdx::ptr_vec_iter<T> last, Compare comp) {
   stdx::ptr_iter_util<stdx::ptr_vec_iter<T> >::partial_sort (first, middle, last, comp);
}
#endif

/** @see std::nth_element */
template <typename PtrRandomAccessIterator>
inline void nth_element (PtrRandomAccessIterator first, PtrRandomAccessIterator nth, PtrRandomAccessIterator last) {
   stdx::ptr_iter_util<PtrRandomAccessIterator>::nth_element (first, nth, last);
}

/** @see std::nth_element */
template <typename PtrRandomAccessIterator, typename Compare>
inline void nth_element (PtrRandomAccessIterator first, PtrRandomAccessIterator nth, PtrRandomAccessIterator last, Compare comp) {
   stdx::ptr_iter_util<PtrRandomAccessIterator>::nth_element (first, nth, last, comp);
}

/** @see std::inplace_merge */
template <typename PtrBidirectionalIterator>
inline void inplace_merge (PtrBidirectionalIterator first, PtrBidirectionalIterator middle, PtrBidirectionalIterator last) {
   stdx::ptr_iter_util<PtrBidirectionalIterator>::inplace_merge (first, middle, last);
}

/** @see std::inplace_merge */
template <typename PtrBidirectionalIterator, typename Compare>
inline void inplace_merge (PtrBidirectionalIterator first, PtrBidirectionalIterator middle,
                   PtrBidirectionalIterator last, Compare comp) {
   stdx::ptr_iter_util<PtrBidirectionalIterator>::inplace_merge (first, middle, last, comp);
}

/** @see std::push_heap */
template <typename PtrRandomAccessIterator>
inline void push_heap (PtrRandomAccessIterator first, PtrRandomAccessIterator last) {
  stdx::ptr_iter_util<PtrRandomAccessIterator>::push_heap (first, last);
}

/** @see std::push_heap */
template <typename PtrRandomAccessIterator, typename Compare>
inline void push_heap (PtrRandomAccessIterator first, PtrRandomAccessIterator last, Compare comp) {
  stdx::ptr_iter_util<PtrRandomAccessIterator>::push_heap (first, last, comp);
}

/** @see std::pop_heap */
template <typename PtrRandomAccessIterator>
inline void pop_heap (PtrRandomAccessIterator first, PtrRandomAccessIterator last) {
  stdx::ptr_iter_util<PtrRandomAccessIterator>::pop_heap (first, last);
}

/** @see std::pop_heap */
#ifndef __GLIBCPP__
template <typename PtrRandomAccessIterator, typename Compare>
inline void pop_heap (PtrRandomAccessIterator first, PtrRandomAccessIterator last, Compare comp) {
     stdx::ptr_iter_util<PtrRandomAccessIterator>::pop_heap (first, last, comp);
}
#else
// lookup bug in g++ std library (overloaded function is ambiguous)
// see also: http://www.open-std.org/jtc1/sc22/wg21/docs/lwg-defects.html#225
// allegedly fixed in g++ 3.4.0
template <typename T, typename Compare>
inline void pop_heap (stdx::ptr_vec_iter<T> first, stdx::ptr_vec_iter<T> last, Compare comp) {
     stdx::ptr_iter_util<stdx::ptr_vec_iter<T> >::pop_heap (first, last, comp);
}
#endif

/** @see std::make_heap */
template <typename PtrRandomAccessIterator>
inline void make_heap (PtrRandomAccessIterator first, PtrRandomAccessIterator last) {
  stdx::ptr_iter_util<PtrRandomAccessIterator>::make_heap (first, last);
}

/** @see std::make_heap */
#ifndef __GLIBCPP__
template <typename PtrRandomAccessIterator, typename Compare>
inline void make_heap (PtrRandomAccessIterator first, PtrRandomAccessIterator last, Compare comp) {
  stdx::ptr_iter_util<PtrRandomAccessIterator>::make_heap (first, last, comp);
}
#else
// lookup bug in g++ std library (overloaded function is ambiguous)
// see also: http://www.open-std.org/jtc1/sc22/wg21/docs/lwg-defects.html#225
// allegedly fixed in g++ 3.4.0
template <typename T, typename Compare>
inline void make_heap (stdx::ptr_vec_iter<T> first, stdx::ptr_vec_iter<T> last, Compare comp) {
  stdx::ptr_iter_util<stdx::ptr_vec_iter<T> >::make_heap (first, last, comp);
}
#endif

/** @see std::sort_heap */
template <typename PtrRandomAccessIterator>
inline void sort_heap (PtrRandomAccessIterator first, PtrRandomAccessIterator last) {
  stdx::ptr_iter_util<PtrRandomAccessIterator>::sort_heap (first, last);
}

/** @see std::sort_heap */
#ifndef __GLIBCPP__
template <typename PtrRandomAccessIterator, typename Compare>
inline void sort_heap (PtrRandomAccessIterator first, PtrRandomAccessIterator last, Compare comp) {
  stdx::ptr_iter_util<PtrRandomAccessIterator>::sort_heap (first, last, comp);
}
#else
// lookup bug in g++ std library (overloaded function is ambiguous)
// see also: http://www.open-std.org/jtc1/sc22/wg21/docs/lwg-defects.html#225
// allegedly fixed in g++ 3.4.0
template <typename T, typename Compare>
inline void sort_heap (stdx::ptr_vec_iter<T> first, stdx::ptr_vec_iter<T> last, Compare comp) {
  stdx::ptr_iter_util<stdx::ptr_vec_iter<T> >::sort_heap (first, last, comp);
}
#endif

/** @see std::next_permutation */
template <typename PtrBidirectionalIterator>
inline bool next_permutation (PtrBidirectionalIterator first, PtrBidirectionalIterator last) {
   return stdx::ptr_iter_util<PtrBidirectionalIterator>::next_permutation (first, last);
}

/** @see std::next_permutation */
template <typename PtrBidirectionalIterator, typename Compare>
inline bool next_permutation (PtrBidirectionalIterator first, PtrBidirectionalIterator last, Compare comp) {
   return stdx::ptr_iter_util<PtrBidirectionalIterator>::next_permutation (first, last, comp);
}

/** @see std::prev_permutation */
template <typename PtrBidirectionalIterator>
inline bool prev_permutation (PtrBidirectionalIterator first, PtrBidirectionalIterator last) {
   return stdx::ptr_iter_util<PtrBidirectionalIterator>::prev_permutation (first, last);
}

/** @see std::prev_permutation */
template <typename PtrBidirectionalIterator, typename Compare>
inline bool prev_permutation (PtrBidirectionalIterator first, PtrBidirectionalIterator last, Compare comp) {
   return stdx::ptr_iter_util<PtrBidirectionalIterator>::prev_permutation (first, last, comp);
}




// ==================================================================
//
// ptr_vector_owner
//
// ==================================================================

/**
 * <h2> ptr_vector_owner </h2>
 * is a <b>scope-guard</b> that takes ownership of dynamically created 
 * elememts in a ptr_vector; <br/>
 * pointed-to objects are deleted when ptr_vector_owner goes out of scope; <br/>
 * Note that ptr_vector is <b>not</b> dependent on ptr_vector_owner!
 */
template <typename T>
class ptr_vector_owner {
public:
   ptr_vector_owner (stdx::ptr_vector<T>& pv) : ptrVec (pv) {}

   ~ptr_vector_owner() {
      assert (check_no_duplicates (ptrVec));
      stdx::ptr_vector<T> tmpVec;
      tmpVec.swap (ptrVec);

      while (! tmpVec.empty()) {
         try {
            delete tmpVec.pop_back();
         } catch (...) {
            assert (false);
         }
      }
   }

private:
   stdx::ptr_vector<T>& ptrVec;

   bool check_no_duplicates (stdx::ptr_vector<T>& ptrVector) {
      std::vector<T*> tmpVec, chkVec;
      tmpVec.reserve (ptrVector.size()) , chkVec.reserve (ptrVector.size());
      ptrVector.swap (tmpVec);
      chkVec = tmpVec; // shallow copy (only pointers copied)
      ptrVector.swap (tmpVec);
      std::sort (chkVec.begin(), chkVec.end());
      typename std::vector<T*>::iterator last = std::unique (chkVec.begin(), chkVec.end());

      return last - chkVec.begin() == ptrVector.end() - ptrVector.begin();
   }

   ptr_vector_owner (const ptr_vector_owner& );
   ptr_vector_owner& operator= (ptr_vector_owner& );
}; // ptr_vector_owner


#if defined (_MSC_VER) && _MSC_VER <= 1200  // MSVC++ 6.0
#pragma message ("hint for VC++ 6.0: put 'using namespace stdx;' after your #include directive(s)")
#endif
} //namespace stdx



#endif // PTR_VECTOR_H
