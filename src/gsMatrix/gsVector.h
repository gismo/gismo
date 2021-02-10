/** @file gsVector.h

    @brief Provides declaration of Vector class.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris
*/

# pragma once

// Assumes that Eigen library has been already included

namespace gismo
{

/** @brief
    A vector with arbitrary coefficient type and fixed or dynamic size.

    This class is based on Eigen::Matrix from the Eigen
    linear algebra library. Most operations from Eigen are supported
    on a gsVector. See therefore also the Eigen documentation,
    http://eigen.tuxfamily.org/dox/.

    \tparam T coefficient type
    \tparam _Rows number of rows: an integer or \c Dynamic

    \ingroup Matrix
*/
template<class T, int _Rows, int _Options>
class gsVector : public gsMatrix<T, _Rows, 1, _Options>
//class gsVector : public gsMatrix<T, _Rows, (_Rows!=-1 ? 1 : -1), _Options>
{
public:
    typedef gsMatrix<T,_Rows,1,_Options> gsBase;
    //typedef gsMatrix<T,_Rows,(_Rows!=-1 ? 1 : -1), _Options> gsBase;

    // Base is the single-column dense matrix class of Eigen
    typedef typename gsBase::Base Base;

    // Self type
    typedef gsVector<T,_Rows, _Options> Self;

    typedef typename Eigen::aligned_allocator<Self> aalloc;

    // The type of the coefficients of the matrix
    typedef T Scalar_t;

    // Type pointing to a block of the vector
    typedef typename Eigen::Block<Base> Block;

    // Type pointing to a block view of the vector
    typedef gsMatrixBlockView<Base> BlockView;

    /// Shared pointer for gsVector
    typedef memory::shared_ptr< gsVector > Ptr;

    /// Unique pointer for gsVector
    typedef memory::unique_ptr< gsVector > uPtr;

    // Type for copying a vector as a permutation matrix
    typedef Eigen::PermutationMatrix<_Rows,Base::SizeAtCompileTime,index_t> Permutation;

    // Type for treating a vector as a permutation matrix
    typedef Eigen::PermutationWrapper<Base> PermutationWrap;

    typedef Eigen::Ref<Base> Ref;

    typedef const Eigen::Ref<const Base> ConstRef;

    // Type for a vector of dimension one less
    typedef gsMatrix< T, ChangeDim<_Rows, -1>::D, ColMajor> Projection_t;

public:

    typedef T * iterator;

    typedef const T * const_iterator;

    T * begin()
    { return this->data(); }

    const T * begin() const
    { return this->data(); }

    T * end()
    { return this->data() + this->size(); }

    const T * end() const
    { return this->data() + this->size(); }

public:

    gsVector() ;

    gsVector(const Base& a) ;

    // implicitly deleted in C++11
    //gsVector(const gsVector& a) : gsBase(a) { }

    explicit gsVector(index_t dimension) ;

    inline operator Ref () { return Ref(*this); }

    inline operator const ConstRef () { return ConstRef(*this); }

    void clear() { this->resize(0); }

    // This constructor allows constructing a gsVector from Eigen expressions
    template<typename OtherDerived>
    gsVector(const Eigen::EigenBase<OtherDerived>& other) : gsBase(other) { }

    // This constructor allows constructing a gsVector from Eigen expressions
    template<typename OtherDerived>
    gsVector(const Eigen::MatrixBase<OtherDerived>& other) : gsBase(other) { }

    // This constructor allows constructing a gsVector from Eigen expressions
    template<typename OtherDerived>
    gsVector(const Eigen::ReturnByValue<OtherDerived>& other) : gsBase(other) { }

    static gsVector<T,2> vec( T x, T y)
    {
        return typename gsVector<T,2>::Base(x, y);
    }

    static gsVector<T,3> vec( T x, T y, T z)
    {
        return typename gsVector<T,3>::Base(x, y, z);
    }

/*
    // Using the assignment operators of Eigen
    // Note: using Base::operator=; is ambiguous in MSVC
#ifdef _MSC_VER
    template <class EigenExpr>
    gsVector& operator= (const EigenExpr & other)
    {
        this->Base::operator=(other);
        return *this;
    }
#else
    using Base::operator=;
#endif
*/
#if !EIGEN_HAS_RVALUE_REFERENCES
    gsVector & operator=(typename Eigen::internal::conditional<
                         -1==_Rows,gsVector, const gsVector &>::type other)
    {
        if (-1==_Rows)
            this->swap(other);
        else
            this->Base::operator=(other);
        return *this;
    }
#endif

    /// \brief Returns the \a i-th element of the vector
    inline T at(index_t i) const { return *(this->data()+i);}

    /// \brief Returns the \a i-th element of the vector
    inline T & at(index_t i)     { return *(this->data()+i);}

    /// \brief Returns the last (bottom) element of the vector
    inline T last() const { return *(this->data()+this->size()-1);}

    /// \brief Returns the last (bottom) element of the vector
    inline T & last() { return *(this->data()+this->size()-1);}

    /// Return a row-block view of the vector with \a rowSizes
    BlockView blockView(const gsVector<index_t> & rowSizes)
    {
        return BlockView(*this, rowSizes);
    }

    PermutationWrap asPermutation() const { return PermutationWrap(*this);}

    /// Removes row \a i from the vector. After the operation the
    /// vector has size one less.
    void removeElement(const index_t i )
    {
        GISMO_ASSERT( i < this->size(), "Invalid vector element." );
        const T * ce = this->data() + this->size();
        for ( T * c = this->data()+i+1; c!= ce; ++c ) *(c-1) = *c;
        this->conservativeResize(this->size()-1,Eigen::NoChange);
    }

}; // class gsVector




/** @brief A fixed-size, statically allocated 3D vector.

    \tparam T coefficient type

    \ingroup Matrix
*/
template<class T>
class gsVector3d : public Eigen::Matrix<T,3,1>
{
public:
    typedef T scalar_t;
    typedef Eigen::Matrix<T,3,1> Base ;

    /// Shared pointer for gsVector3d
    typedef memory::shared_ptr< gsVector3d > Ptr;

public:

    gsVector3d();

    gsVector3d(scalar_t x, scalar_t y, scalar_t z = 0 );

    // implcitly declared deleted in C++11
    //gsVector3d(const gsVector3d& a) : Base(a) { }

    gsVector3d(const Base& a) ;

    /// This constructor allows constructing a gsVector3d from Eigen expressions
    template<typename OtherDerived>
    gsVector3d(const Eigen::MatrixBase<OtherDerived>& other) : Base(other) { }

    T angle(const gsVector3d<T> & other)
    {
        return math::acos(this->normalized().dot(other.normalized()));
    }

    /*
    // Using the assignment operators of Eigen
    // Note: using Base::operator=; is ambiguous in MSVC
#ifdef _MSC_VER
    template <class EigenExpr>
    gsVector3d& operator= (const EigenExpr & other)
    {
        this->Base::operator=(other);
        return *this;
    }
#else
    using Base::operator=;
#endif
    */

    // implicitly deleted in C++11
    gsVector3d & operator=(const gsVector3d & other)
    {
        this->Base::operator=(other);
        return *this;
    }

    inline T   at (index_t i) const { return (*this)(i,0); }
    inline T & at (index_t i)       { return (*this)(i,0); }
    inline T   x () const   { return (*this)(0); }
    inline T & x ()         { return (*this)(0); }
    inline T   y () const   { return (*this)(1); }
    inline T & y ()         { return (*this)(1); }
    inline T   z () const   { return (*this)(2); }
    inline T & z ()         { return (*this)(2); }

}; // class gsVector3d


template<class T, int _Rows, int _Options> inline
gsVector<T,_Rows,_Options>::gsVector() : gsBase() { }

template<class T, int _Rows, int _Options> inline
gsVector<T,_Rows,_Options>::gsVector(const Base& a): gsBase(a) { }

template<class T, int _Rows, int _Options> inline
gsVector<T,_Rows,_Options>::gsVector(index_t dimension): gsBase(dimension,1) { }

template<class T> inline
gsVector3d<T>::gsVector3d() : Base() { }

template<class T> inline
gsVector3d<T>::gsVector3d(scalar_t x, scalar_t y,scalar_t z )
        {
            (*this)(0,0)=x;
            (*this)(1,0)=y;
            (*this)(2,0)=z;
        }

template<class T> inline
gsVector3d<T>::gsVector3d(const Base& a): Base(a) { }

// template<class T>
// template<typename OtherDerived> inline
// gsVector3d<T>::gsVector3d(const Eigen::MatrixBase<OtherDerived>& other) : Base(other) { }


// template<class T> inline
//     inline T   gsVector3d<T>::x () const { return (*this)(0); }
// template<class T> inline
//     inline T & gsVector3d<T>::x () { return (*this)(0); }
// template<class T> inline
//     inline T   gsVector3d<T>::y () const { return (*this)(1); }
// template<class T> inline
//     inline T & gsVector3d<T>::y () { return (*this)(1); }
// template<class T> inline
//     inline T   gsVector3d<T>::z () const { return (*this)(2); }
// template<class T> inline
//     inline T & gsVector3d<T>::z () { return (*this)(2); }




} // namespace gismo
