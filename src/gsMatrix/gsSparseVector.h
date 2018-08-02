/** @file gsSparseVector.h

    @brief Provides declaration of SparseVector class (wrapping Eigen)

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
    Sparse vector class, based on Eigen::SparseVector.

    \tparam T coefficient type
    \tparam _Option zero is ColMajor order.
    \tparam _Index index type

    \ingroup Matrix
*/

template<typename T, int _Options, typename _Index>
class gsSparseVector : public Eigen::SparseVector<T,_Options,_Index>
{
public:

    typedef Eigen::SparseVector<T,_Options,_Index> Base;
    
    // Type pointing to a block of the sparse vector
    typedef typename Eigen::Block<Base> Block;
    
    // Type pointing to a block view of the sparse vector
    typedef gsMatrixBlockView<Base> BlockView;
    
    /// Shared pointer for gsSparseVector
    typedef memory::shared_ptr< gsSparseVector > Ptr;
    
public:
    gsSparseVector() : Base() { }
    gsSparseVector(_Index rows) : Base(rows) { }

    /// This constructor allows constructing a gsSparseVector from
    /// Eigen expressions
    template<typename OtherDerived>
    gsSparseVector(const Eigen::EigenBase<OtherDerived>& other)  : Base(other) { }

    /// This constructor allows constructing a gsSparseVector from
    /// another sparse expression
    template<typename OtherDerived> 
    gsSparseVector(const Eigen::MatrixBase<OtherDerived>& other)  : Base(other) { } 
    
    /// This constructor allows constructing a gsSparseVector from
    /// another sparse expression
    template<typename OtherDerived> 
    gsSparseVector(const Eigen::SparseMatrixBase<OtherDerived>& other)  : Base(other) { } 

    /// This constructor allows constructing a gsSparseVector from
    /// another sparse expression
    template<typename OtherDerived> 
    gsSparseVector(const Eigen::ReturnByValue<OtherDerived>& other)  : Base(other) { } 
    
    ~gsSparseVector() { }
    
#if !EIGEN_HAS_RVALUE_REFERENCES
    // Using the assignment operators of Eigen
    // Note: using Base::operator=; is ambiguous in MSVC
#ifdef _MSC_VER
    template <class EigenExpr>
    gsSparseVector& operator= (const EigenExpr & other) 
    {
        this->Base::operator=(other);
        return *this;
    }
#else
    using Base::operator=;
#endif
    
#else
    
    // Avoid default keyword for MSVC<2013 
    // https://msdn.microsoft.com/en-us/library/hh567368.aspx
    gsSparseVector(const gsSparseVector& other) : Base(other)
    { Base::operator=(other); }

    gsSparseVector& operator= (const gsSparseVector & other)
    { Base::operator=(other); return *this; }
        
    gsSparseVector(gsSparseVector&& other)
    { operator=(std::forward<gsSparseVector>(other)); }

    gsSparseVector & operator=(gsSparseVector&& other)
    {
        this->swap(other);
        other.clear();
        return *this;
    }    

#endif

    void clear()
    {
        this->resize(0);
        this->data().squeeze();
    }

    inline T   at (_Index i ) const { return this->coeff(i); }
    inline T & at (_Index i ) { return this->coeffRef(i); }

    inline T    operator () (_Index i) const { return this->coeff(i); }
    inline T  & operator () (_Index i) { return this->coeffRef(i); }

    inline T    operator [] (_Index i) const { return this->coeff(i); }
    inline T  & operator [] (_Index i) { return this->coeffRef(i); }

    /// Clone function. Used to make a copy of the matrix
    gsSparseVector * clone() const
    { return new gsSparseVector(*this); }

}; // class gsSparseVector




} // namespace gismo
