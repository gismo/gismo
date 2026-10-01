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

    typedef typename Eigen::SparseVector<T,_Options,_Index>::InnerIterator InnerIterator;
    
    class iterator : public InnerIterator
    {
    public:
        iterator() = default;
        iterator(const gsSparseVector & sv) : InnerIterator(sv) { }

        inline T& operator[](size_t i)
        { return const_cast<T&>(*(this->m_values+i)); }
    };

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
    
    // Copy and move constructors
    gsSparseVector(const gsSparseVector& other) = default;

    gsSparseVector(gsSparseVector&& other)
    { this->swap(other); other.clear(); }

    using Base::operator=;

    gsSparseVector& operator= (const gsSparseVector & other)
    { Base::operator=(other); return *this; }
        
    gsSparseVector & operator=(gsSparseVector&& other)
    {
        this->swap(other);
        other.clear();
        return *this;
    }

    void clear()
    {
        this->resize(0);
        this->data().squeeze();
    }

    iterator begin() const { return iterator(*this); }

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

namespace Eigen { namespace internal {

template<typename T, int _Options, typename _Index>
struct traits<gismo::gsSparseVector<T,_Options,_Index> > :
    Eigen::internal::traits<Eigen::SparseVector<T,_Options,_Index> > { };

template<typename T, int _Options, typename _Index>
struct evaluator<gismo::gsSparseVector<T,_Options,_Index> > :
    evaluator<Eigen::SparseVector<T,_Options,_Index> >
{
    typedef gismo::gsSparseVector<T,_Options,_Index> XprType;
    typedef evaluator<Eigen::SparseVector<T,_Options,_Index> > Base;
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE evaluator() = default;
    EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE explicit evaluator(const XprType& m)
        : Base(static_cast<const Eigen::SparseVector<T,_Options,_Index>&>(m)) {}
};

} } // namespace Eigen::internal
