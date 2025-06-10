/** @file gsRationalTHBSplineBasis.h

    @brief Provides declaration of RationalTHBSplineBasis abstract interface.

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): A. Mantzaflaris, C. Karampatzakis
*/

#pragma once

#include <gsCore/gsRationalBasis.h>
#include <gsHSplines/gsRationalTHBSpline.h>
#include <gsCore/gsForwardDeclarations.h>
#include <gsHSplines/gsTHBSplineBasis.h>


namespace gismo
{

/** \brief 
    A rational Truncated Hierarchical B-Spline basis.

    This is the rational version of gsTHBSplineBasis.

    \tparam d dimension of the parameter domain
    \tparam T coefficient type

    \ingroup basis
    \ingroup HSplines
*/
template<short_t d, class T>
class gsRationalTHBSplineBasis : public gsRationalBasis<gsTHBSplineBasis<d,T>>
{

public: 

    /// @brief Base type
    typedef gsRationalBasis<gsTHBSplineBasis<d,T>> Base;

    /// @brief Source basis type
    typedef gsTHBSplineBasis<d,T> Src_t;

    /// @brief Coordinate basis type
    typedef Src_t Basis_t;

    /// @brief Coefficient type
    typedef T Scalar_t;

    typedef typename gsHTensorBasis<d,T>::tensorBasis tensorBasis;

    /// @brief Associated geometry type
    typedef gsRationalTHBSpline<d,T> GeometryType;

    /// @brief Associated Boundary basis type
    typedef gsRationalTHBSplineBasis< (1<d?d-1:d),T> BoundaryBasisType;

    /// @brief Shared pointer for gsRationalTHBSplineBasis
    typedef memory::shared_ptr< gsRationalTHBSplineBasis > Ptr;

    /// @brief Unique pointer for gsRationalTHBSplineBasis
    typedef memory::unique_ptr< gsRationalTHBSplineBasis > uPtr;
    
    // typedef typename Base::iterator iterator;
    // typedef typename Base::const_iterator const_iterator;

public:
    // Constructors forwarded from the base class
    gsRationalTHBSplineBasis() : Base() { };

    gsRationalTHBSplineBasis(const Src_t& basis ) : Base(basis) { }

    gsRationalTHBSplineBasis( Src_t* basis, gsMatrix<T> w ) : Base(basis, give(w)) { }

    gsRationalTHBSplineBasis( gsConstantBasis<T>* basis, gsMatrix<T>)
    {GISMO_ERROR("!!"); }

    gsRationalTHBSplineBasis(const gsRationalTHBSplineBasis & o) : Base(o) { }

    GISMO_CLONE_FUNCTION(gsRationalTHBSplineBasis)
  
    GISMO_MAKE_GEOMETRY_NEW

    /// @brief Prints the object as a string.
    std::ostream &print(std::ostream &os) const
    {
        os << "Rational THB-Spline basis: dim=" << this->dim()<< ", size="<< this->size() << ".";
        // for ( unsigned i = 0; i!=d; ++i )
        //     os << "\n  Direction "<< i <<": "<< this->m_src->component(i).knots() <<" ";
        os << "\n";
        return os;
    }

    GISMO_UPTR_FUNCTION_DEF(BoundaryBasisType, boundaryBasis, boxSide const &)
    {
        typename Src_t::BoundaryBasisType::uPtr bb = m_src->boundaryBasis(n1);
        gsMatrix<index_t> ind = m_src->boundary(n1);

        gsMatrix<T> ww( ind.size(),1);
        for ( index_t i=0; i<ind.size(); ++i)
            ww(i,0) = m_weights( (ind)(i,0), 0);

        return new BoundaryBasisType(bb.release(), give(ww));// note: constructor consumes the pointer
    }
public:
    void refine_withCoefs(gsMatrix<T> & coefs, gsMatrix<T> const & boxes)
    {
        auto tmp = m_src->clone();
        coefs = m_weights.asDiagonal() * coefs;
        tmp->refine_withCoefs(coefs, boxes);
        m_src->refine_withCoefs(m_weights, boxes);
        coefs.array().colwise() /= m_weights.col(0).array();
    }

    std::vector<index_t> asElements(gsMatrix<T> const & boxes, int refExt) const
    {
        return m_src->asElements(boxes, refExt);
    }

    /// The number of basis functions in the direction of the k-th parameter component
    // void size_cwise(gsVector<index_t,d> & result) const
    // {
    //     // call the function of the underlying basis
    //     m_src->size_cwise(result);
    // }

    // /// Returns the strides for all dimensions
    // void stride_cwise(gsVector<index_t,d> & result) const
    // {
    //     // call the function of the underlying basis
    //     m_src->stride_cwise(result);
    // }

    // void swapDirections(const unsigned i, const unsigned j)
    // {
    //     gsVector<index_t, d> sz;
    //     size_cwise(sz);

    //     // First swap the weights
    //     swapTensorDirection(i, j, sz, m_weights);

    //     // Then swap the basis components
    //     m_src->swapDirections(i, j);
    // }

    // void uniformRefine_withCoefs(gsMatrix<T>& coefs, int numKnots=1, int mul=1, int dir=-1)
    // {
    //     GISMO_ASSERT( coefs.rows() == this->size() && m_weights.rows() == this->size(),
    //                   "Invalid dimensions" );

    //     gsSparseMatrix<T, RowMajor> transfer;
    //     if (dir==-1)
    //     {
    //         m_src->uniformRefine_withTransfer(transfer, numKnots, mul);

    //         coefs     = transfer * ( m_weights.asDiagonal() * coefs);
    //         m_weights = transfer * m_weights;
    //         // Alternative way
    //         // gsBasis<T> * tmp = m_src->clone();
    //         // tmp->uniformRefine_withCoefs(coefs, numKnots);
    //         // delete tmp;
    //         // m_src->uniformRefine_withCoefs(m_weights, numKnots);

    //         // back to affine coefs
    //         coefs.array().colwise() /= m_weights.col(0).array();
    //         // equiv:
    //         // for (int i = 0; i < coefs.rows(); ++i)
    //         //    coefs.row(i) /= m_weights.at(i);
    //     }
    //     else
    //     {
    //         GISMO_ASSERT( dir >= 0 && static_cast<unsigned>(dir) < d,
    //                       "Invalid basis component "<< dir <<" requested for degree elevation" );

    //         gsVector<index_t,d> sz;
    //         m_src->size_cwise(sz);
    //         m_src->component(dir).uniformRefine_withTransfer( transfer, numKnots, mul );

    //         const index_t coefs_cols = coefs.cols();
    //         const index_t weights_cols = m_weights.cols();

    //         coefs = m_weights.asDiagonal() * coefs; //<<<<-----this goes wrong!!
    //         swapTensorDirection(0, dir, sz, coefs);
    //         coefs.resize( sz[0], coefs_cols * sz.template tail<static_cast<short_t>(d-1)>().prod() );
    //         coefs     = transfer * coefs;

    //         swapTensorDirection(0, dir, sz, m_weights);
    //         m_weights.resize( sz[0], weights_cols * sz.template tail<static_cast<short_t>(d-1)>().prod() );
    //         m_weights = transfer * m_weights;

    //         sz[0] = coefs.rows();

    //         coefs.resize( sz.prod(), coefs_cols );
    //         m_weights.resize( sz.prod(), weights_cols );
    //         swapTensorDirection(0, dir, sz, coefs);
    //         swapTensorDirection(0, dir, sz, m_weights);

    //         coefs.array().colwise() /= m_weights.col(0).array();
    //     }
    // }

#ifdef __DOXYGEN__
    /// @brief Gives back the boundary basis at boxSide s
    typename BoundaryBasisType::uPtr boundaryBasis(boxSide const & s);
#endif

    // GISMO_UPTR_FUNCTION_DEF(BoundaryBasisType, boundaryBasis, boxSide const &)
    // {
    //     typename Src_t::BoundaryBasisType::uPtr bb = m_src->boundaryBasis(n1);
    //     gsMatrix<index_t> ind = m_src->boundary(n1);
        
    //     gsMatrix<T> ww( ind.size(),1);
    //     for ( index_t i=0; i<ind.size(); ++i)
    //         ww(i,0) = m_weights( (ind)(i,0), 0);
        
    //     return new BoundaryBasisType(bb.release(), give(ww));// note: constructor consumes the pointer
    // }

    // void matchWith(const boundaryInterface & bi, const gsBasis<T> & other,
    //                gsMatrix<index_t> & bndThis, gsMatrix<index_t> & bndOther) const
    // {
    //     this->matchWith(bi,other,bndThis,bndOther,0);
    // }

    // // see gsBasis for documentation
    // void matchWith(const boundaryInterface & bi, const gsBasis<T> & other,
    //                gsMatrix<index_t> & bndThis, gsMatrix<index_t> & bndOther, index_t offset) const
    // {
    //     if ( const gsRationalTHBSplineBasis<d,T> * _other = dynamic_cast<const gsRationalTHBSplineBasis<d,T> *>(&other) )
    //         m_src->matchWith(bi,_other->source(),bndThis,bndOther,offset);
    //     else if ( const gsTensorBasis<d,T> * __other = dynamic_cast<const gsTensorBasis<d,T> *>(&other) )
    //         m_src->matchWith(bi,*__other,bndThis,bndOther,offset);
    //     else
    //         gsWarn<<"Cannot match with "<<other<<"\n";
    // }


protected:
    using Base::m_src;
    using Base::m_weights;

};


} // namespace gismo

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsRationalTHBSplineBasis.hpp)
#endif