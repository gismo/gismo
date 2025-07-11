/** @file gsTensorNurbsBasis.h

    @brief Provides declaration of TensorNurbsBasis abstract interface.

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): A. Mantzaflaris
*/

#pragma once

#include <gsCore/gsRationalBasis.h>
#include <gsNurbs/gsTensorBSplineBasis.h>
#include <gsTensor/gsTensorTools.h>
#include <gsNurbs/gsNurbsBasis.h>

namespace gismo
{

/** \brief 
    A tensor product Non-Uniform Rational B-spline (NURBS) basis.

    This is the rational version of gsTensorBSplineBasis.

    \tparam d dimension of the parameter domain
    \tparam T coefficient type
    \tparam KnotVectorType  the knot vector type the underlying NURBS bases use

    \ingroup basis
    \ingroup Nurbs
*/
template<short_t d, class T>
class gsTensorNurbsBasis : public gsRationalBasis<typename gsBSplineTraits<d,T>::Basis>
{

public: 
    typedef gsKnotVector<T> KnotVectorType;

    typedef memory::unique_ptr<gsGeometry<T> > gsGeoPtr;

    /// @brief Base type
    typedef gsRationalBasis<typename gsBSplineTraits<d,T>::Basis> Base;

    /// @brief Family type
    typedef gsBSplineBasis<T>  Family_t;

    /// @brief Source basis type
    typedef typename gsBSplineTraits<d,T>::Basis Src_t;

    /// @brief Coordinate basis type
    typedef typename Src_t::Basis_t Basis_t;

    /// @brief Coefficient type
    typedef T Scalar_t;

    /// @brief Associated geometry type
    typedef typename gsBSplineTraits<d,T>::RatGeometry GeometryType;

    /// @brief Associated Boundary basis type
    typedef typename gsBSplineTraits<d-1,T>::RatBasis BoundaryBasisType;

    /// @brief Shared pointer for gsTensorNurbsBasis
    typedef memory::shared_ptr< gsTensorNurbsBasis > Ptr;

    /// @brief Unique pointer for gsTensorNurbsBasis
    typedef memory::unique_ptr< gsTensorNurbsBasis > uPtr;
    
    //typedef typename Base::iterator iterator;
    //typedef typename Base::const_iterator const_iterator;

public:

    explicit gsTensorNurbsBasis(std::vector<KnotVectorType> KV, gsMatrix<T> w)
    : Base(new Src_t(give(KV)), give(w)) { }

    // Constructors forwarded from the base class
    gsTensorNurbsBasis() : Base() { };

    gsTensorNurbsBasis( Src_t* basis ) : Base(basis) { }

    gsTensorNurbsBasis( const Src_t & basis ) : Base(basis) { }

    gsTensorNurbsBasis( Src_t* basis, gsMatrix<T> w ) : Base(basis, give(w)) { }

    gsTensorNurbsBasis(const gsTensorNurbsBasis & o) : Base(o) { }

    GISMO_CLONE_FUNCTION(gsTensorNurbsBasis)

    gsGeoPtr makeGeometry( gsMatrix<T> coefs ) const override;

public:

    /// @brief Prints the object as a string.
    std::ostream &print(std::ostream &os) const override
    {
        os << "TensorNurbsBasis: dim=" << this->dim()<< ", size="<< this->size() << ".";
        for ( unsigned i = 0; i!=d; ++i )
            os << "\n  Direction "<< i <<": "<< this->m_src->component(i).knots() <<" ";
        os << "\n";
        return os;
    }

    gsKnotVector<T> & knots (int i)
    { return m_src->knots(i); }

    const gsKnotVector<T> & knots (int i) const
    { return m_src->knots(i); }

    // knot \a k of direction \a i
    T knot(int i, int k) const
    { return m_src->knot(i, k); }

    /// The number of basis functions in the direction of the k-th parameter component
    void size_cwise(gsVector<index_t,d> & result) const
    {
        // call the function of the underlying basis
        m_src->size_cwise(result);
    }

    /// Returns the strides for all dimensions
    void stride_cwise(gsVector<index_t,d> & result) const
    {
        // call the function of the underlying basis
        m_src->stride_cwise(result);
    }

    void swapDirections(const unsigned i, const unsigned j)
    {
        gsVector<index_t, d> sz;
        size_cwise(sz);

        // First swap the weights
        swapTensorDirection(i, j, sz, m_weights);

        // Then swap the basis components
        m_src->swapDirections(i, j);
    }

    void uniformRefine_withCoefs(gsMatrix<T>& coefs, int numKnots = 1, int mul = 1, short_t const dir = -1) override
    {
        GISMO_ASSERT( coefs.rows() == this->size() && m_weights.rows() == this->size(),
                      "Invalid dimensions" );

        gsSparseMatrix<T, RowMajor> transfer;
        if (dir==-1)
        {
            m_src->uniformRefine_withTransfer(transfer, numKnots, mul);

            coefs     = transfer * ( m_weights.asDiagonal() * coefs);
            m_weights = transfer * m_weights;
            // Alternative way
            // gsBasis<T> * tmp = m_src->clone();
            // tmp->uniformRefine_withCoefs(coefs, numKnots);
            // delete tmp;
            // m_src->uniformRefine_withCoefs(m_weights, numKnots);

            // back to affine coefs
            coefs.array().colwise() /= m_weights.col(0).array();
            // equiv:
            // for (int i = 0; i < coefs.rows(); ++i)
            //    coefs.row(i) /= m_weights.at(i);
        }
        else
        {
            GISMO_ASSERT( dir >= 0 && static_cast<unsigned>(dir) < d,
                          "Invalid basis component "<< dir <<" requested for uniform refinement." );

            gsVector<index_t,d> sz;
            m_src->size_cwise(sz);
            m_src->component(dir).uniformRefine_withTransfer( transfer, numKnots, mul );

            const index_t coefs_cols = coefs.cols();
            const index_t weights_cols = m_weights.cols();

            coefs = m_weights.asDiagonal() * coefs; //<<<<-----this goes wrong!!
            swapTensorDirection(0, dir, sz, coefs);
            coefs.resize( sz[0], coefs_cols * sz.template tail<static_cast<short_t>(d-1)>().prod() );
            coefs     = transfer * coefs;

            swapTensorDirection(0, dir, sz, m_weights);
            m_weights.resize( sz[0], weights_cols * sz.template tail<static_cast<short_t>(d-1)>().prod() );
            m_weights = transfer * m_weights;

            sz[0] = coefs.rows();

            coefs.resize( sz.prod(), coefs_cols );
            m_weights.resize( sz.prod(), weights_cols );
            swapTensorDirection(0, dir, sz, coefs);
            swapTensorDirection(0, dir, sz, m_weights);

            coefs.array().colwise() /= m_weights.col(0).array();
        }
    }

    /// @brief Gives back the boundary basis at boxSide s
    gsBasis<real_t> * boundaryBasis_impl(const boxSide & s) const override;

    void matchWith(const boundaryInterface & bi, const gsBasis<T> & other,
                   gsMatrix<index_t> & bndThis, gsMatrix<index_t> & bndOther) const
    {
        this->matchWith(bi,other,bndThis,bndOther,0);
    }

    // see gsBasis for documentation
    void matchWith(const boundaryInterface & bi, const gsBasis<T> & other,
                   gsMatrix<index_t> & bndThis, gsMatrix<index_t> & bndOther, index_t offset) const override
    {
        if ( const gsTensorNurbsBasis<d,T> * _other = dynamic_cast<const gsTensorNurbsBasis<d,T> *>(&other) )
            m_src->matchWith(bi,_other->source(),bndThis,bndOther,offset);
        else if ( const gsTensorBasis<d,T> * __other = dynamic_cast<const gsTensorBasis<d,T> *>(&other) )
            m_src->matchWith(bi,*__other,bndThis,bndOther,offset);
        else
            gsWarn<<"Cannot match with "<<other<<"\n";
    }


protected:
    using Base::m_src;
    using Base::m_weights;

};

} // namespace gismo

// *****************************************************************
#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsTensorNurbsBasis.hpp)
#else
#ifdef gsTensorNurbsBasis_EXPORT
#include GISMO_HPP_HEADER(gsTensorNurbsBasis.hpp)
#undef  EXTERN_CLASS_TEMPLATE
#define EXTERN_CLASS_TEMPLATE CLASS_TEMPLATE_INST
#endif
namespace gismo
{
EXTERN_CLASS_TEMPLATE gsTensorNurbsBasis<2,real_t>;
EXTERN_CLASS_TEMPLATE gsTensorNurbsBasis<3,real_t>;
EXTERN_CLASS_TEMPLATE gsTensorNurbsBasis<4,real_t>;
}
#endif
// *****************************************************************
