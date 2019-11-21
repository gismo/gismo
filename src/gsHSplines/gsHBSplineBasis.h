/** @file gsHBSplineBasis.h

    @brief Provides declaration of HBSplineBasis class.

    This file is part of the G+Smo library.
    
    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): G. Kiss, A. Mantzaflaris, J. Speh
*/

#pragma once

#include <gsHSplines/gsHTensorBasis.h>
#include <gsHSplines/gsHBSpline.h>

namespace gismo
{
    /** 
     * \brief
     * A hierarchical B-spline basis of parametric dimension \em d.
     *
     * See \cite Kraft1997 for the theory behind this kind of basis.
     * 
     * \tparam d the dimension of the parameter domain
     * \tparam T coefficient type
     *
     * \ingroup basis
     * \ingroup HSplines
    */ 
    
template<short_t d, class T>
class gsHBSplineBasis : public gsHTensorBasis<d,T>
{
public:
    /// @brief Associated geometry type
    typedef gsHBSpline<d,T> GeometryType;
    
    typedef typename gsHTensorBasis<d,T>::CMatrix CMatrix;

    typedef typename gsHTensorBasis<d,T>::tensorBasis tensorBasis;

    typedef typename 
    util::conditional<d==1, gsConstantBasis<T>, gsHBSplineBasis<static_cast<short_t>(d-1),T>
                      >::type BoundaryBasisType;

    /// @brief Shared pointer for gsHBSplineBasis
    typedef memory::shared_ptr< gsHBSplineBasis > Ptr;

    /// @brief Unique pointer for gsHBSplineBasis
    typedef memory::unique_ptr< gsHBSplineBasis > uPtr;

public:

    gsHBSplineBasis() { }

    /// @brief Constructor out of a tensor BSpline Basis
    gsHBSplineBasis(gsBasis<T> const&  tbasis)
        : gsHTensorBasis<d,T>(tbasis) 
    {
        // initialize(); // is done in the base constructor
    }
    
    gsHBSplineBasis( gsTensorBSplineBasis<d,T> const&  tbasis,
                     std::vector<index_t> & boxes)
        : gsHTensorBasis<d,T>(tbasis, boxes) 
    {
        // initialize(); // is done in the base constructor
    }
    
    gsHBSplineBasis( gsTensorBSplineBasis<d,T> const&  tbasis,
                     gsMatrix<T> const & boxes)
        : gsHTensorBasis<d,T>(tbasis, boxes) 
    {
        // initialize(); // is done in the base constructor
    }

#ifdef __DOXYGEN__
    /// @brief Gives back the boundary basis at boxSide s
    typename BoundaryBasisType::uPtr boundaryBasis(boxSide const & s);
#endif
    GISMO_UPTR_FUNCTION_DEF(BoundaryBasisType, boundaryBasis, boxSide const &)
    {
        return basisSlice(n1.direction(),n1.parameter());
    }

public:
    /// @brief Gives back the basis at a slice in \a dir_fixed at \a par
    BoundaryBasisType * basisSlice(index_t dir_fixed,T par ) const;

    short_t domainDim() const { return d; }
    
    void eval_into(const gsMatrix<T> & u, gsMatrix<T>& result) const;

    void deriv_into(const gsMatrix<T> & u, gsMatrix<T>& result) const;

    void deriv2_into(const gsMatrix<T> & u, gsMatrix<T>& result) const;
    
    void evalSingle_into(index_t i, const gsMatrix<T> & u, gsMatrix<T>& result) const;
    
    void derivSingle_into(index_t i, const gsMatrix<T> & u, gsMatrix<T>& result) const;
    
    void deriv2Single_into(index_t i, const gsMatrix<T> & u, gsMatrix<T>& result) const;

    GISMO_CLONE_FUNCTION(gsHBSplineBasis)

    /// @brief Prints the object as a string.
    std::ostream &print(std::ostream &os) const;
    /// @brief returns transfer matrices betweend the levels of the given hierarchical spline
    void transferbyLvl(std::vector<gsSparseMatrix<T> >& result);

    GISMO_MAKE_GEOMETRY_NEW
    
private:
    
    /// @brief Initialize the characteristic and coefficient matrices and the
    /// internal bspline representations
    void initialize();

    gsSparseMatrix<T> coarsening(const std::vector<gsSortedVector<index_t> >& old, const std::vector<gsSortedVector<index_t> >& n, const gsSparseMatrix<T,RowMajor> & transfer) const;
    gsSparseMatrix<T> coarsening_direct( const std::vector<gsSortedVector<index_t> >& old, const std::vector<gsSortedVector<index_t> >& n, const std::vector<gsSparseMatrix<T,RowMajor> >& transfer) const;

    gsSparseMatrix<T> coarsening_direct2( const std::vector<gsSortedVector<index_t> >& old, const std::vector<gsSortedVector<index_t> >& n, const std::vector<gsSparseMatrix<T,RowMajor> >& transfer) const;

    using gsHTensorBasis<d,T>::m_bases;
    using gsHTensorBasis<d,T>::m_xmatrix;
    using gsHTensorBasis<d,T>::m_xmatrix_offset;
    using gsHTensorBasis<d,T>::m_deg;
    
}; // class gsHBSplineBasis


} // namespace gismo

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsHBSplineBasis.hpp)
#endif
