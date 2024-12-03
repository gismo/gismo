/** @file gsComposedBasis.h

    @brief

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s):
        H.M. Verhelst   (2019-..., TU Delft)
        A. Mantzaflaris (2019-..., Inria)
*/

#pragma once

#include <gsCore/gsBasis.h>
#include <gsCore/gsComposedGeometry.h>
#include <gsCore/gsFunction.h>


namespace gismo
{

template <class T>
class gsComposedBasis : public gsBasis<T>
{

    // /// Geometry Type
    typedef gsComposedGeometry<T> GeometryType;

    GISMO_OVERRIDE_CLONE_FUNCTION(gsComposedBasis)

    GISMO_OVERRIDE_MAKE_GEOMETRY_NEW

    typedef typename gsBasis<T>::domainIter domainIter;

    typedef T ScalarType;

public:
    typedef memory::shared_ptr< gsComposedBasis > Ptr;
    typedef memory::unique_ptr< gsComposedBasis > uPtr;

    typedef gsBasis<T>      BasisT;
    typedef gsFunction<T>   CompositionT;

public:

    gsComposedBasis();

    gsComposedBasis(const CompositionT * composition, const BasisT * basis);

    gsComposedBasis(const CompositionT & composition, const BasisT & basis);

    gsComposedBasis(typename CompositionT::Ptr composition,
                    typename BasisT::Ptr basis);

    short_t domainDim() const override;
    short_t targetDim() const override;

    short_t maxDegree() const override;

    gsMatrix<T> support() const override;

    gsMatrix<T> support(const index_t & i) const override;

    void active_into(const gsMatrix<T> & u, gsMatrix<index_t>& result) const override;

    // void evalAllDers_into(const gsMatrix<T> & u, int n,
    //                         std::vector<gsMatrix<T> >& result,
    //                         bool sameElement) const
    // {
    //     gsMatrix<T> coords = m_composition->eval(u);
    //     this->_applyBounds(coords);
    //     m_basis->evalAllDers_into(coords,n,result,sameElement);
    // }

    void eval_into(const gsMatrix<T>& u, gsMatrix<T>& result) const override;

    void evalSingle_into(index_t i, const gsMatrix<T>& u, gsMatrix<T>& result) const override;

    void deriv_into(const gsMatrix<T>& u, gsMatrix<T>& result) const override;

    void derivSingle_into(index_t i, const gsMatrix<T>& u, gsMatrix<T>& result) const override;

    void deriv2_into(const gsMatrix<T>& u, gsMatrix<T>& result) const override;

    // void control_deriv_into(const gsMatrix<T> & points, gsMatrix<T> & result)
    // {
    //     // The number of rows is the target dimension times the number of controls
    //     // The number of cols is the number of points
    //     result.resize(targetDim()*m_composition->nControls(),points.cols());

    //     // Pre-compute the coordinates of the composition, the derivatives of G and the derivatives of the composition
    //     gsMatrix<T> c, dc, dG;
    //     m_composition->eval_into(points,c);
    //     m_composition->control_deriv_into(points,dc);
    //     m_geom->deriv_into(c,dG);

    //     // Store some sizes
    //     index_t nControls = m_composition->nControls();
    //     index_t dd = m_geom->domainDim();
    //     index_t td = m_geom->targetDim();

    //     // Loop over the points
    //     for (index_t k=0; k!=points.cols(); k++)
    //     {
    //         // We need to compute dG/dpi = dG/dc * dc/dpi
    //         gsAsMatrix<T> DG = result.reshapeCol(k,nControls,td);
    //         DG = dc.reshapeCol(k,nControls,dd) * dG.reshapeCol(k,dd,td);
    //     }
    // }

    // memory::unique_ptr<gsGeometry<T> > makeGeometry(gsMatrix<T> coefs) const override
    // {
    //     return memory::unique_ptr<gsGeometry<T> >(new gsGeometry<T>(*this, give(coefs)));
    //     // GISMO_NO_IMPLEMENTATION;
    // }

    ////// Pass throughs of basis

    /// See \ref gsBasis for documentation
    short_t degree(short_t i) const override;

    /// See \ref gsBasis for documentation
    gsMatrix<index_t> boundaryOffset(boxSide const & s, index_t offset) const override;

    /// See \ref gsBasis for documentation
    void matchWith(const boundaryInterface & bi, const gsBasis<T> & other,
                    gsMatrix<index_t> & bndThis, gsMatrix<index_t> & bndOther, index_t offset = 0) const override;

    /// See \ref gsBasis for documentation
    domainIter makeDomainIterator() const override;

    /// See \ref gsBasis for documentation
    virtual domainIter makeDomainIterator(const boxSide & s) const override;

    /// See \ref gsBasis for documentation
    std::string detail() const override;

    /// See \ref gsBasis for documentation
    size_t numElements(boxSide const & s = 0) const override;

    /// See \ref gsBasis for documentation
    index_t size() const override;

    /// See \ref gsBasis for documentation
    void anchors_into(gsMatrix<T> & result) const override;

    /// See \ref gsBasis for documentation
    void connectivity(const gsMatrix<T> & nodes, gsMesh<T> & mesh) const override;

    /// See \ref gsBasis for documentation
    void uniformRefine(int numKnots = 1, int mul=1, int dir=-1) override;

    /// See \ref gsBasis for documentation
    void uniformRefine_withCoefs(gsMatrix<T>& coefs, int numKnots = 1, int mul = 1, int dir=-1) override;

    // See \ref gsBasis for documentation
    void degreeElevate(short_t const & i = 1, short_t const dir = -1) override;


    //// NEEDS TO EVALUATE THE INVERSE MAP, NOT m_composition!!
    void mapMesh(gsMesh<T> & mesh) const;


    /// Return the composition
    const CompositionT & composition() const;
          CompositionT & composition()      ;
    /// Return the basis
    const BasisT & basis() const;

    /// See \ref gsBasis for the documentation of this function
    std::ostream &print(std::ostream &os) const override;

private:
    /// Applies the bounds to the coordinates
    void _applyBounds(gsMatrix<T> & coords) const;

protected:
    typename CompositionT::Ptr   m_composition;
    typename BasisT::Ptr         m_basis;

}; // class gsComposedBasis

} // namespace gismo

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsComposedBasis.hpp)
#endif

