/** @file gsCompositeIncrSmoothnessGeom.h

    @brief Provides declaration of Basis abstract interface.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): F. Buchegger
*/

#pragma once

#include <gsCore/gsMultiPatch.h>

#include <gsMSplines/gsMappedSpline.h>
#include <gsMSplines/gsMappedSingleBasis.h>
#include <gsMSplines/gsCompositeIncrSmoothnessBasis.h>

namespace gismo
{
template<short_t d,class T>
class gsCompositeIncrSmoothnessGeom : public gsMappedSpline<d,T>
{

private:
    typedef gsMappedSingleBasis<d,T> Basis;
    typedef gsMappedSpline<d,T> Base;
public:
    /// Default empty constructor
    gsCompositeIncrSmoothnessGeom() : Base() { }

    /// Construct B-Spline by basis and coefficient matrix
    gsCompositeIncrSmoothnessGeom( const gsCompositeIncrSmoothnessBasis<d,T> & basis, const gsMatrix<T> & coefs );

    gsCompositeIncrSmoothnessGeom( gsMultiPatch<T> const & mp,int incrSmoothness = -1,int minEVDistance = -1 );

    gsCompositeIncrSmoothnessGeom( gsMultiPatch<T> const & mp,std::vector<patchCorner> C0List,int incrSmoothness = -1,int minEVDistance = -1 );

    ~gsCompositeIncrSmoothnessGeom() {   } //destructor

    using Base::m_compBasis;

    void setCornerC0(patchCorner const & pc);

    void smoothCornerEdge(const patchCorner& pc,const patchSide& ps,bool updateBasis=true);

    void smoothEverything();

    virtual void uniformRefine(int numKnots=1, int mul=1);

    void refine(const index_t patch, const gsMatrix<T> &boxes);

    void refineElements(const index_t patch, std::vector<index_t> const & boxes);

    void refineElements(std::vector<std::vector<index_t> > const & boxes);

    void uniformRefineAndSmooth(int numKnots=1);

    void refineAndSmooth(const index_t patch, const gsMatrix<T> &boxes);

    void refineElementsAndSmooth(const index_t patch, std::vector<index_t> const & boxes);

}; // class gsCompositeIncrSmoothnessGeom

} // namespace gismo
