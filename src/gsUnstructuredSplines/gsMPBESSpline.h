/** @file gsMPBESSpline.h

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
#include <gsMSplines/gsMPBESBasis.h>

namespace gismo
{
template<short_t d,class T>
class gsMPBESSpline : public gsMappedSpline<d,T>
{

private:
    typedef gsMappedSingleBasis<d,T> Basis;
    typedef gsMappedSpline<d,T> Base;
public:
    /// Default empty constructor
    gsMPBESSpline() : Base() { }

    /// Construct B-Spline by basis and coefficient matrix
    gsMPBESSpline( const gsMPBESBasis<d,T> & basis, const gsMatrix<T> & coefs );

    gsMPBESSpline( gsMultiPatch<T> const & mp,int incrSmoothness = -1,int minEVDistance = -1 );

    gsMPBESSpline( gsMultiPatch<T> const & mp,std::vector<patchCorner> C0List,int incrSmoothness = -1,int minEVDistance = -1 );

    ~gsMPBESSpline() {   } //destructor

    using Base::m_mbases;
    using Base::m_global;

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

}; // class gsMPBESSpline

} // namespace gismo
