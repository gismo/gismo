/** @file gsRemapInterface.h

    @brief Provides a mapping from the patch side of geometry one to the corresponding patch side of geometry two

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Seiler, R. Schneckenleitner
    Created on: 2018-06-12
*/


#pragma once

#include <gsCore/gsAffineFunction.h>
#include <gsModeling/gsCurveFitting.h>
#include <gsModeling/gsFitting.h>

namespace gismo {

/// Provides a mapping from the patch side of geometry one to the corresponding patch side of geometry two
template <class T>
class gsRemapInterface : public gsFunction<T>
{
public:
    /// Shared pointer for gsRemapInterface
    typedef memory::shared_ptr< gsRemapInterface > Ptr;

    /// Unique pointer for gsRemapInterface
    typedef memory::unique_ptr< gsRemapInterface > uPtr;

    /// Unique pointer to domain iterator
    typedef memory::unique_ptr< gsDomainIterator<T> > domainIterUPtr;

    /// Constructor which takes a multipatch and a boundary interface, useful if the interface is fully matching
    gsRemapInterface(const gsMultiPatch<T> & mp, const gsMultiBasis<T> & basis, const boundaryInterface & bi);

    /// Constructor for the class which takes two geometries
    gsRemapInterface(const gsGeometry<T> & g1, const gsGeometry<T> & g2, const gsBasis<T>* basis1, const gsBasis<T>* basis2) : m_g1(g1), m_g2(g2), m_b1(*basis1), m_b2(*basis2)
    {
        GISMO_ASSERT(m_g1.geoDim() == m_g2.geoDim(), "The two given geometries do not have the same dimension!");
        findInterface();
        constructReparam();
        if(!m_isMatching)
            constructBreaks();
   }

    /// Helper to compute the closest point to lti on the other patch via Newton's method
    gsMatrix<T> closestPoint(const gsMatrix<T> b_null, const gsGeometry<T> & R, const gsMatrix<T> & lti);

    // rename: getPointsOnInterface() --> check eval_into, then evaluate both g1 and g2 and check equality

    const std::vector<T> getPointsOnInterface() const;

    const typename gsFunction<T>::Ptr & giveInterfaceMap() const { return m_fittedInterface; }
    //const boxSide & giveSide() const { return m_side; }

    /// Computes parametric values on the boundary from patch1 to the corresponding boundary values on patch2
    virtual void eval_into(const gsMatrix<T>& u, gsMatrix<T>& result) const;

    /// Returns arameter dimension of the domain
    virtual short_t domainDim() const { return m_g1.geoDim(); }

    /// Returns a domain iterator
    domainIterUPtr makeDomainIterator() const;

    /// @brief Returns true iff the discretization is matching.
    ///
    /// In this case, the mapping is only affine-linear.
    bool isMatching() const { return m_isMatching; }

    /// Returns the break points
    const gsMatrix<T> & breakPoints() const { return m_breakpoints; }

private:
    // flag if the interfaces are matching
    // if true then an affine map is created -> faster since no inversions etc. must be performed
    bool m_isMatching;

    // flag which says whether the orientation of the second side is fliped
    // this is important for reparameterization, especially for constructing the breakpoints
    bool m_flipSide2;

    // The geometries to consider
    const gsGeometry<T> & m_g1;
    const gsGeometry<T> & m_g2;

    // The bases to consider
    const gsBasis<T> & m_b1;
    const gsBasis<T> & m_b2;

    gsMatrix<T> m_breakpoints;
    typename gsFunction<T>::Ptr m_fittedInterface;

    // Store which boundary is the interface for patch 1 and patch 2, respectively
    patchSide m_side1;
    patchSide m_side2;

    // Store the parameter bounds for both patches
    // A single matrix has the structure [lower, upper]^T
    std::pair<gsMatrix<T>, gsMatrix<T> > m_parameterbounds;

    // Member to enrich a matrix of 1D points to a matrix of m_domain.geoDim() points
    void enrichToVector(const short_t boundarySide, const gsGeometry<T> & geo, const gsMatrix<T> & intervals, gsMatrix<T> & pts);

    // Find the interface between the two incoming patches
    void findInterface();
    void findInterface(const boundaryInterface& bi);

    // Check if the incoming patches are matching or not
    bool checkIfMatching();

    // Check if the incoming evaluation points are out of bounds because of rounding errors
    gsMatrix<T> checkIfInBound(const gsMatrix<T> & u) const;

    // Change dir direction of the parameterization of the patches
    void changeDir(const boundaryInterface & bi);

    // Constructs the reparametrization \a m_reparamInterfaceMap
    void constructReparam();

    // Cconstructs the breakpoints \a m_breakpoints
    void constructBreaks();

}; // End gsRemapInterface


} // End namespace gismo

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsRemapInterface.hpp)
#endif
