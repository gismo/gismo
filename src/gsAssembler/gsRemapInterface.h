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


template <class T>
class gsRemapInterface : public gsFunction<T>
{
public:
    /// Shared pointer for gsAffineFunction
    typedef memory::shared_ptr< gsRemapInterface > Ptr;

    /// Unique pointer for gsAffineFunction
    typedef memory::unique_ptr< gsRemapInterface > uPtr;

    /// definition of a domain iterator
    typedef memory::unique_ptr< gsDomainIterator<T> > domainIter;

    // Default constructor
    //gsRemapInterface(): m_g1(NULL), m_g2(NULL), m_b1(NULL), m_b2(NULL) {}

    // Constructor which takes a multipatch and a boundary interface, useful if the interface is fully matching
    gsRemapInterface(const gsMultiPatch<T> & mp, const gsMultiBasis<T> & basis, const boundaryInterface & bi); // in most cases they are just the other way around

    // Constructor for the class which takes two geometries
    gsRemapInterface(const gsGeometry<T> & g1, const gsGeometry<T> & g2, const gsBasis<T>* basis1, const gsBasis<T>* basis2) : m_g1(g1), m_g2(g2), m_b1(*basis1), m_b2(*basis2)
    {
        GISMO_ASSERT(m_g1.geoDim() == m_g2.geoDim(), "The two given geometries do not have the same dimension!");

        findInterface();

        // if the two patches match, maybe for testing purposes
        //m_domain.computeTopology();
        //gsInfo << "detail: " << m_domain.detail() << "\n";
    }

    // Destructor for the class
    ~gsRemapInterface() {};

    // Member for constructing the reparametrization
    // Fills in m_reparamInterfaceMap
    void constructReparam();

    // Memeber function for constructing the breakpoints
    // Fills in m_breakpoints
    void constructBreaks();

    // Helper to compute the closest point to lti on the other patch via Newton's method
    gsMatrix<T> closestPoint(const gsMatrix<T> b_null, const gsGeometry<T> & R, const gsMatrix<T> & lti);

    // rename: getPointsOnInterface() --> check eval_into, then evaluate both g1 and g2 and check equality

    const std::vector<T> getPointsOnInterface() const;

    const typename gsFunction<T>::Ptr & giveInterfaceMap() const { return m_fittedInterface; }
    //const boxSide & giveSide() const { return m_side; }

    // Member to compute parametric values on the boundary from patch1 to the corresponding boundary values on patch2
    virtual void eval_into(const gsMatrix<T>& u, gsMatrix<T>& result) const;

    // Overloaded member to get the parameter dimension of the domain
    virtual int domainDim() const { return m_g1.geoDim(); }

    domainIter makeDomainIterator() const
    {
        if (m_isMatching) return m_b1.makeDomainIterator(m_side1);

        gsTensorDomainBoundaryIterator<T> * tdi = new gsTensorDomainBoundaryIterator<T> (m_b1, m_side1);

        std::vector<T> newBreaks = getPointsOnInterface();
        gsInfo << "newBreaks: \n";
        for(index_t i = 0; i < m_breakpoints.cols(); i++)
            gsInfo << newBreaks[i] << "\t";

        gsInfo << "\n";

        // the input must be the direction which is moving
        //tdi->setBreaks(newBreaks, m_side1.direction()); -> gives the fixed direction
        // workaround: only works for 2 dimensions

        if (m_side1.direction() == 1)
            tdi->setBreaks(newBreaks, 0);
        else //m_side1.direction() == 0
            tdi->setBreaks(newBreaks, 1);


        return domainIter(tdi);
    }

    // Method to do all the necessary stuff, i.e., computing the new parametrization and computing the breakpoints for the
    // Gauss quadrature for the first geometry
    // can be moved to the constructor??
    void init()
    {
        constructReparam();

        if(!m_isMatching)
            constructBreaks();
    }

    // methods for IETIAssembler to compute the faceaverages
    bool isMatching() const { return m_isMatching; }
    const gsMatrix<T> & breakPoints() const { return m_breakpoints; }

private:
    // flag if the interfaces are matching
    // if true then an affine map is created -> faster since no inversions etc. must be performed
    bool m_isMatching;

    // The geometries to consider
    const gsGeometry<T> & m_g1;
    const gsGeometry<T> & m_g2;

    // The basis to consider
    const gsBasis<T> & m_b1;
    const gsBasis<T> & m_b2;

    gsMatrix<T> m_breakpoints;
    typename gsFunction<T>::Ptr m_fittedInterface;
    //boundaryInterface m_interface; // maybe not used anymore

    // Member to store which boundary is the interface for patch 1 and patch 2, respectively
    patchSide m_side1;
    patchSide m_side2;

    // Member to store the parameter bounds for both patches
    // A single matrix has the structure [lower, upper]^T
    std::pair<gsMatrix<T>, gsMatrix<T> > m_parameterbounds;

    // Member to enrich a matrix of 1D points to a matrix of m_domain.geoDim() points
    void enrichToVector(const int boundarySide, const gsGeometry<T> & geo, const gsMatrix<T> & intervals, gsMatrix<T> & pts);

    // Member to find the interface between the two incoming patches
    void findInterface();
    void findInterface(const boundaryInterface& bi);

    //Member to check if the incoming patches are matching or not
    bool checkIfMatching();

}; // End gsRemapInterface


} // End namespace gismo

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsRemapInterface.hpp)
#endif
