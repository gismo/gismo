/** @file gsCPPInterface.h

    @brief Provides a mapping between the corresponding sides of two patches sharing an interface,
    by means of a Closest Point Projection.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s):  A. Mantzaflaris, C. Karampatzakis
    Created on: 2022-06-28
*/

#pragma once

#include <gsCore/gsMultiPatch.h>
#include <gsUtils/gsSortedVector.h>
#include <gsNurbs/gsTensorNurbsBasis.h>
#include <gsModeling/gsCurveFitting.h>

namespace gismo {

template<class T>
gsCPPInterface<T>::gsCPPInterface(const gsMultiPatch<T>   & mp,
                                  const gsMultiBasis<T>   & basis,
                                  const boundaryInterface & bi,
                                  const gsOptionList      & opt)
    : m_slaveGeom( &(mp[bi.first().patch ])),
      m_masterGeom( &(mp[bi.second().patch ])),
      m_masterBdr((mp[bi.second().patch]).boundary(bi.second())),
      m_slaveBasis(&(basis[bi.first().patch ])),
      m_masterBasis(&(basis[bi.second().patch])),
      m_boundaryInterface(bi),
      m_Tolerance(opt.getReal("Tolerance"))
{
    
    GISMO_ASSERT( m_slaveGeom->geoDim()==m_masterBdr->geoDim(), "gsCPPInterface: Dimensions do not agree." );

    m_fixedDir   = m_boundaryInterface.second().direction(); 
    // Index of the parametric direction which is normal to the boundary of the master surface
    m_fixedParam = m_boundaryInterface.second().parameter();
    // Value of the parameter in the fixed direction ( usually 0 or 1).
    // This is constant on the entire master boundary

    // Vector with the ordered indices of directions that are 'free' i.e. not the m_fixedDir
    for (index_t j=0; j<m_masterBdr->domainDim(); j++)
    {
        if (j != m_fixedDir ) { m_freeDirs.push_back(j); }
    }
}


template <class T>
void gsCPPInterface<T>::eval_into(const gsMatrix<T> & u, gsMatrix<T> & result) const
{
    GISMO_ASSERT( u.rows() == domainDim(), "gsCPPInterface::eval_into: "
        "The rows of the evaluation points must be equal to the dimension of the domain." );

    gsVector<T> masterParams;
    /* -- version with isInside check
    gsVector<bool> isInside;
    gsMatrix<T> slavePts, preIm;

    // Get the evaluation of u on the slave surface
    slavePts = m_slaveGeom->eval( u );
    //check whether slave points lie within master geometry
    m_masterGeom->locateOn(slavePts,isInside,preIm,false,m_Tolerance);

    // Pick the relevant points only, wipe out the rest
    index_t npts = 0;
    for (index_t i=u.cols()-1; i>=0; --i)
        if (isInside[i])
            u.col(i).swap(u.col(npts++));
    u.conservativeResize(u.rows(), npts);
    result.resize(u.rows(), npts);

    for (index_t i=0; i<npts; i++) // For every set( column ) of parameters
    {
        // Find the parameters of the point of the master surface boundary,
        // which is closest to the slavePoint
        m_masterBdr->closestPointTo(slavePts.col(i), masterParams, m_Tolerance);

        // masterParams are 1 less than the masters' domain dimensions
        // since only the boundary is considered in the CPP procedure
        // but we need to return domainDim, parameters thus the following loop
        for (size_t j=0; j<m_freeDirs.size(); j++)
        {
            result( m_freeDirs[j], i ) = masterParams(j);
        }
        result(  m_fixedDir, i ) = m_fixedParam; // this param is known a priori
    }

    */

    //--- closest point version
    result.resizeLike(u);
    gsVector<T> slavePoint;

    for (index_t i=0; i<u.cols(); i++) // For every set( column ) of parameters
    {
        // Get the evaluation of u on the slave surface
        slavePoint = m_slaveGeom->eval( u.col(i) );

        // Find the parameters of the point of the master surface boundary,
        // which is closest to the slavePoint
        m_masterBdr->closestPointTo( slavePoint, masterParams, m_Tolerance);

        // masterParams are 1 less than then masters domain dimensions
        // since only the boundary is considered in the CPP procedure
        // but we need to return domainDim, parameters thus the following loop
        for (size_t j=0; j<m_freeDirs.size(); j++)
        {
            result( m_freeDirs[j], i ) = masterParams(j);
        }
        result(  m_fixedDir, i ) = m_fixedParam; // this param is known a priori
    }
    //--- closest point version
}


template <class T>
typename gsDomainIterator<T>::uPtr gsCPPInterface<T>::makeDomainIterator() const
{
    gsTensorDomainBoundaryIterator<T> * tdi = new gsTensorDomainBoundaryIterator<T> (*m_slaveBasis, m_boundaryInterface.first());
    for (index_t i=0; i<domainDim(); ++i)
    {
        if (i!=m_boundaryInterface.first().direction())
            tdi->setBreaks(m_breakpoints[i],i);
    }
    return typename gsDomainIterator<T>::uPtr(tdi);
}


template <class T>
std::ostream& gsCPPInterface<T>::print(std::ostream & os) const
{
    os << "gsCPPInterface:"
       << "\n    First  (slave)  side:        " << m_boundaryInterface.first()
       << "\n    Second (master) side:        " << m_boundaryInterface.second();

    os << "\n";
    return os;
}


} // End namespace gismo
