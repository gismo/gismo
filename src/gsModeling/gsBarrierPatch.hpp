/** @file gsBarrierPatch.hpp

@brief This software facilitates the creation of analysis-suitable
    parameterizations from given boundary representations. Serving as a
    reference implementation, it embodies the methods and concepts detailed
    in Ye Ji's doctoral research.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): Ye Ji, H.M. Verhelst
*/

#pragma once

#include "gsModeling/gsBarrierCore.h"
#include "gsModeling/gsBarrierCore.hpp"

namespace gismo
{
/**
 * \brief Computes a patch parametrization given a set of boundary geometries.
 * Parametrization is not guaranteed to be non-singular. Works for planar surfaces and volumes.
 *
 * \tparam d domain dimension
 * \tparam T Coefficient type
 *
 * \param bRep Input boundary representation
 * \param initialMethod Specify initialization method
 * \param filename Name of input data file
 *
 * \ingroup Modeling
*/

/// Constructs the object using a given multi-patch and a degree of freedom mapper.
template<short_t d, typename T>
gsBarrierPatch<d, T>::gsBarrierPatch(const gsMultiPatch<T> &mp,
                                     const gsDofMapper &mapper)
    : m_mp(mp), m_mb(mp)
{
  // Input data sanity check
  if (mp.empty() || !mapper.freeSize()) {
    throw std::invalid_argument("Invalid input to gsBarrierPatch constructor");
  }

  // Set the default options
  this->defaultOptions();

  // Set the mapper
  this->setMapper(mapper);
}

/// Constructs the object using a given multi-patch and an optional patch-wise flag (default true).
template<short_t d, typename T>
gsBarrierPatch<d, T>::gsBarrierPatch(const gsMultiPatch<T> &mp, bool patchWise)
    : m_mp(mp), m_mb(mp)
{
  // Input data sanity check
  if (mp.empty()) {
    throw std::invalid_argument("Invalid multi-patch input to gsBarrierPatch constructor");
  }

  // Set the default options
  this->defaultOptions();

  // Depending on patchWise flag, choose the appropriate mapper method
  if (patchWise) {
    // make mapper that also fixes interfaces
    m_freeInterface = 0;
  } else {
    // make mapper that free interfaces
    _makeMapperGlobalPatches();
  }
}

// // vertex constructor
// template<short_t d, typename T>
// gsBarrierPatch<d, T>::gsBarrierPatch( const gsMultiPatch<T> &mp, std::vector<patchCorner> corners, index_t offset)
// {
//     this->defaultOptions();
//     // make mapper that fixes a ring around vertices
//     // _makeMapper();
// }

// // interface constructor
// template<short_t d, typename T>
// gsBarrierPatch<d, T>::gsBarrierPatch( const gsMultiPatch<T> &mp, std::vector<boundaryInterface> ifaces, index_t offset)
// {
//     this->defaultOptions();
//     // make mapper that fixes an offset around sides
//     // _makeMapper();
// }

template<short_t d, typename T>
void gsBarrierPatch<d, T>::compute()
{
  // Preprocessing: Scale the computational domain to [0,1]^d for better
  // numerical stability and to have consistent convergence criteria.
  const double boxsize = 1.0;
  gsMatrix<T> boundingBox;
  m_mp.boundingBox(boundingBox);
  gsVector<T, d> bbmin = boundingBox.col(0);
  gsVector<T, d> bbmax = boundingBox.col(1);
  T maxside = (bbmax-bbmin).maxCoeff();

  gsVector<T, d> scaleFactor;
  scaleFactor.setConstant(boxsize/maxside);

  // Apply transformation (translation and scaling) to the patches
  for (auto &patch:m_mp) {
    patch->translate(-bbmin);
    patch->scale(scaleFactor);
  }

  // Start parameterization construction.
  gsStopwatch timer;
  switch (m_freeInterface)
  {
    case 0:
    {
      // Construct analysis-suitable parameterization piecewise.
      for (auto iptch = 0; iptch < m_mp.nPatches(); ++iptch)
      {
        gsInfo << "Parameterizing " << iptch
               << "-th patch. Total number of patches is "
               << m_mp.nPatches() << ".\n";

        gsDofMapper currMapper = _makeMapperOnePatch(m_mp.patch(iptch));
        gsMultiPatch<T> optCurrPatch = gsBarrierCore<d, T>::compute(m_mp.patch(iptch),
                                                                    currMapper,
                                                                    m_options);

        m_mp.patch(iptch).setCoefs(optCurrPatch.patch(0).coefs());
      }
      break;
    }
    case 1:
    {
      // Construct analysis-suitable parameterization with moving interfaces.
      m_mp = gsBarrierCore<d, T>::compute(m_mp, m_mapper, m_options);
      break;
    }
  }

  gsInfo << "\nParameterization construction completed. Running time: " << timer.stop() << ".\n";

  // Restore scale to the original size of the input model.
  for (auto &isf:scaleFactor) {
    isf = 1./isf;
  }

  for (auto &ptch:m_mp) {
    ptch->scale(scaleFactor);
    ptch->translate(bbmin);
  }
}

//template<short_t d, typename T>
//    void gsBarrierPatch<d, T>::compute()
//    {
//        // preprocessing: scale the computational domain to [0,10]^d for better
//        // numerical stability and to have a consistent convergence criteria
//        const double boxsize = 1.0;
//        gsMatrix<T> boundingBox;
//        m_mp.boundingBox(boundingBox);
//        gsVector<T, d> bbmin = boundingBox.col(0);
//        gsVector<T, d> bbmax = boundingBox.col(1);
//        real_t maxside = (bbmax-bbmin).maxCoeff();
//
//        gsVector<T, d> scaleFactor;
//        scaleFactor.setConstant(boxsize/maxside);
//        for (auto &ptch:m_mp) {
//            ptch->translate(-bbmin);
//            ptch->scale(scaleFactor);
//        }
//
//        // start parameterization construction
//        gsStopwatch timer;
//        switch (m_freeInterface)
//        {
//            case 0:
//            {
//                // construct analysis-suitable parameterization piecewisely
//                for (auto iptch = 0; iptch < m_mp.nPatches(); ++iptch)
//                {
//                    gsInfo << "I am parameterizing " << iptch
//                           << "-th patch, total number of patches is "
//                           << m_mp.nPatches() << "\n";
//                    gsDofMapper currMapper = _makeMapperOnePatch(
//                            m_mp.patch(iptch));
//                    gsMultiPatch<T> optCurrPatch = gsBarrierCore<d, T>::compute(
//                            m_mp.patch(iptch), currMapper, m_options);
//                    m_mp.patch(iptch).setCoefs(optCurrPatch.patch(0).coefs());
//                }
//                break;
//            }
//            case 1:
//            {
//                // construct analysis-suitable parameterization with moving interfaces
//                m_mp = gsBarrierCore<d, T>::compute(m_mp, m_mapper, m_options);
//                break;
//            }
//        }
//
//        real_t runtime = timer.stop();
//        gsInfo  << "\n"
//                << "Multi-patch parameterization construction completed! Running time is "
//                << runtime << "\n";
//
//        // restore scale to the original size of the input model
//        for (auto &isf:scaleFactor)
//            isf = 1./isf;
//        for (auto &ptch:m_mp) {
//            ptch->scale(scaleFactor);
//            ptch->translate(bbmin);
//        }
//    }

    /// This function sets the default options of the BarrierPatch object by
    // calling the defaultOptions method of gsBarrierCore.
    template<short_t d, typename T>
    void gsBarrierPatch<d, T>::defaultOptions()
    {
        m_options = gsBarrierCore<d, T>::defaultOptions();
    }

/// Creates a mapper.
template<short_t d, typename T>
void gsBarrierPatch<d, T>::_makeMapper()
{
  // Initiate the mapper for the design Degrees of Freedom (DoFs)
  m_mapper.init(m_mb, m_mp.targetDim());

  // For each patch, we mark the boundary in each target dimension
  for (size_t iptch = 0; iptch != m_mp.nPatches(); iptch++)
  {
    gsMatrix<index_t> idx = m_mp.basis(iptch).allBoundary();
    for (index_t idim = 0; idim != m_mp.targetDim(); idim++)
    {
      m_mapper.markBoundary(iptch, idx, idim);
    }
  }

  // Finalize the mapper after marking all boundaries
  m_mapper.finalize();

  logMapperInformation();
}

/// Creates a mapper for local patches.
template<short_t d, typename T>
void gsBarrierPatch<d, T>::_makeMapperLocalPatches()
{
  // Initiate the mapper for the design Degrees of Freedom (DoFs)
  // In this method, we set all the inner control points as optimization variables
  // However, it is also possible to set only a part of them as optimization variables
  m_mapper.init(m_mb, m_mp.targetDim());

  // Iterate over each patch to mark the boundary in each target dimension
  for (gsMultiPatch<>::const_biterator bit = m_mp.bBegin();
       bit != m_mp.bEnd(); ++bit)
  {
    gsMatrix<index_t> idx = m_mb.basis(bit->patch).allBoundary();
    for (index_t idim = 0; idim != m_mp.targetDim(); ++idim)
    {
      m_mapper.markBoundary(bit->patch, idx, idim);
    }
  }

  // Finalize the mapper after marking all boundaries
  m_mapper.finalize();

  logMapperInformation();
}

template<short_t d, typename T>
void gsBarrierPatch<d, T>::_makeMapperGlobalPatches()
{
  // Initiate the mapper for the design Degrees of Freedom (DoFs)
  // We're setting all the inner control points as optimization variables
  // However, it's also possible to set only a part of them as optimization variables
  m_mapper.init(m_mb, m_mp.targetDim());

  // Ensure C^0 continuity at the interface by mapping equivalent degrees of freedom
  for (gsBoxTopology::const_iiterator it = m_mb.topology().iBegin();
       it != m_mb.topology().iEnd(); ++it) // C^0 at the interface
  {
    m_mb.matchInterface(*it, m_mapper);
  }

  // Mark the boundary in each target dimension for each patch
  for (gsMultiPatch<>::const_biterator bit = m_mp.bBegin();
       bit != m_mp.bEnd(); ++bit)
  {
    gsMatrix<index_t> idx = m_mb.basis(bit->patch).boundary(bit->index());
    for (index_t idim = 0; idim != m_mp.targetDim(); ++idim)
    {
      m_mapper.markBoundary(bit->patch, idx, idim);
    }
  }

  // Finalize the mapper after marking all boundaries
  m_mapper.finalize();

  logMapperInformation();
}

template<short_t d, typename T>
gsDofMapper gsBarrierPatch<d, T>::_makeMapperOnePatch(const gsGeometry<T> &currPatch) const
{
  // Initialize the mapper for the design Degrees of Freedom (DoFs)
  // By default, all the inner control points are set as optimization variables
  // However, it's also possible to set only a part of them as optimization variables
  gsDofMapper mapper;
  gsMultiBasis<T> mb(currPatch);
  mapper.init(mb, currPatch.targetDim());

  // Mark the boundary in each target dimension
  gsMatrix<index_t> idx = currPatch.basis().allBoundary();
  for (index_t idim = 0; idim != m_mp.targetDim(); ++idim)
  {
    mapper.markBoundary(0, idx, idim);
  }

  // Finalize the mapper after marking all boundaries
  mapper.finalize();

  // Return the generated mapper
  return mapper;
}

template<short_t d, typename T>
void gsBarrierPatch<d, T>::logMapperInformation() {
  gsDebug << "#Number of free  variables: " << m_mapper.freeSize() << "\n";
  gsDebug << "#Number of fixed variables: " << m_mapper.boundarySize() << "\n";
  gsDebug << "#Total number of variables: " << m_mapper.size() << "\n\n\n";
}

}// namespace gismo
