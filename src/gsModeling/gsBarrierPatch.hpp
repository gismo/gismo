/** @file gsBarrierPatch.hpp

@brief Provides patch construction from boundary data by using barrier method. It is a
reference implementation of the following paper. If you make use of the code or the
idea/algorithm in your work, please cite our paper:
Ji, Y., Yu, Y. Y., Wang, M. Y., & Zhu, C. G. (2021).
Constructing high-quality planar NURBS parameterization for
isogeometric analysis by adjustment control points and weights.
Journal of Computational and Applied Mathematics, 396, 113615.
(https://www.sciencedirect.com/science/article/pii/S0377042721002375)

This file is part of the G+Smo library.

This Source Code Form is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this
file, You can obtain one at http://mozilla.org/MPL/2.0/.

Author(s): Ye Ji, H.M. Verhelst
*/

#pragma once

#include <gsModeling/gsBarrierCore.h>
#include <gsModeling/gsBarrierCore.hpp>

namespace gismo
{
/**
\brief Computes a patch parametrization given a set of
boundary geometries.  Parametrization is not guaranteed to be
non-singular. Works for surface, volumes.

\tparam d dimension of physical domain
\tparam T Coefficient type

\param bRep Input boundary representation
\param initialMethod Specify initialization method
\param filename Name of input data file

\ingroup Modeling
*/

// default constructor
    template<short_t d, typename T>
    gsBarrierPatch<d, T>::gsBarrierPatch(const gsMultiPatch<T> &mp,
                                         const gsDofMapper &mapper)
            :
            m_mp(mp),
            m_mb(mp)
    {
        // Input data Sanity Check
        this->defaultOptions();
        this->setMapper(mapper);
    }

// // multipatch constructor
// // patch-wise constructor
    template<short_t d, typename T>
    gsBarrierPatch<d, T>::gsBarrierPatch(const gsMultiPatch<T> &mp,
                                         bool patchWise)
            :
            m_mp(mp),
            m_mb(mp)
    {
        this->defaultOptions();
        if (patchWise)
//            _makeMapperLocalPatches();
            m_freeInterface = 0;
        else
        {
            _makeMapperGlobalPatches();
        }

        // if patchWise==true
        // make mapper that also fixes interfaces --> current _makeMapper
        // if false
        // fix boundaries and match interfaces
        // _makeMapper();
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
        // preprocessing: scale the computational domain to [0,10]^d for better
        // numerical stability and to have a consistent convergence criteria
        const double boxsize = 1.0;
        gsMatrix<T> boundingBox;
        m_mp.boundingBox(boundingBox);
        gsVector<T, d> bbmin = boundingBox.col(0);
        gsVector<T, d> bbmax = boundingBox.col(1);
        real_t maxside = (bbmax-bbmin).maxCoeff();

        gsVector<T, d> scaleFactor;
        scaleFactor.setConstant(boxsize/maxside);
        for (auto &ptch:m_mp) {
            ptch->translate(-bbmin);
            ptch->scale(scaleFactor);
        }

        // start parameterization construction
        gsStopwatch timer;
        switch (m_freeInterface)
        {
            case 0:
            {
                // construct analysis-suitable parameterization piecewisely
                for (auto iptch = 0; iptch < m_mp.nPatches(); ++iptch)
                {
                    gsInfo << "I am parameterizing " << iptch
                           << "-th patch, total number of patches is "
                           << m_mp.nPatches() << "\n";
                    gsDofMapper currMapper = _makeMapperOnePatch(
                            m_mp.patch(iptch));
                    gsMultiPatch<T> optCurrPatch = gsBarrierCore<d, T>::compute(
                            m_mp.patch(iptch), currMapper, m_options);
                    m_mp.patch(iptch).setCoefs(optCurrPatch.patch(0).coefs());
                }
                break;
            }
            case 1:
            {
                // construct analysis-suitable parameterization with moving interfaces
                m_mp = gsBarrierCore<d, T>::compute(m_mp, m_mapper, m_options);
                break;
            }
        }

        real_t runtime = timer.stop();
        gsInfo  << "\n"
                << "Multi-patch parameterization construction completed! Running time is "
                << runtime << "\n";

        // restore scale to the original size of the input model
        for (auto &isf:scaleFactor)
            isf = 1./isf;
        for (auto &ptch:m_mp) {
            ptch->scale(scaleFactor);
            ptch->translate(bbmin);
        }
    }

    template<short_t d, typename T>
    void gsBarrierPatch<d, T>::defaultOptions()
    {
        m_options = gsBarrierCore<d, T>::defaultOptions();
    }

    template<short_t d, typename T>
    void gsBarrierPatch<d, T>::_makeMapper()
    {
        ////! [Make mapper for the design DoFs]
        // Now, we set all the inner control points as optimization variables
        // It is also possible to set only a part of them as optimization variables
        m_mapper.init(m_mb, m_mp.targetDim());
        for (size_t iptch = 0; iptch != m_mp.nPatches(); iptch++)
        {
            gsMatrix<index_t> idx = m_mp.basis(iptch).allBoundary();
            for (index_t idim = 0; idim != m_mp.targetDim(); idim++)
            {
                m_mapper.markBoundary(iptch, idx, idim);
            }
        }
        m_mapper.finalize();

        gsInfo << "#Numb of free  variables is " << m_mapper.freeSize() << "\n";
        gsInfo << "#Numb of fixed variables is " << m_mapper.boundarySize()
               << "\n";
        gsInfo << "#Numb of total variables is " << m_mapper.size() << "\n\n\n";
    }

    template<short_t d, typename T>
    void gsBarrierPatch<d, T>::_makeMapperLocalPatches()
    {
        ////! [Make mapper for the design DoFs]
        // Now, we set all the inner control points as optimization variables
        // It is also possible to set only a part of them as optimization variables
        m_mapper.init(m_mb, m_mp.targetDim());
        for (gsMultiPatch<>::const_biterator bit = m_mp.bBegin();
             bit != m_mp.bEnd(); ++bit)
        {
            gsMatrix<index_t> idx = m_mb.basis(bit->patch).allBoundary();
            for (index_t idim = 0; idim != m_mp.targetDim(); ++idim)
                m_mapper.markBoundary(bit->patch, idx, idim);
        }
        m_mapper.finalize();

        gsInfo << "#Numb of free  variables is " << m_mapper.freeSize() << "\n";
        gsInfo << "#Numb of fixed variables is " << m_mapper.boundarySize()
               << "\n";
        gsInfo << "#Numb of total variables is " << m_mapper.size() << "\n\n\n";
    }

    template<short_t d, typename T>
    void gsBarrierPatch<d, T>::_makeMapperGlobalPatches()
    {
        ////! [Make mapper for the design DoFs]
        // Now, we set all the inner control points as optimization variables
        // It is also possible to set only a part of them as optimization variables
        m_mapper.init(m_mb, m_mp.targetDim());

        for (gsBoxTopology::const_iiterator it = m_mb.topology().iBegin();
             it != m_mb.topology().iEnd(); ++it) // C^0 at the interface
            m_mb.matchInterface(*it, m_mapper);

        for (gsMultiPatch<>::const_biterator bit = m_mp.bBegin();
             bit != m_mp.bEnd(); ++bit)
        {
            gsMatrix<index_t> idx = m_mb.basis(bit->patch).boundary(
                    bit->index());
            for (index_t idim = 0; idim != m_mp.targetDim(); ++idim)
                m_mapper.markBoundary(bit->patch, idx, idim);
        }

        m_mapper.finalize();
        gsInfo << "#Numb of free  variables is " << m_mapper.freeSize() << "\n";
        gsInfo << "#Numb of fixed variables is " << m_mapper.boundarySize()
               << "\n";
        gsInfo << "#Numb of total variables is " << m_mapper.size() << "\n\n\n";
    }

    template<short_t d, typename T>
    gsDofMapper
    gsBarrierPatch<d, T>::_makeMapperOnePatch(const gsGeometry<T> &currPatch)
    {
        ////! [Make mapper for the design DoFs]
        // Now, we set all the inner control points as optimization variables
        // It is also possible to set only a part of them as optimization variables
        gsDofMapper mapper;
        gsMultiBasis<T> mb(currPatch);
        mapper.init(mb, currPatch.targetDim());
        gsMatrix<index_t> idx = currPatch.basis().allBoundary();
        for (index_t idim = 0; idim != m_mp.targetDim(); ++idim)
            mapper.markBoundary(0, idx, idim);
        mapper.finalize();
        return mapper;
    }


    template<short_t d, typename T>
    gsMultiPatch<T> &gsBarrierPatch<d, T>::result() { return m_mp; }

}// namespace gismo
