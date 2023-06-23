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

#include <gismo.h>
#include <gsHLBFGS/gsHLBFGS.h>

using namespace gismo;

namespace gismo
{

    /**
      \brief Computes a patch parametrization given a set of
      boundary geometries.  Parametrization is not guaranteed to be
      non-singular. Works for surface, volumes, or any dimension. (not
      yet, only works for planar domain now)

      \tparam d domain dimension
      \tparam T Coefficient type

      \ingroup Modeling
      */

    // It is a barrier function-based method. The main step ff
    template<short_t d, typename T>
    class gsBarrierPatch
    {
    public:

        gsBarrierPatch(const gsMultiPatch<T> &mp, const gsDofMapper &mapper);

        gsBarrierPatch(const gsMultiPatch<T> &mp, bool patchWise = false);
        // gsBarrierPatch(const gsMultiPatch<T> &mp, patchCorner corner);
        // gsBarrierPatch(const gsMultiPatch<T> &mp, boundaryInterface iface);

    public:
        void setMapper(const gsDofMapper &mapper) { m_mapper = mapper; }

        gsMultiPatch<T> &result();

        /// compute analysis-suitable by different methods
        void compute();

        gsOptionList &options() { return m_options; }

        void defaultOptions();

    protected:
        void _makeMapper();

        gsDofMapper _makeMapperOnePatch(const gsGeometry<T>&currPatch);

        void _makeMapperGlobalPatches();

        void _makeMapperLocalPatches();

    private:

        mutable gsExprEvaluator<T> m_evaluator;
        mutable gsExprAssembler<T> m_assembler;
        gsDofMapper m_mapper;

        mutable gsMultiPatch<T> m_mp;
        gsMultiBasis<T> m_mb;
        gsMultiPatch<T> m_bRep;
//        T m_area; // area of computational domain
//        T m_scaledArea = 1.0;

        std::string m_filename;

        gsVector<T, d> m_boundingBoxLeftBottomCorner;
        gsVector<T, d> m_scalingVec;

//        const T m_boxsize = 1.0;
//        real_t m_scalingFactor = 1.0;

        size_t m_freeInterface = 1;
        gsOptionList m_options;
    };
}// namespace gismo

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsBarrierPatch.hpp)
#endif
