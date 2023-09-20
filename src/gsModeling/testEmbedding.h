//
// Created by Jingya Li on 04/09/2023.
//

#pragma once

#include <gsCore/gsGeometry.h>
#include <gsCore/gsMultiPatch.h>
#include <gsNurbs/gsBSpline.h>
#include <gsNurbs/gsNurbs.h>

namespace gismo {

/** \brief Class that performs an L2 projection

    \tparam T coefficient type
 */

    template <class T>
    struct testEmbedding
    {

//    protected:
//        typedef gsExprAssembler<>::geometryMap geometryMap;
//        typedef gsExprAssembler<>::space       space;
//        typedef gsExprAssembler<>::solution    solution;
//        typedef gsExprAssembler<>::element     element;

    public:
        /**
         * @brief      Projects a \a source geometry onto \a basis and returns it in
         *             \a result
         *
         * @param[in]  basis     The basis to project on
         * @param[in]  geometry  The geometry
         * @param      result    The new geometry
         */
         static T EmbedCurvesOnSurface(const gsGeometry<T> & curve,
                                      const gsGeometry<T> & surface,
                                       gsGeometry<T> & result,
                                      int numEval);


    };

} // gismo

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(testEmbedding.hpp)
#endif
