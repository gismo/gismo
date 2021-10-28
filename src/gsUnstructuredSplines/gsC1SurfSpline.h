/** @file gsC1SurfSpline.h

    @brief Provides declaration of Basis abstract interface.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): F. Buchegger
*/

/*
    TO DO
 */

#pragma once

#include<gsIO/gsOptionList.h>
#include<gsUnstructuredSplines/gsContainerBasisBase.h>

namespace gismo
{

    template<short_t d,class T>
    class gsC1SurfSpline : public gsContainerBasisBase<d,T>
    {
    public:

        // gsContainerBasisBase:
        // - Interior space: [0] : inner,
        // - Edge spaces:    [1] : west, [2] : east, [3] : south, [4] : north,
        // - Vertex spaces:  [5] : southwest, [6] : southeast, [7] : northwest, [8] : northeast
        using Base = gsContainerBasisBase<d,T>;

        gsC1SurfSpline(gsMultiPatch<T> & patches, gsMultiBasis<T> & multiBasis)
                :
                Base(patches, multiBasis)
        {
            this->defaultOptions();
        };

    public:
        // To be overwritten in inheriting classes

        void init();
        void compute();

    private:
        void defaultOptions();

    protected:
        // Container[patch][side]
        std::vector<std::vector<index_t>> rowContainer;

        // Put here the members of the shared functions
        using Base::m_patches;
        using Base::m_multiBasis;
        using Base::m_options;
        using Base::m_matrix;
        using Base::m_bases;
    };

}

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsC1SurfSpline.hpp)
#endif
