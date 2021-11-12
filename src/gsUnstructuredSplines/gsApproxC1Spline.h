/** @file gsApproxC1Spline.h

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
class gsApproxC1Spline : public gsContainerBasisBase<d,T>
{
public:

    // gsContainerBasisBase:
    // - Interior space: [0] : inner,
    // - Edge spaces:    [1] : west, [2] : east, [3] : south, [4] : north,
    // - Vertex spaces:  [5] : southwest, [6] : southeast, [7] : northwest, [8] : northeast
    using Base = gsContainerBasisBase<d,T>;

    gsApproxC1Spline(gsMultiPatch<T> & patches, gsMultiBasis<T> & multiBasis)
    :
    Base(patches, multiBasis)
    {
        defaultOptions();
    };

public:
    // To be overwritten in inheriting classes

    void init();
    void compute();

    void defaultOptions();

private:

    // Helper functions
    void createPlusMinusSpace(gsKnotVector<T> & kv1, gsKnotVector<T> & kv2,
                              gsKnotVector<T> & kv1_patch, gsKnotVector<T> & kv2_patch,
                              gsKnotVector<T> & kv1_result, gsKnotVector<T> & kv2_result);

    void createPlusMinusSpace(gsKnotVector<T> & kv1,
                              gsKnotVector<T> & kv1_patch,
                              gsKnotVector<T> & kv1_result, gsKnotVector<T> & kv2_result);

    void createGluingDataSpace(gsKnotVector<T> & kv1, gsKnotVector<T> & kv2,
                               gsKnotVector<T> & kv1_patch, gsKnotVector<T> & kv2_patch,
                               gsKnotVector<T> & kv_result);

    void createLocalEdgeSpace(gsKnotVector<T> & kv_plus, gsKnotVector<T> & kv_minus,
                              gsKnotVector<T> & kv_gD_1, gsKnotVector<T> & kv_gD_2,
                              gsKnotVector<T> & kv_1, gsKnotVector<T> & kv_2,
                              gsKnotVector<T> & kv1_result, gsKnotVector<T> & kv2_result);

    void createLocalEdgeSpace(gsKnotVector<T> & kv_plus, gsKnotVector<T> & kv_minus,
                              gsKnotVector<T> & kv_1,
                              gsKnotVector<T> & kv1_result);

protected:
    // Data members
    index_t p_tilde, r_tilde;

    // Put here the members of the shared functions
    using Base::m_patches;
    using Base::m_multiBasis;
    using Base::m_options;
    using Base::m_matrix;
    using Base::m_bases;
};

}

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsApproxC1Spline.hpp)
#endif
