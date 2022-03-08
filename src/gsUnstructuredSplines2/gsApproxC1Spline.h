/** @file gsApproxC1Spline.h

    @brief Construct the approx c1 spline space.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): P. Weinmueller
*/

#pragma once

#include<gsIO/gsOptionList.h>

#include <gsMSplines/gsMappedBasis.h>
#include<gsUnstructuredSplines2/gsContainerBasisBase.h>

//#include <gsUnstructuredSplines/gsApproxC1Utils.h>

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
        this->defaultOptions();
    };

public:
    // To be overwritten in inheriting classes

    void init();
    void compute();

    void update(gsMappedBasis<2,real_t> & bb2, bool compute = true)
    {
        if (compute) {
            this->init();
            this->compute();
        }

        m_matrix = m_matrix.transpose();
        gsMultiBasis<> dbasis_temp;
        this->getMultiBasis(dbasis_temp);

        bb2.init(dbasis_temp,m_matrix);
    }

    void defaultOptions();

private:
    // Helper functions


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
