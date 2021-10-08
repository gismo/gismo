/** @file gsC1SurfGD.h

    @brief Compute the gluing data for one interface.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Farahat
*/

#pragma once

#include <gsUnstructuredSplines/gsC1SurfGluingDataAssembler.h>
#include <gsIO/gsOptionList.h>


namespace gismo {

    template<class T>
    class gsC1SurfGD {

    public:

        gsC1SurfGD() {}

        gsC1SurfGD(gsMultiPatch <T> const &mp, gsMultiBasis <T> &mb)
                : m_mp(mp), m_mb(mb) {
        }


        gsC1SurfGD(gsMultiPatch <T> const &mp, gsMultiBasis <T> const &mb, index_t uv, bool isBoundary,
                     gsG1OptionList const &optionList)
                : m_mp(mp), m_mb(mb), m_uv(uv), m_isBoundary(isBoundary), m_optionList(optionList) {
            m_gamma = 1.0;
        }

        const gsBSpline <T> get_alpha_tilde() const { return alpha_tilde; }

        const gsBSpline <T> get_beta_tilde() const { return beta_tilde; }

    protected:

        gsBSpline <T> alpha_tilde;
        gsBSpline <T> beta_tilde;

        // The geometry for a single interface in the right parametrizations
        gsMultiPatch <T> m_mp;
        gsMultiBasis <T> m_mb;
        index_t m_uv;
        bool m_isBoundary;
        gsG1OptionList m_optionList;

        real_t m_gamma;


    };

}
