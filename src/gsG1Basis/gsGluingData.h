//
// Created by afarahat on 3/26/20.
//

#pragma once

#include <gsG1Basis/gsG1ASGluingDataAssembler.h>
# include <gsG1Basis/gsG1OptionList.h>


namespace gismo
{

template<class T>
class gsGluingData
{

public:

    gsGluingData()
    { }

    gsGluingData(gsMultiPatch<T> const & mp, gsMultiBasis<T> const & mb, index_t uv, bool isBoundary, gsG1OptionList const & optionList)
        : m_mp(mp), m_mb(mb), m_uv(uv), m_isBoundary(isBoundary), m_optionList(optionList)
    {
        m_gamma = 1.0;
    }

    const gsBSpline<T> get_alpha_tilde() const {return alpha_tilde; }
    const gsBSpline<T> get_beta_tilde() const {return beta_tilde; }

protected:

    gsBSpline<T> alpha_tilde;
    gsBSpline<T> beta_tilde;

    // The geometry for a single interface in the right parametrizations
    gsMultiPatch<T> m_mp;
    gsMultiBasis<T> m_mb;
    index_t m_uv;
    bool m_isBoundary;
    gsG1OptionList m_optionList;

    real_t m_gamma;


};

} // namespace gismo
