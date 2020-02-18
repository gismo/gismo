/** @file gsGluingData.h

    @brief Compute the gluing data for one interface.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): P. Weinmueller
*/

#pragma once

//#include <gsG1Basis/gsGlobalGDAssembler.h>


namespace gismo
{

template<class T>
class gsGluingData
{
public:
    gsGluingData()
    { }

    gsGluingData(gsMultiPatch<T> const & mp,
                 gsMultiBasis<T> const & mb,
                 gsOptionList const & optionList)
        : m_mp(mp), m_mb(mb), m_optionList(optionList)
    {
        m_gamma = 1.0;

        p_tilde = m_optionList.getInt("p_tilde");
        r_tilde = m_optionList.getInt("r_tilde");

        if (m_optionList.getSwitch("local"))
            gsInfo << "Is not yet implemented \n";
        else
            setGlobalGluingData();


        gsInfo << "MP: " << m_mp.patch(0).coefs() << "\n";
        gsInfo << "MP: " << m_mp.patch(1).coefs() << "\n";
    }

    // Not exact!!! Interpolation via Greville abscissae // But only necessary for condition test!
    void beta_exact();

    // Computed the gluing data globally
    void setGlobalGluingData();

    void conditionTest()
    {
        gsMatrix<> points(1,1000);
        points.setRandom();
        points = points.array().abs();

        beta_exact();

        gsMatrix<> temp;

        temp = alpha_tilde_L.eval(points).cwiseProduct(beta_tilde_R.eval(points))
            - alpha_tilde_R.eval(points).cwiseProduct(beta_tilde_L.eval(points))
            + beta_bar.eval(points);

        gsInfo << "Conditiontest gluing data : " << temp.array().abs().maxCoeff() << "\n\n";

    }
protected:
    // The geometry for a single interface in the right parametrizations
    gsMultiPatch<T> m_mp;
    gsMultiBasis<T> m_mb;
    gsOptionList m_optionList;

    real_t m_gamma;

    // Spline space for the gluing data (p_tilde,r_tilde,k)
    index_t p_tilde, r_tilde;

protected:
    // Global Gluing data
    gsBSpline<T> beta_bar;
    gsBSpline<T> alpha_tilde_L, alpha_tilde_R;
    gsBSpline<T> beta_tilde_L, beta_tilde_R;

}; // class gsGluingData

template <class T>
void gsGluingData<T>::beta_exact()
{

} // beta_exact

template<class T>
void gsGluingData<T>::setGlobalGluingData()
{

} // setGlobalGluingData

} // namespace gismo

