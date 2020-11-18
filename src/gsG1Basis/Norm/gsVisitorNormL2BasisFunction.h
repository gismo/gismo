/** @file gsNormL2.h

    @brief Computes the Semi H1 norm, needs for the parallel computing.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): P. Weinmüller
*/

#include <gsCore/gsGeometryEvaluator.h>

#pragma once

namespace gismo
{


template <class T>
class gsVisitorNormL2BasisFunction
{
public:

    gsVisitorNormL2BasisFunction()
    {
        f2param = false;

    }

    void initialize(const gsBasis<T> & basis,
                    gsQuadRule<T> & rule,
                    unsigned      & evFlags) // replace with geoEval ?
    {
        // Setup Quadrature
        const unsigned d = basis.dim();
        gsVector<index_t> numQuadNodes( d );
        for (unsigned i = 0; i < d; ++i)
            numQuadNodes[i] = basis.degree(i) + 1;

        // Setup Quadrature
        rule = gsGaussRule<T>(numQuadNodes);// harmless slicing occurs here

        // Set Geometry evaluation flags
        evFlags = NEED_MEASURE|NEED_VALUE|NEED_GRAD_TRANSFORM;
    }

    // Evaluate on element.
    void evaluate(gsGeometryEvaluator<T> & geoEval,
                  gsMatrix<T>  & quNodes,
                  const gsMultiPatch<T> * patchesPtr,
                  gsMultiBasis<T> * basis,
                  const gsBSplineBasis<> & basis_plus,
                  const gsBSplineBasis<> & basis_minus,
                  const gsApproxGluingData<T> & gD,
                  const index_t uv,
                  index_t pn)
    {
        gsMatrix<unsigned> actives;
        gsMatrix<T> bGrads;

        f1ders.setZero(1,quNodes.cols());
        f1ders = patchesPtr->patch(pn).eval(quNodes);

        rhsGrads.setZero(2,quNodes.cols());
        rhsGrads = patchesPtr->patch(pn).deriv(quNodes);

        gsBSplineBasis<T> bsp_temp = dynamic_cast<gsBSplineBasis<> & >(basis->basis(0).component(1-uv));

        real_t p = basis->basis(0).component(1-uv).maxDegree();
        real_t tau_1 = bsp_temp.knots().at(p + 1); // p + 2

        gsMatrix<T> alpha, beta,
            N_0, N_1,
            N_j_minus, N_i_plus,
            der_N_i_plus, result_test;

        gsMatrix<T> der_alpha, der_beta,
            der_N_0, der_N_1,
            der_N_j_minus,
            der2_N_i_plus;

        // pn+2 falls non-two Patch
        basis_plus.evalSingle_into(pn,quNodes.row(uv),N_i_plus); // v
        basis_plus.derivSingle_into(pn,quNodes.row(uv),der_N_i_plus);

        basis_minus.evalSingle_into(pn,quNodes.row(uv),N_j_minus);
        basis_minus.derivSingle_into(pn,quNodes.row(uv),der_N_j_minus);

        basis->basis(0).component(1-uv).evalSingle_into(0,quNodes.row(1-uv),N_0); // v
        basis->basis(0).component(1-uv).evalSingle_into(1,quNodes.row(1-uv),N_1); // v

        basis->basis(0).component(1-uv).derivSingle_into(1,quNodes.row(1-uv),der_N_1); // v

        gD.get_beta_S_tilde(1-uv).eval_into(quNodes.row(uv),beta);
        gD.get_alpha_S_tilde(1-uv).eval_into(quNodes.row(uv),alpha);

        gD.get_alpha_S_tilde(1-uv).deriv_into(quNodes.row(uv),der_alpha);

        gsMatrix<> temp = beta.cwiseProduct(der_N_i_plus);
        //f2ders = N_i_plus.cwiseProduct(N_0 + N_1) - temp.cwiseProduct(N_1) * tau_1 / p; // beta
        f2ders = ( uv == 0 ? -1 : +1 ) * alpha.cwiseProduct(N_j_minus.cwiseProduct(N_1)) * tau_1 / p;
        //f2ders = ( uv == 0 ? -1 : +1 ) * N_j_minus.cwiseProduct(N_1) * tau_1 / p;


        rhsGrads2.setZero(2,quNodes.cols());
        rhsGrads2.row(1-uv) = ( uv == 0 ? -1 : +1 ) * alpha.cwiseProduct(N_j_minus.cwiseProduct(der_N_1)) * tau_1 / p;
        rhsGrads2.row(uv) = ( uv == 0 ? -1 : +1 ) * (der_alpha.cwiseProduct(N_j_minus)+alpha.cwiseProduct(der_N_j_minus)).cwiseProduct(N_1) * tau_1 / p;
    }


    // assemble on element
    inline T compute(gsDomainIterator<T>    & geo,
                     gsGeometryEvaluator<T> & geoEval,
                     gsVector<T> const      & quWeights,
                     T & accumulated)
    {
        T sum(0.0);
        for (index_t k = 0; k < quWeights.rows(); ++k) // loop over quadrature nodes
        {
            const T weight = quWeights[k];
            sum += weight * (f1ders.col(k) - f2ders.col(k)).squaredNorm();
        }
        accumulated += sum;
        return sum;
    }


protected:

    gsMatrix<T> f1ders, f2ders, rhsGrads, rhsGrads2;

    bool f2param;

};






}







