/** @file gsNormL2.h

    @brief Computes the Semi H2 norm, needs for the parallel computing.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Farahat
*/


#pragma once

namespace gismo
{


template <class T>
class gsG1ASVisitorResidualSeminormH2
{
public:

    gsG1ASVisitorResidualSeminormH2()
    {
    }

    void initialize(const std::vector<gsMultiBasis<>> * basis,
                    gsQuadRule<T> & rule,
                    unsigned      & evFlags) // replace with geoEval ?
    {
        // Setup Quadrature
        const unsigned d = basis->at(0).dim();
        gsVector<index_t> numQuadNodes( d );
        for (unsigned i = 0; i < d; ++i)
            numQuadNodes[i] = basis->at(0).degree(i) + 1;

        // Setup Quadrature
        rule = gsGaussRule<T>(numQuadNodes);// harmless slicing occurs here

        // Set Geometry evaluation flags
        evFlags = NEED_MEASURE|NEED_VALUE|NEED_GRAD_TRANSFORM|NEED_2ND_DER|NEED_JACOBIAN;
    }

    // Evaluate on element.
    void evaluate(gsGeometryEvaluator<T> & geoEval,
                  const std::vector<gsSparseMatrix<>> * sol_sparse,
                  const std::vector<gsMultiBasis<>> * basis_vec,
                  std::vector<gsG1System<real_t>> & sys_vec,
                  gsMatrix<T>            & quNodes)
    {
        gsMatrix<unsigned> actives;
        gsMatrix<T> derivData, deriv2Data;

        qN = quNodes;

        basis_vec->at(0).basis(geoEval.id()).active_into(quNodes.col(0), actives);

        // Evaluate basis functions on element
        basis_vec->at(0).basis(geoEval.id()).deriv_into(quNodes,derivData);
        basis_vec->at(0).basis(geoEval.id()).deriv2_into(quNodes,deriv2Data);
        geoEval.evaluateAt(quNodes);

        f1ders.setZero(2,actives.rows());
        f1ders2.setZero(3,actives.rows());
        for (index_t i = 0; i < sol_sparse->at(0).rows(); i++)
            for (index_t j = 0; j < actives.rows(); j++)
            {
                f1ders += sol_sparse->at(0).at(i,sys_vec.at(0).get_numBasisFunctions()[geoEval.id()] + actives.at(j)) * derivData.block(2*j,0,2,f1ders.dim().second);
                f1ders2 += sol_sparse->at(0).at(i,sys_vec.at(0).get_numBasisFunctions()[geoEval.id()] + actives.at(j)) * deriv2Data.block(3*j,0,3,f1ders.dim().second);
            }

        basis_vec->at(1).basis(geoEval.id()).active_into(quNodes.col(0), actives);

        // Evaluate basis functions on element
        basis_vec->at(1).basis(geoEval.id()).deriv_into(quNodes,derivData);
        basis_vec->at(1).basis(geoEval.id()).deriv2_into(quNodes,deriv2Data);
        geoEval.evaluateAt(quNodes);

        f2ders.setZero(2,actives.rows());
        f2ders2.setZero(3,actives.rows());
        for (index_t i = 0; i < sol_sparse->at(1).rows(); i++)
            for (index_t j = 0; j < actives.rows(); j++)
            {
                f2ders += sol_sparse->at(1).at(i,sys_vec.at(1).get_numBasisFunctions()[geoEval.id()] + actives.at(j)) * derivData.block(2*j,0,2,f1ders.dim().second);
                f2ders2 += sol_sparse->at(1).at(i,sys_vec.at(1).get_numBasisFunctions()[geoEval.id()] + actives.at(j)) * deriv2Data.block(3*j,0,3,f1ders.dim().second);
            }


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

            gsMatrix<T> Jk = geoEval.jacobian(k);
            gsMatrix<T> G = Jk.transpose() * Jk;

            T weight = quWeights[k] * sqrt(G.determinant());


                sum += weight * ( (f1ders2 - f2ders2).squaredNorm());

        }
        accumulated += sum;
        return sum;
    }

private:


    gsMatrix<T> f1ders, f1ders2, f2ders, f2ders2;

    std::vector<gsMatrix<T> > basisData;
    gsMatrix<T> basisGrads;
    gsMatrix<T> basis2ndDerivs;
    index_t numActive;
    gsMatrix<T> qN;
};

}






