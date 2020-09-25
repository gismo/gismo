/** @file gsNormL2.h

    @brief Computes the L2 norm, needs for the parallel computing.

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
class gsG1ASVisitorResidualNormL2
{

public:

    gsG1ASVisitorResidualNormL2(index_t p = 2)
    {
        f2param = false;
        m_p = p;
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
        evFlags = NEED_MEASURE| NEED_VALUE | NEED_JACOBIAN;
    }

    // Evaluate on element.
    void evaluate(gsGeometryEvaluator<T> & geoEval,
                  const std::vector<gsSparseMatrix<>> * sol_sparse,
                  const std::vector<gsMultiBasis<>> * basis_vec,
                  std::vector<gsG1System<real_t>> & sys_vec,
    gsMatrix<T>            & quNodes)
    {
        gsMatrix<unsigned> actives;
        gsMatrix<T> basisData;

        basis_vec->at(0).basis(geoEval.id()).active_into(quNodes.col(0), actives);

        // Evaluate basis functions on element
        basis_vec->at(0).basis(geoEval.id()).eval_into(quNodes,basisData);

        f1vals.setZero(1,actives.rows());
        for (index_t i = 0; i < sol_sparse->at(0).rows(); i++) // -1 bcs of interior solution
            for (index_t j = 0; j < actives.rows(); j++)
                f1vals += sol_sparse->at(0).at(i,sys_vec.at(0).get_numBasisFunctions()[geoEval.id()] + actives.at(j)) * basisData.row(j);

        // Compute geometry related values
        geoEval.evaluateAt(quNodes);

        basis_vec->at(1).basis(geoEval.id()).active_into(quNodes.col(0), actives);

        // Evaluate basis functions on element
        basis_vec->at(1).basis(geoEval.id()).eval_into(quNodes,basisData);

        f2vals.setZero(1,actives.rows());
        for (index_t i = 0; i < sol_sparse->at(1).rows(); i++) // -1 bcs of interior solution
            for (index_t j = 0; j < actives.rows(); j++)
                f2vals += sol_sparse->at(1).at(i,sys_vec.at(1).get_numBasisFunctions()[geoEval.id()] + actives.at(j)) * basisData.row(j);

        // Compute geometry related values
        geoEval.evaluateAt(quNodes);
    }

    // assemble on element
    T compute(gsDomainIterator<T>    & ,
              gsGeometryEvaluator<T> & geoEval,
              gsVector<T> const      & quWeights,
              T & accumulated )
    {
        T sum(0.0);
        for (index_t k = 0; k < quWeights.rows(); ++k) // loop over quadrature nodes
        {


            gsMatrix<T> Jk = geoEval.jacobian(k);
            gsMatrix<T> G = Jk.transpose() * Jk;
            real_t detG = G.determinant();
            T weight = quWeights[k] * sqrt(detG);

            sum += weight * ( f1vals.col(k) - f2vals.col(k) ).squaredNorm();
        }
        accumulated += sum;
        return sum;
    }

private:

    index_t m_p;

    gsMatrix<T> f1vals, f2vals;

    bool f2param;

};






}







