/** @file gsNormL2.h

    @brief Computes the L2 norm, needs for the parallel computing.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): P. Weinmüller
*/

#pragma once
#include <gsG1Basis/gsG1System.h>

namespace gismo
{


template <class T>
class gsVisitorNormL2
{

public:

    gsVisitorNormL2(index_t p = 2)
    {
        f2param = false;
        m_p = p;
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
        evFlags = NEED_MEASURE| NEED_VALUE | NEED_JACOBIAN;
    }

    // Evaluate on element.
    void evaluate(gsGeometryEvaluator<T> & geoEval,
                  const gsFunction<T>    & _func2,
                  const std::vector<gsMultiBasis<T>> * basis,
                  const gsSparseMatrix<T> * sol_sparse,
                  gsG1System<T> & g1System,
                  gsMatrix<T>            & quNodes,
                  bool isogeometric)
    {
        gsMatrix<unsigned> actives;
        gsMatrix<T> basisData;

        basis->at(0).basis(geoEval.id()).active_into(quNodes.col(0), actives);

        // Evaluate basis functions on element
        basis->at(0).basis(geoEval.id()).eval_into(quNodes,basisData);

        f1vals.setZero(1,quNodes.cols());
        if (!isogeometric)
        {
            gsMatrix<unsigned> actives2;
            gsMatrix<T> basisData2;

            basis->at(1).basis(geoEval.id()).active_into(quNodes.col(0), actives2);

            // Evaluate basis functions on element
            basis->at(1).basis(geoEval.id()).eval_into(quNodes,basisData2);

            for (index_t i = 0; i < g1System.get_numInterfaceFunctions().last(); i++)
                for (index_t j = 0; j < actives2.rows(); j++)
                    f1vals += sol_sparse->at(i,g1System.get_numBasisFunctionsInterface()[geoEval.id()] + actives2.at(j)) * basisData2.row(j);

            for (index_t i = g1System.get_numBoundaryVertexFunctions()[1]; i < g1System.get_numBoundaryVertexFunctions()[2]; i++)
                for (index_t j = 0; j < actives2.rows(); j++)
                    f1vals += sol_sparse->at(i,g1System.get_numBasisFunctionsInterface()[geoEval.id()] + actives2.at(j)) * basisData2.row(j);

            for (index_t i = g1System.get_numBoundaryVertexFunctions()[3]; i < g1System.get_numBoundaryVertexFunctions()[4]; i++)
                for (index_t j = 0; j < actives2.rows(); j++)
                    f1vals += sol_sparse->at(i,g1System.get_numBasisFunctionsInterface()[geoEval.id()] + actives2.at(j)) * basisData2.row(j);



            for (index_t i = g1System.get_numInterfaceFunctions().last(); i < g1System.get_numBoundaryVertexFunctions()[1]; i++) // -1 bcs of interior solution
                for (index_t j = 0; j < actives.rows(); j++)
                    f1vals += sol_sparse->at(i,g1System.get_numBasisFunctions()[geoEval.id()] + actives.at(j)) * basisData.row(j);

            for (index_t i = g1System.get_numBoundaryVertexFunctions()[2]; i < g1System.get_numBoundaryVertexFunctions()[3]; i++)
                for (index_t j = 0; j < actives.rows(); j++)
                    f1vals += sol_sparse->at(i,g1System.get_numBasisFunctionsInterface()[geoEval.id()] + actives.at(j)) * basisData.row(j);

            for (index_t i = g1System.get_numBoundaryVertexFunctions()[4]; i < g1System.get_numBoundaryVertexFunctions().last(); i++)
                for (index_t j = 0; j < actives.rows(); j++)
                    f1vals += sol_sparse->at(i,g1System.get_numBasisFunctionsInterface()[geoEval.id()] + actives.at(j)) * basisData.row(j);
        }
        else
        {
            for (index_t i = 0; i < g1System.get_numInterfaceFunctions().last(); i++)
                for (index_t j = 0; j < actives.rows(); j++)
                    f1vals += sol_sparse->at(i,g1System.get_numBasisFunctionsInterface()[geoEval.id()] + actives.at(j)) * basisData.row(j);

            for (index_t i = g1System.get_numInterfaceFunctions().last(); i < sol_sparse->rows() - 1; i++) // -1 bcs of interior solution
                for (index_t j = 0; j < actives.rows(); j++)
                    f1vals += sol_sparse->at(i,g1System.get_numBasisFunctions()[geoEval.id()] + actives.at(j)) * basisData.row(j);
        }


        for (index_t j = 0; j < actives.rows(); j++) // interior solution
            f1vals += sol_sparse->at(sol_sparse->rows()-1,g1System.get_numBasisFunctions()[geoEval.id()] + actives.at(j)) * basisData.row(j);


        // Compute geometry related values
        geoEval.evaluateAt(quNodes);

        // Evaluate second function
        _func2.eval_into( f2param ? quNodes : geoEval.values() , f2vals);
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
            gsMatrix<T> G_inv = G.cramerInverse();

            real_t detG = G.determinant();


            //const T weight = quWeights[k] * geoEval.measure(k);
            const T weight = quWeights[k] * sqrt(detG);
            switch (m_p)
            {
                case 0: // infinity norm
                    // .template lpNorm<Eigen::Infinity>();
                    sum = (f1vals - f2vals).array().abs().maxCoeff();
                    accumulated = math::max(accumulated, sum);
                    return sum;
                    break;
                case 1:
                    sum += weight * ( f1vals.col(k) - f2vals.col(k) ).template lpNorm<1>();
                    break;
                case 2:
                    sum += weight * ( f1vals.col(k) - f2vals.col(k) ).squaredNorm();
                    break;
                default:
                    sum += weight * ( f1vals.col(k) - f2vals.col(k) ).array().abs().pow(m_p).sum();
                    //sum += weight * ( f1vals.col(k) - f2vals.col(k) ).template lpNorm<p>().squared();
            }
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







