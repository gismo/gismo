/** @file gsNormL2.h

    @brief Computes the Semi H2 norm, needs for the parallel computing.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): P. Weinmüller
*/


#pragma once

namespace gismo
{


template <class T>
class gsVisitorSeminormH2
{
public:

    gsVisitorSeminormH2()
    {

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
        evFlags = NEED_MEASURE|NEED_VALUE|NEED_GRAD_TRANSFORM|NEED_2ND_DER;
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
        gsMatrix<T> derivData, deriv2Data;

        basis->at(0).basis(geoEval.id()).active_into(quNodes.col(0), actives);

        // Evaluate basis functions on element
        basis->at(0).basis(geoEval.id()).deriv_into(quNodes,derivData);
        basis->at(0).basis(geoEval.id()).deriv2_into(quNodes,deriv2Data);

        f1ders.setZero(2,quNodes.cols());
        f1ders2.setZero(3,quNodes.cols());
        if (!isogeometric)
        {
            gsMatrix<unsigned> actives2;
            gsMatrix<T> derivData2, deriv2Data2;

            basis->at(1).basis(geoEval.id()).active_into(quNodes.col(0), actives2);

            // Evaluate basis functions on element
            basis->at(1).basis(geoEval.id()).deriv_into(quNodes,derivData2);
            basis->at(1).basis(geoEval.id()).deriv2_into(quNodes,deriv2Data2);


            for (index_t i = 0; i < g1System.get_numInterfaceFunctions().last(); i++)
                for (index_t j = 0; j < actives2.rows(); j++)
                {
                    f1ders += sol_sparse->at(i,g1System.get_numBasisFunctionsInterface()[geoEval.id()] + actives2.at(j)) * derivData2.block(2*j,0,2,f1ders.dim().second);
                    f1ders2 += sol_sparse->at(i,g1System.get_numBasisFunctionsInterface()[geoEval.id()] + actives2.at(j)) * deriv2Data2.block(3*j,0,3,f1ders.dim().second);
                }

            for (index_t i = g1System.get_numBoundaryVertexFunctions()[1]; i < g1System.get_numBoundaryVertexFunctions()[2]; i++)
                for (index_t j = 0; j < actives2.rows(); j++)
                {
                    f1ders += sol_sparse->at(i,g1System.get_numBasisFunctionsInterface()[geoEval.id()] + actives2.at(j)) * derivData2.block(2*j,0,2,f1ders.dim().second);
                    f1ders2 += sol_sparse->at(i,g1System.get_numBasisFunctionsInterface()[geoEval.id()] + actives2.at(j)) * deriv2Data2.block(3*j,0,3,f1ders.dim().second);
                }

            for (index_t i = g1System.get_numBoundaryVertexFunctions()[3]; i < g1System.get_numBoundaryVertexFunctions()[4]; i++)
                for (index_t j = 0; j < actives2.rows(); j++)
                {
                    f1ders += sol_sparse->at(i,g1System.get_numBasisFunctionsInterface()[geoEval.id()] + actives2.at(j)) * derivData2.block(2*j,0,2,f1ders.dim().second);
                    f1ders2 += sol_sparse->at(i,g1System.get_numBasisFunctionsInterface()[geoEval.id()] + actives2.at(j)) * deriv2Data2.block(3*j,0,3,f1ders.dim().second);
                }



            for (index_t i = g1System.get_numInterfaceFunctions().last(); i < g1System.get_numBoundaryVertexFunctions()[1]; i++)
                for (index_t j = 0; j < actives.rows(); j++)
                {
                    f1ders += sol_sparse->at(i,g1System.get_numBasisFunctions()[geoEval.id()] + actives.at(j)) * derivData.block(2*j,0,2,f1ders.dim().second);
                    f1ders2 += sol_sparse->at(i,g1System.get_numBasisFunctions()[geoEval.id()] + actives.at(j)) * deriv2Data.block(3*j,0,3,f1ders.dim().second);
                }

            for (index_t i = g1System.get_numBoundaryVertexFunctions()[2]; i < g1System.get_numBoundaryVertexFunctions()[3]; i++)
                for (index_t j = 0; j < actives.rows(); j++)
                {
                    f1ders += sol_sparse->at(i,g1System.get_numBasisFunctions()[geoEval.id()] + actives.at(j)) * derivData.block(2*j,0,2,f1ders.dim().second);
                    f1ders2 += sol_sparse->at(i,g1System.get_numBasisFunctions()[geoEval.id()] + actives.at(j)) * deriv2Data.block(3*j,0,3,f1ders.dim().second);
                }

            for (index_t i = g1System.get_numBoundaryVertexFunctions()[4]; i < g1System.get_numBoundaryVertexFunctions().last(); i++)
                for (index_t j = 0; j < actives.rows(); j++)
                {
                    f1ders += sol_sparse->at(i,g1System.get_numBasisFunctions()[geoEval.id()] + actives.at(j)) * derivData.block(2*j,0,2,f1ders.dim().second);
                    f1ders2 += sol_sparse->at(i,g1System.get_numBasisFunctions()[geoEval.id()] + actives.at(j)) * deriv2Data.block(3*j,0,3,f1ders.dim().second);
                }
        }
        else
        {
            for (index_t i = 0; i < g1System.get_numInterfaceFunctions().last(); i++)
                for (index_t j = 0; j < actives.rows(); j++)
                {
                    f1ders += sol_sparse->at(i,g1System.get_numBasisFunctionsInterface()[geoEval.id()] + actives.at(j)) * derivData.block(2*j,0,2,f1ders.dim().second);
                    f1ders2 += sol_sparse->at(i,g1System.get_numBasisFunctionsInterface()[geoEval.id()] + actives.at(j)) * deriv2Data.block(3*j,0,3,f1ders.dim().second);
                }


            for (index_t i = g1System.get_numInterfaceFunctions().last(); i < sol_sparse->rows() - 1; i++)
                for (index_t j = 0; j < actives.rows(); j++)
                {
                    f1ders += sol_sparse->at(i,g1System.get_numBasisFunctions()[geoEval.id()] + actives.at(j)) * derivData.block(2*j,0,2,f1ders.dim().second);
                    f1ders2 += sol_sparse->at(i,g1System.get_numBasisFunctionsInterface()[geoEval.id()] + actives.at(j)) * deriv2Data.block(3*j,0,3,f1ders.dim().second);
                }

        }


        for (index_t j = 0; j < actives.rows(); j++)
        {
            f1ders += sol_sparse->at(sol_sparse->rows()-1,g1System.get_numBasisFunctions()[geoEval.id()] + actives.at(j)) * derivData.block(2*j,0,2,f1ders.dim().second);
            f1ders2 += sol_sparse->at(sol_sparse->rows()-1,g1System.get_numBasisFunctions()[geoEval.id()] + actives.at(j)) * deriv2Data.block(3*j,0,3,f1ders.dim().second);
        }


        // get the gradients to columns
        f1ders.resize(quNodes.rows(), quNodes.cols() );
        //f1ders2.transposeInPlace();

        // Evaluate second function (defined of physical domain)
        geoEval.evaluateAt(quNodes);
        _func2.deriv2_into(geoEval.values(), f2ders2);
    }

    // assemble on element
    inline T compute(gsDomainIterator<T>    & ,
                     gsGeometryEvaluator<T> & geoEval,
                     gsVector<T> const      & quWeights,
                     T & accumulated)
    {
        T sum(0.0);
        for (index_t k = 0; k < quWeights.rows(); ++k) // loop over quadrature nodes
        {
            // Transform the gradients
            geoEval.transformDeriv2Hgrad(k, f1ders, f1ders2, f1pders2);
            f1pders2.transposeInPlace();
            //if ( f2Param )
            //f2ders.col(k)=geoEval.gradTransforms().block(0, k*d,d,d) * f2ders.col(k);// to do: generalize
            short_t parDim = geoEval.parDim();
            int rest = f1pders2.rows()-parDim;
            const T weight = quWeights[k] *  geoEval.measure(k);
            sum += weight * ((f1pders2.topRows(parDim) - f2ders2.col(k).topRows(parDim)).squaredNorm() +
                2*(f1pders2.bottomRows(rest) - f2ders2.col(k).bottomRows(rest)).squaredNorm());
        }
        accumulated += sum;
        return sum;
    }

private:


    gsMatrix<T> f1ders, f1ders2, f2ders2;
    gsMatrix<T> f1pders2;

};

}






