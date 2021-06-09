/** @file gsNormL2.h
    @brief Computes the Semi H1 norm, needs for the parallel computing.
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
    class gsVisitorG1Norm
    {
    public:

        gsVisitorG1Norm()
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
            md.flags = NEED_MEASURE|NEED_VALUE|NEED_GRAD_TRANSFORM|NEED_2ND_DER;
        }

        // Evaluate on element.
        void evaluate(const gsGeometry<T> & geometry,
                      const gsMultiBasis<T> * basisPtr,
                      const gsMultiPatch<T> * mpSolPtr,
                      const gsMatrix<T> * g1Sol,
                      const gsFunctionWithDerivatives<T> & exactSol,
                      gsG1MultiBasis<T> & g1MultiBasis,
                      gsMatrix<T> & quNodes)
        {
            md.points = quNodes;

            gsMatrix<index_t> actives;
            gsMatrix<T> bGrads;

            //basisPtr->basis(geometry.id()).active_into(quNodes.col(0), actives);

            // Evaluate basis functions on element
            //basisPtr->basis(geometry.id()).deriv_into(quNodes,bGrads);

            f1.setZero(1,quNodes.cols());
            f1ders.setZero(2,quNodes.cols());
            f1ders2.setZero(3,quNodes.cols());

            geometry.computeMap(md);

            // Interior solution
            f1 = mpSolPtr->patch(geometry.id()).eval(quNodes);
            f1ders = mpSolPtr->patch(geometry.id()).deriv(quNodes);
            f1ders2 = mpSolPtr->patch(geometry.id()).deriv2(quNodes);

            // G1 solution
            std::vector<gsMatrix<T>> eval_deriv_g1;
            //g1MultiBasis.eval_into(quNodes, eval_g1, geometry.id());
            g1MultiBasis.evalAllDers_into(quNodes, 2, eval_deriv_g1, geometry.id());

            for (index_t i = 0; i < eval_deriv_g1[0].rows(); i++)
            {
                f1 += eval_deriv_g1[0].row(i) * (*g1Sol)(i,0); // Only 1 dim solution!!!
                f1ders += eval_deriv_g1[1].block(2*i,0,2,quNodes.cols()) * (*g1Sol)(i,0);
                f1ders2 += eval_deriv_g1[2].block(3*i,0,3,quNodes.cols()) * (*g1Sol)(i,0);
                // TODO Add here bdy
            }

            //for (index_t j = 0; j < actives.rows(); j++)
            //    f1ders += sol_sparse->at(sol_sparse->rows() -1,g1System.get_numBasisFunctions()[geoEval.id()] + actives.at(j)) * bGrads.block(2*j,0,2,f1ders.dim().second);

            // Evaluate second function
            //geometry.evaluateAt(quNodes);
            exactSol.eval_into( f2param ? quNodes : md.values[0] , f2);
            exactSol.getDeriv().eval_into( f2param ? quNodes : md.values[0] , f2ders);
            exactSol.getDeriv2().eval_into( f2param ? quNodes : md.values[0] , f2ders2);


        }


        // assemble on element
        inline void compute(gsDomainIterator<T>    & geo,
                         gsVector<T> const      & quWeights,
                         T & accumulatedL2,
                         T & accumulatedH1,
                         T & accumulatedH2)
        {
            T sumL2(0.0);
            T sumH1(0.0);
            T sumH2(0.0);

            for (index_t k = 0; k < quWeights.rows(); ++k) // loop over quadrature nodes
            {
                // Transform the gradients
                transformGradients(md, k, f1ders, f1pders);
                //geometry.transformGradients(k, f1ders, f1pders);

                transformDeriv2Hgrad(md, k, f1ders, f1ders2, f1pders2);
                f1pders2.transposeInPlace();

                short_t parDim = 2;
                int rest = f1pders2.rows()-parDim;

                const T weight = quWeights[k] *  md.measure(k);

                // for each k: put the gradients into the columns (as in f1pders)
                //
                sumL2 += weight * (f1.col(k) - f2.col(k)).squaredNorm();
                sumH1 += weight * (f1pders - f2ders.col(k)).squaredNorm();
                sumH2 += weight * ((f1pders2.topRows(parDim) - f2ders2.col(k).topRows(parDim)).squaredNorm() +
                                   2*(f1pders2.bottomRows(rest) - f2ders2.col(k).bottomRows(rest)).squaredNorm());

            }
            accumulatedL2 += sumL2;
            accumulatedH1 += sumH1;
            accumulatedH2 += sumH2;
        }


    protected:

        gsMatrix<T> f1, f2, f1ders, f1ders2, f2ders, f2ders2;
        gsMatrix<T> f1pders, f1pders2; // f2pders only needed if f2param = true

        bool f2param;

        gsMapData<T> md;

    };






}