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
class gsVisitorSeminormH1
{
public:

    gsVisitorSeminormH1()
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
                  const gsFunction<T>    & _func2,
                  const gsBasis<T> & basis,
                  const gsSparseMatrix<T> * sol_sparse,
                  const gsVector<T> & numBasisFunctions,
                  gsMatrix<T>            & quNodes)
    {
        gsMatrix<unsigned> actives;
        gsMatrix<T> bGrads;

        basis.active_into(quNodes.col(0), actives);

        // Evaluate basis functions on element
        basis.deriv_into(quNodes,bGrads);

        f1ders.setZero(2,actives.rows());
        for (index_t i = 0; i < sol_sparse->rows(); i++)
            for (index_t j = 0; j < actives.rows(); j++)
                f1ders += sol_sparse->at(i,numBasisFunctions[geoEval.id()] + actives.at(j)) * bGrads.block(2*j,0,2,f1ders.dim().second);

        // Evaluate second function
        geoEval.evaluateAt(quNodes);

        _func2.eval_into( geoEval.values(), f2ders);

//        _func2.deriv_into( f2param ? quNodes : geoEval.values() , f2ders); // Not working and useless

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

            T weight = quWeights[k];

            // Surface error computation -> if (paramDim + 1 == targetDim)
            if( Jk.dim().second +1 == Jk.dim().first )
            {
                gsMatrix<T> G = Jk.transpose() * Jk;
                gsMatrix<T> G_inv = G.cramerInverse();

                f1pders = Jk * G_inv * f1ders.col(k); // Computed gradient

                f2pders =  f2ders.col(k); // Exact gradient

                weight *= sqrt(G.determinant());
            }
            else
            {
                // Transform the gradients
                geoEval.transformGradients(k, f1ders, f1pders);

                // Transform the gradients, if func2 is defined on the parameter space (f2param = true)
                if(f2param)
                    geoEval.transformGradients(k, f2ders, f2pders);

                // old
                //if ( f2param )
                //f2ders.col(k)=geoEval.gradTransforms().block(0, k*d,d,d) * f2ders.col(k);// to do: generalize

                weight *= geoEval.measure(k);
            }

            // for each k: put the gradients into the columns (as in f1pders)
            gsMatrix<T> f2dersk = f2ders.col(k);
            f2dersk.resize(Jk.dim().first, 1); // pardim(), targetDim() // TODO

            sum += weight * (f1pders - f2dersk).squaredNorm();

        }
        accumulated += sum;
        return sum;
    }


protected:

    gsMatrix<T> f1ders, f2ders;
    gsMatrix<T> f1pders, f2pders; // f2pders only needed if f2param = true

    bool f2param;

};






}







