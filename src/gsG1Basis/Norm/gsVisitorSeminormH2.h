/** @file gsNormL2.h

    @brief Computes the Semi H2 norm, needs for the parallel computing.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): P. Weinmüller, A. Farahat
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
        evFlags = NEED_MEASURE|NEED_VALUE|NEED_GRAD_TRANSFORM|NEED_2ND_DER|NEED_JACOBIAN;
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
        gsMatrix<T> derivData, deriv2Data;

        qN = quNodes;

        basis.active_into(quNodes.col(0), actives);

        // Evaluate basis functions on element
        basis.deriv_into(quNodes,derivData);
        basis.deriv2_into(quNodes,deriv2Data);
        geoEval.evaluateAt(quNodes);

        gsMatrix<T> Jk = geoEval.jacobian(0);

        // Surface error computation -> if (paramDim + 1 == targetDim)
        if( Jk.dim().second + 1 == Jk.dim().first )
        {
            basis.evalAllDers_into(quNodes, 2, basisData);
            basisGrads = basisData[1];
            basis2ndDerivs = basisData[2];
            numActive = actives.rows();
        }

        f1ders.setZero(2,actives.rows());
        f1ders2.setZero(3,actives.rows());
        for (index_t i = 0; i < sol_sparse->rows(); i++)
            for (index_t j = 0; j < actives.rows(); j++)
            {
                f1ders += sol_sparse->at(i,numBasisFunctions[geoEval.id()] + actives.at(j)) * derivData.block(2*j,0,2,f1ders.dim().second);
                f1ders2 += sol_sparse->at(i,numBasisFunctions[geoEval.id()] + actives.at(j)) * deriv2Data.block(3*j,0,3,f1ders.dim().second);
            }


        // get the gradients to columns
        f1ders.resize(quNodes.rows(), quNodes.cols() );
        //f1ders2.transposeInPlace();

        // Evaluate second function (defined of physical domain)

        _func2.eval_into(geoEval.values(), f2ders2);

    }

    // assemble on element
    inline T compute(gsDomainIterator<T>    & ,
                     gsGeometryEvaluator<T> & geoEval,
                     gsVector<T> const      & quWeights,
                     T & accumulated,
                     const gsGeometry<T>    & geo)
    {
        T sum(0.0);
        for (index_t k = 0; k < quWeights.rows(); ++k) // loop over quadrature nodes
        {

            gsMatrix<T> Jk = geoEval.jacobian(k);

            T weight = quWeights[k];

            // Surface error computation -> if (paramDim + 1 == targetDim)
            if( Jk.dim().second + 1 == Jk.dim().first )
            {
                gsMatrix<T> geoMapDeriv1 = geo.deriv(qN); // First derivative of the geometric mapping with respect to the parameter coordinates
                gsMatrix<T> geoMapDeriv2 = geo.deriv2(qN); // Second derivative of the geometric mapping with respect to the parameter coordinates


//            FIRST FUNDAMENTAL FORM: G = J^T * J
//
//            G = | G11   G12|
//                | G21   G22|
//
//            INVERSE OF THE FIRST FUNDAMENTAL FORM
//
//                      1    | G22  -G12|      1
//            G^-1 = ------- |          | = ------- G* ^-1
//                    det(G) | -G21  G11|    det(G)

                // First fundamental form
                real_t G11 = ( geoMapDeriv1(0, k) * (geoMapDeriv1(0, k)) +
                               geoMapDeriv1(2, k) * (geoMapDeriv1(2, k)) +
                               geoMapDeriv1(4, k) * (geoMapDeriv1(4, k)));

//          G12 = G21
                real_t G12 = ( geoMapDeriv1(0, k) * (geoMapDeriv1(1, k)) +
                                geoMapDeriv1(2, k) * (geoMapDeriv1(3, k)) +
                                geoMapDeriv1(4, k) * (geoMapDeriv1(5, k)));

                real_t G22 = ( geoMapDeriv1(1, k) * (geoMapDeriv1(1, k)) +
                                geoMapDeriv1(3, k) * (geoMapDeriv1(3, k)) +
                                geoMapDeriv1(5, k) * (geoMapDeriv1(5, k)));

                // Derivative of the first fundamental form
                real_t DuG11 = 2 * ( geoMapDeriv2(0, k) * (geoMapDeriv1(0, k)) +
                                      geoMapDeriv2(3, k) * (geoMapDeriv1(2, k)) +
                                      geoMapDeriv2(6, k) * (geoMapDeriv1(4, k)) );


                real_t DvG11 = 2 * ( geoMapDeriv2(2, k) * (geoMapDeriv1(0, k)) +
                                      geoMapDeriv2(5, k) * (geoMapDeriv1(2, k)) +
                                      geoMapDeriv2(8, k) * (geoMapDeriv1(4, k)) );

//          DuG12 = DuG21
                real_t DuG12 = ( geoMapDeriv2(0, k) * (geoMapDeriv1(1, k)) +
                                  geoMapDeriv2(2, k) * (geoMapDeriv1(0, k)) +
                                  geoMapDeriv2(3, k) * (geoMapDeriv1(3, k)) +
                                  geoMapDeriv2(5, k) * (geoMapDeriv1(2, k)) +
                                  geoMapDeriv2(6, k) * (geoMapDeriv1(5, k)) +
                                  geoMapDeriv2(8, k) * (geoMapDeriv1(4, k)) );

//          DvG12 = DvG21
                real_t DvG21 = ( geoMapDeriv2(2, k) * (geoMapDeriv1(1, k)) +
                                  geoMapDeriv2(1, k) * (geoMapDeriv1(0, k)) +
                                  geoMapDeriv2(5, k) * (geoMapDeriv1(3, k)) +
                                  geoMapDeriv2(4, k) * (geoMapDeriv1(2, k)) +
                                  geoMapDeriv2(8, k) * (geoMapDeriv1(5, k)) +
                                  geoMapDeriv2(7, k) * (geoMapDeriv1(4, k)) );


                real_t DuG22 = 2 * ( geoMapDeriv2(2, k) * (geoMapDeriv1(1, k)) +
                                      geoMapDeriv2(5, k) * (geoMapDeriv1(3, k)) +
                                      geoMapDeriv2(8, k) * (geoMapDeriv1(5, k)) );

                real_t DvG22 = 2 *( geoMapDeriv2(1, k) * (geoMapDeriv1(1, k)) +
                                     geoMapDeriv2(4, k) * (geoMapDeriv1(3, k)) +
                                     geoMapDeriv2(7, k) * (geoMapDeriv1(5, k)) );

                real_t detG = G11 * (G22) - G12 * (G12);

//          1 / sqrt( det( G ) )
                real_t sqrtDetG_inv;

//          1 / ( 2 * det( G )^( 3/2 ) )
                real_t sqrtDetG_inv_derivative;


                sqrtDetG_inv = 1 / sqrt( detG );
                sqrtDetG_inv_derivative = 1 / ( 2 * detG * sqrt( detG ) );

                real_t Du_SqrtDetGinv = sqrtDetG_inv_derivative * (
                                        2 * G12 * ( DuG12 )  -
                                        G22 * ( DuG11 ) -
                                        G11 * ( DuG22 ) );

                real_t Dv_SqrtDetGinv = sqrtDetG_inv_derivative * (
                                        2 * G12 * ( DvG21 )  -
                                        G22 * ( DvG11 ) -
                                        G11 * ( DvG22 ) );

//          div ( sqrt( det( G ) ) * ( 1 / det( G ) * G* ^-1 * grad( u ) ) )
                gsMatrix<T> surfParametricLaplace(1, 1);
                surfParametricLaplace.setZero();

                gsMatrix<T> exactLap(1, 1);
                exactLap.setZero();

                for(index_t i = 0; i < numActive; i++)
                {
//              1 / sqrt^4( det( G ) ) *
//              [
//              Du( 1 / sqrt( det( G ) ) ) * G* ^-1 * grad( u ) ) +
//              Dv( 1 / sqrt( det( G ) ) ) * G* ^-1 * grad( u ) ) +
//              1 / sqrt( det( G ) ) * ( div ( G* ^-1 * grad( u ) ) )
//              ]
                    surfParametricLaplace(0, 0) += ( Du_SqrtDetGinv * (
                                                    G22 * ( basisGrads( i * 2, k) ) -
                                                    G12 * ( basisGrads( i * 2 + 1, k) ) )
                                                    +
                                                    Dv_SqrtDetGinv * (
                                                    G11 * ( basisGrads( i * 2 + 1, k) ) -
                                                    G12 * ( basisGrads( i * 2, k) ) )
                                                    +
                                                    sqrtDetG_inv * (
                                                    DuG22 * ( basisGrads( i * 2, k) ) +
                                                    G22 * ( basis2ndDerivs( i * 3, k) ) -
                                                    DuG12 * ( basisGrads( i * 2 + 1, k) ) -
                                                    G12 * ( basis2ndDerivs( i * 3 + 2, k) ) -
                                                    G12 * ( basis2ndDerivs( i * 3 + 2, k) ) +
                                                    DvG11 * ( basisGrads( i * 2 + 1, k) ) +
                                                    G11 * ( basis2ndDerivs( i * 3 + 1, k) ) -
                                                    DvG21 * ( basisGrads( i * 2, k) ) )
                                                    ) * sqrtDetG_inv;

                }

                weight *= sqrt(detG);

                exactLap(0, 0) = f2ders2(0, k) + f2ders2(1, k) + f2ders2(2, k);
                sum += weight * (surfParametricLaplace - exactLap).squaredNorm();

            }
            else
            {
                // Transform the gradients
                geoEval.transformDeriv2Hgrad(k, f1ders, f1ders2, f1pders2);
                f1pders2.transposeInPlace();
                //if ( f2Param )
                //f2ders.col(k)=geoEval.gradTransforms().block(0, k*d,d,d) * f2ders.col(k);// to do: generalize
                short_t parDim = geoEval.parDim();
                int rest = f1pders2.rows() - parDim;
                weight *= geoEval.measure(k);
                sum += weight * ((f1pders2.topRows(parDim) - f2ders2.col(k).topRows(parDim)).squaredNorm() +
                    2 * (f1pders2.bottomRows(rest) - f2ders2.col(k).bottomRows(rest)).squaredNorm());
            }
        }
        accumulated += sum;
        return sum;
    }

private:


    gsMatrix<T> f1ders, f1ders2, f2ders2;
    gsMatrix<T> f1pders2;

    std::vector<gsMatrix<T> > basisData;
    gsMatrix<T> basisGrads;
    gsMatrix<T> basis2ndDerivs;
    index_t numActive;
    gsMatrix<T> qN;
};

}






