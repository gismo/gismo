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
        {
            numQuadNodes[i] = basis.degree(i) + 1;
        }
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

        f1ders.setZero(2,quNodes.cols());
        f1ders2.setZero(3,quNodes.cols());
        for (index_t i = 0; i < sol_sparse->rows(); i++)
            for (index_t j = 0; j < actives.rows(); j++)
            {
                f1ders += sol_sparse->at(i,numBasisFunctions[geoEval.id()] + actives.at(j)) * derivData.block(2*j,0,2,f1ders.dim().second);
                f1ders2 += sol_sparse->at(i,numBasisFunctions[geoEval.id()] + actives.at(j)) * deriv2Data.block(3*j,0,3,f1ders.dim().second);
            }


        // get the gradients to columns
//        f1ders.resize(quNodes.rows(), quNodes.cols() );

//        f1ders2.transposeInPlace();

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
                real_t G12 = (  geoMapDeriv1(0, k) * (geoMapDeriv1(1, k)) +
                                geoMapDeriv1(2, k) * (geoMapDeriv1(3, k)) +
                                geoMapDeriv1(4, k) * (geoMapDeriv1(5, k)));

                real_t G22 = (  geoMapDeriv1(1, k) * (geoMapDeriv1(1, k)) +
                                geoMapDeriv1(3, k) * (geoMapDeriv1(3, k)) +
                                geoMapDeriv1(5, k) * (geoMapDeriv1(5, k)));

                // Derivative of the first fundamental form
                real_t DuG11 = 2 * (  geoMapDeriv2(0, k) * (geoMapDeriv1(0, k)) +
                                      geoMapDeriv2(3, k) * (geoMapDeriv1(2, k)) +
                                      geoMapDeriv2(6, k) * (geoMapDeriv1(4, k)) );


                real_t DvG11 = 2 * (  geoMapDeriv2(2, k) * (geoMapDeriv1(0, k)) +
                                      geoMapDeriv2(5, k) * (geoMapDeriv1(2, k)) +
                                      geoMapDeriv2(8, k) * (geoMapDeriv1(4, k)) );

//          DuG12 = DuG21
                real_t DuG12 = (  geoMapDeriv2(0, k) * (geoMapDeriv1(1, k)) +
                                  geoMapDeriv2(2, k) * (geoMapDeriv1(0, k)) +
                                  geoMapDeriv2(3, k) * (geoMapDeriv1(3, k)) +
                                  geoMapDeriv2(5, k) * (geoMapDeriv1(2, k)) +
                                  geoMapDeriv2(6, k) * (geoMapDeriv1(5, k)) +
                                  geoMapDeriv2(8, k) * (geoMapDeriv1(4, k)) );

//          DvG12 = DvG21
                real_t DvG21 = (  geoMapDeriv2(2, k) * (geoMapDeriv1(1, k)) +
                                  geoMapDeriv2(1, k) * (geoMapDeriv1(0, k)) +
                                  geoMapDeriv2(5, k) * (geoMapDeriv1(3, k)) +
                                  geoMapDeriv2(4, k) * (geoMapDeriv1(2, k)) +
                                  geoMapDeriv2(8, k) * (geoMapDeriv1(5, k)) +
                                  geoMapDeriv2(7, k) * (geoMapDeriv1(4, k)) );


                real_t DuG22 = 2 * (  geoMapDeriv2(2, k) * (geoMapDeriv1(1, k)) +
                                      geoMapDeriv2(5, k) * (geoMapDeriv1(3, k)) +
                                      geoMapDeriv2(8, k) * (geoMapDeriv1(5, k)) );

                real_t DvG22 = 2 *(  geoMapDeriv2(1, k) * (geoMapDeriv1(1, k)) +
                                     geoMapDeriv2(4, k) * (geoMapDeriv1(3, k)) +
                                     geoMapDeriv2(7, k) * (geoMapDeriv1(5, k)) );


                gsMatrix<T> G = Jk.transpose() * Jk;
                gsMatrix<T> G_inv = G.cramerInverse();

                real_t detG = G11 * (G22) - G12 * (G12);
//                real_t detG_inv = 1 / ( detG );

//                real_t Du_detG_Inv = detG_inv * detG_inv * ( 2 * G12 * DuG12 -
//                                                                G22 * DuG11 -
//                                                                G11 * DuG22 );
//
//                real_t Dv_detG_Inv = detG_inv * detG_inv * ( 2 * G12 * DvG21 -
//                                                                G22 * DvG11 -
//                                                                G11 * DvG22 );

                gsMatrix<T> sqrt4DetG_inv(1, 1);

//          1 / sqrt( det( G ) )
                gsMatrix<T> sqrtDetG_inv(1, 1);

//          1 / ( 2 * det( G )^( 3/2 ) )
                gsMatrix<T> sqrtDetG_inv_derivative(1, 1);

//          Creating the vector of the determinant of the first fundamental form

                sqrtDetG_inv(0, 0) = 1 / sqrt( detG );
                sqrt4DetG_inv(0, 0) = 1 / ( sqrt( sqrt( detG ) ) );
                sqrtDetG_inv_derivative(0, 0) = 1 / ( 2 * detG * sqrt( detG ) );

                gsMatrix<T> Du_SqrtDetGinv = sqrtDetG_inv_derivative * (
                                                            2 * G12 * ( DuG12 )  -
                                                                G22 * ( DuG11 ) -
                                                                G11 * ( DuG22 ) );

                gsMatrix<T> Dv_SqrtDetGinv = sqrtDetG_inv_derivative * (
                                                            2 * G12 * ( DvG21 )  -
                                                                G22 * ( DvG11 ) -
                                                                G11 * ( DvG22 ) );
////          Computing the divergence of the first fundamental form
//                gsMatrix<> div_G_inv(1, 2);
//                div_G_inv.setZero();
//
//                div_G_inv(0, 0) = Du_detG_Inv * G22 - Dv_detG_Inv * G12;
//                div_G_inv(0, 0) += ( ( DuG22 - DvG21 ) / detG );
//
//                div_G_inv(0, 1) = Dv_detG_Inv * G11 - Du_detG_Inv * G12;
//                div_G_inv(0, 1) += ( ( DvG11 - DuG12 ) / detG );
//
//
////          Computing the divergence of the jacobian
//                gsMatrix<> div_Jk_transpose(1, 3);
//                div_G_inv.setZero();
//
//                div_Jk_transpose(0, 0) = geoMapDeriv2(0, 0) + geoMapDeriv2(1, 0);
//                div_Jk_transpose(0, 1) = geoMapDeriv2(3, 0) + geoMapDeriv2(4, 0);
//                div_Jk_transpose(0, 2) = geoMapDeriv2(6, 0) + geoMapDeriv2(7, 0);

//                gsMatrix<T> surfParametricLaplace(6, 1);
//                surfParametricLaplace.setZero();
//
////          Computing the transformation of the hessian matrix from parameter to surface
//                gsMatrix<> grad_par_first = f1ders.col(k);
//
//                gsMatrix<> hessian_par_first(2, 2);
//                hessian_par_first.setZero();
//                hessian_par_first(0, 0) = f1ders2( 0, k);
//                hessian_par_first(1, 1) = f1ders2( 1, k);
//                hessian_par_first(0, 1) = f1ders2( 2, k);
//                hessian_par_first(1, 0) = f1ders2( 2, k);
//
//
//// ---------------------------------------------------------------------------------------------------------------------
//
////          Computing the tranformartion of the gradient basis functions
////                gsMatrix<> hessian_fromGrad_first(3, 3);
////                hessian_fromGrad_first.setZero();
//                gsMatrix<> hessian_fromGrad_first = Jk * G_inv * div_G_inv.transpose() * grad_par_first.transpose() * Jk.transpose();
//                           hessian_fromGrad_first += Jk * G_inv * G_inv * grad_par_first * div_Jk_transpose;
//
//                hessian_fromGrad_first.setZero();
//
////          Computing the tranformation of the hessian basis functions
//                gsMatrix<> hessian_phys_first = Jk * G_inv * hessian_par_first * G_inv * Jk.transpose();
//
//                surfParametricLaplace(0, 0) += (hessian_fromGrad_first(0, 0) + hessian_phys_first(0, 0));
//                surfParametricLaplace(1, 0) += (hessian_fromGrad_first(1, 1) + hessian_phys_first(1, 1));
//                surfParametricLaplace(2, 0) += (hessian_fromGrad_first(2, 2) + hessian_phys_first(2, 2));
//                surfParametricLaplace(3, 0) += (hessian_fromGrad_first(0, 1) + hessian_phys_first(0, 1));
//                surfParametricLaplace(4, 0) += (hessian_fromGrad_first(0, 2) + hessian_phys_first(0, 2));
////                surfParametricLaplace(5, 0) += (hessian_fromGrad_first(1, 0) + hessian_phys_first(1, 0));
//                surfParametricLaplace(5, 0) += (hessian_fromGrad_first(1, 2) + hessian_phys_first(1, 2));
////                surfParametricLaplace(7, 0) += (hessian_fromGrad_first(2, 0) + hessian_phys_first(2, 0));
////                surfParametricLaplace(8, 0) += (hessian_fromGrad_first(2, 1) + hessian_phys_first(2, 1));



                gsMatrix<T> surfParametricLaplace(1, 1);
                surfParametricLaplace.setZero();

                surfParametricLaplace = ( (    G22 * ( DvG11 ) -
                                            2 * G12 * ( DvG21 ) +
                                                G11 * ( DvG22 ) ) * (
                                                G12 * ( f1ders(0, k) ) -
                                                G11 * ( f1ders(1, k) ) ) -
                                                ( G22 * ( DuG11 ) -
                                            2 * G12 * ( DuG12 ) +
                                                G11 * ( DuG22 ) ) * (
                                                G22 * ( f1ders(0, k) ) -
                                                G12 * ( f1ders(1, k) ) ) +
                                            2 * ( G12 * ( G12 ) -
                                                G11 * ( G22 ) ) * (
                                                DvG21 * ( f1ders(0, k) ) -
                                                DvG11 * ( f1ders(1, k) ) -
                                                G11 * ( f1ders2(1, k) ) +
                                                G12 * ( f1ders2(2, k) ) ) +
                                            2 * ( G12 * ( G12 ) -
                                                G11 * ( G22 ) ) * (
                                                DuG12 * ( f1ders(1, k) ) -
                                                DuG22 * ( f1ders(0, k) ) +
                                                G12 * ( f1ders2(2, k) ) -
                                                G22 * ( f1ders2(0, k) ) ) ) * (
                                                sqrtDetG_inv_derivative );


                surfParametricLaplace *= sqrtDetG_inv;
                weight *= sqrt(detG);

                sum += weight * ( (surfParametricLaplace - f2ders2.col(k)).squaredNorm()
//                                +  (surfParametricLaplace.bottomRows(6) - f2ders2.col(k).bottomRows(6)).squaredNorm());
//                                +  2 * (surfParametricLaplace.bottomRows(3) - f2ders2.col(k).bottomRows(3)).squaredNorm()
                                );


            }
            else
            {
                // Transform the gradients
                geoEval.transformDeriv2Hgrad(k, f1ders, f1ders2, f1pders2);
                f1pders2.transposeInPlace();
                //if ( f2Param )
                //f2ders.col(k)=geoEval.gradTransforms().block(0, k*d,d,d) * f2ders.col(k);// to do: generalize
//                short_t parDim = geoEval.parDim();
//                int rest = f1pders2.rows() - parDim;
                weight *= geoEval.measure(k);

                //Using the equivalent norm of the laplacian
                gsMatrix<T> surfParametricLaplace(1, 1);
                surfParametricLaplace.setZero();
                surfParametricLaplace = f1pders2.row(0) + f1pders2.row(1);

                sum += weight * ((surfParametricLaplace - f2ders2.col(k)).squaredNorm() );

//                sum += weight * ((f1pders2.topRows(parDim) - f2ders2.col(k).topRows(parDim)).squaredNorm() +
//                    2 * (f1pders2.bottomRows(rest) - f2ders2.col(k).bottomRows(rest)).squaredNorm() );
            }
        }
        accumulated += sum;


        return sum;
    }

private:


    gsMatrix<T> f1ders, f1ders2, f2ders2;
    gsMatrix<T> f1pders2;
    gsMatrix<T> qN;
};

}






