/** @file gsVisitorG1Basis.h

    @brief Visitor for the G1 Basis.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): P. Weinmüller
*/

#pragma once

namespace gismo
{
    template <short_t d, class T>
    class gsC1ArgyrisVertexBasisProjectionVisitor
    {
    public:

    gsC1ArgyrisVertexBasisProjectionVisitor()
    {
    }

    void initialize(const gsBasis<T>       & basis, //
                    gsQuadRule<T>    & rule)
    {
        gsVector<index_t> numQuadNodes( basis.dim() );
        for (int i = 0; i < basis.dim(); ++i) // to do: improve
            numQuadNodes[i] = basis.degree(i) + 1;

        // Setup Quadrature
        rule = gsGaussRule<T>(numQuadNodes);// NB!

        localMat.resize(6);
        localRhs.resize(6);

        rhsVals.resize(6);
    }

    // Evaluate on element.
    inline void evaluate(const gsGeometry<T>    & geo,
                         gsBasis<T>       & basis,
                         std::vector<gsBSplineBasis<T>>       & basis_plus,
                         std::vector<gsBSplineBasis<T>>       & basis_minus,
                         std::vector<gsBSplineBasis<T>>       & basis_geo,
                         gsApproxGluingData<d, T> approxGluingData,
                         gsMatrix<T>            & quNodes,
                         const real_t & sigma,
                         const std::vector<bool> & kindOfEdge,
                         const gsOptionList optionList)
    {
        // Compute the active basis functions
        // Assumes actives are the same for all quadrature points on the elements
        basis.active_into(quNodes.col(0), actives);

        // Evaluate basis functions on element
        basis.eval_into(quNodes, basisData);

        numActive = actives.rows();

        // Computing the basis functions at the vertex
        gsMatrix<> Phi(6,6);
        Phi.setIdentity();

        Phi.row(1) *= sigma;
        Phi.row(2) *= sigma;
        Phi.row(3) *= sigma * sigma;
        Phi.row(4) *= sigma * sigma;
        Phi.row(5) *= sigma * sigma;

        // Computing c, c+ and c-
        // Point zero
        gsMatrix<> zero;
        zero.setZero(2,1);

        std::vector<gsMatrix<>> c_0, c_1;
        std::vector<gsMatrix < >> c_0_plus, c_1_plus, c_2_plus;
        std::vector<gsMatrix < >> c_0_plus_deriv, c_1_plus_deriv, c_2_plus_deriv;
        std::vector<gsMatrix < >> c_0_minus, c_1_minus;
        for (index_t i = 0; i < 2; i++) // i == 0 == u , i == 1 == v
        {
            gsMatrix<> b_0, b_1;
            gsMatrix<> b_0_plus, b_1_plus, b_2_plus;
            gsMatrix<> b_0_plus_deriv, b_1_plus_deriv, b_2_plus_deriv;
            gsMatrix<> b_0_minus, b_1_minus;

            gsBSplineBasis<T> bsp_temp = dynamic_cast<gsBSplineBasis<> & >(basis_geo[1-i]);
            real_t p = bsp_temp.degree();
            real_t h_geo = bsp_temp.knots().at(p + 1);

            basis_geo[1-i].evalSingle_into(0, quNodes.row(i),b_0); // first
            basis_geo[1-i].evalSingle_into(1, quNodes.row(i),b_1); // second

            basis_plus[i].evalSingle_into(0, quNodes.row(i),b_0_plus);
            basis_plus[i].evalSingle_into(1, quNodes.row(i),b_1_plus);
            basis_plus[i].evalSingle_into(2, quNodes.row(i),b_2_plus);

            basis_plus[i].derivSingle_into(0, quNodes.row(i),b_0_plus_deriv);
            basis_plus[i].derivSingle_into(1, quNodes.row(i),b_1_plus_deriv);
            basis_plus[i].derivSingle_into(2, quNodes.row(i),b_2_plus_deriv);

            basis_minus[i].evalSingle_into(0, quNodes.row(i),b_0_minus);
            basis_minus[i].evalSingle_into(1, quNodes.row(i),b_1_minus);

            c_0.push_back(b_0 + b_1);
            c_1.push_back((h_geo / p) * b_1);

            c_0_minus.push_back(b_0_minus + b_1_minus);
            c_1_minus.push_back(h_geo/ (p-1) * b_1_minus);

            gsMatrix<> der_b_1_plus_0, der2_b_1_plus_0, der2_b_2_plus_0;
            basis_plus[i].derivSingle_into(1, zero.row(i), der_b_1_plus_0);
            basis_plus[i].deriv2Single_into(1, zero.row(i), der2_b_1_plus_0);
            basis_plus[i].deriv2Single_into(2, zero.row(i), der2_b_2_plus_0);

            real_t factor_c_1_plus = 1/der_b_1_plus_0(0,0);
            real_t factor2_c_1_plus = -der2_b_1_plus_0(0,0)/(der_b_1_plus_0(0,0)*der2_b_2_plus_0(0,0));
            real_t factor_c_2_plus = 1/der2_b_2_plus_0(0,0);

            c_0_plus.push_back(b_0_plus + b_1_plus + b_2_plus);
            c_1_plus.push_back(factor_c_1_plus * b_1_plus + factor2_c_1_plus * b_2_plus);
            c_2_plus.push_back(factor_c_2_plus * b_2_plus );

            c_0_plus_deriv.push_back(b_0_plus_deriv + b_1_plus_deriv + b_2_plus_deriv);
            c_1_plus_deriv.push_back(factor_c_1_plus * b_1_plus_deriv + factor2_c_1_plus * b_2_plus_deriv);
            c_2_plus_deriv.push_back(factor_c_2_plus * b_2_plus_deriv);

            // TODO IF CASE
            /*
            if ( p == 3)
            {
                // WORKS ONLY FOR p=3 AND r=1
                c_0_plus.push_back(b_0_plus + b_1_plus + b_2_plus);
                c_1_plus.push_back((h_geo / p) * (b_1_plus + 3 * b_2_plus));
                c_2_plus.push_back((h_geo * h_geo / (p * (p - 1))) * 2 * b_2_plus);

                c_0_plus_deriv.push_back(b_0_plus_deriv + b_1_plus_deriv + b_2_plus_deriv);
                c_1_plus_deriv.push_back((h_geo / p) * (b_1_plus_deriv + 3 * b_2_plus_deriv));
                c_2_plus_deriv.push_back((h_geo * h_geo / (p * (p - 1))) * 2 * b_2_plus_deriv);
            }
            else
            {
                c_0_plus.push_back(b_0_plus + b_1_plus + b_2_plus);
                c_1_plus.push_back((h_geo / p) * (b_1_plus + 2 * b_2_plus));
                c_2_plus.push_back((h_geo * h_geo / (p * (p - 1))) * b_2_plus);

                c_0_plus_deriv.push_back(b_0_plus_deriv + b_1_plus_deriv + b_2_plus_deriv);
                c_1_plus_deriv.push_back((h_geo / p) * (b_1_plus_deriv + 2 * b_2_plus_deriv));
                c_2_plus_deriv.push_back((h_geo * h_geo / (p * (p - 1))) * b_2_plus_deriv);

            }
         */
        }

        std::vector<gsMatrix<>> alpha, beta, alpha_0, beta_0, alpha_deriv, beta_deriv;

        gsMatrix < T > temp_mat;
        if (kindOfEdge[0])
        {
            approxGluingData.alphaS(0).eval_into(quNodes.row(0),temp_mat); // 1-dir == PatchID
            alpha.push_back(temp_mat); // u

            approxGluingData.alphaS(0).eval_into(zero.row(0),temp_mat); // 1-dir == PatchID
            alpha_0.push_back(temp_mat); // u

            approxGluingData.alphaS(0).deriv_into(zero.row(0),temp_mat); // 1-dir == PatchID
            alpha_deriv.push_back(temp_mat); // u

            approxGluingData.betaS(0).eval_into(quNodes.row(0),temp_mat); // 1-dir == PatchID
            beta.push_back(temp_mat); // u

            approxGluingData.betaS(0).eval_into(zero.row(0),temp_mat); // 1-dir == PatchID
            beta_0.push_back(temp_mat); // u

            approxGluingData.betaS(0).deriv_into(zero.row(0),temp_mat); // 1-dir == PatchID
            beta_deriv.push_back(temp_mat); // u
        }
        else
        {
            temp_mat.setOnes(1, quNodes.cols());
            alpha.push_back(temp_mat); // u

            temp_mat.setOnes(1, zero.cols());
            alpha_0.push_back(temp_mat); // u

            temp_mat.setZero(1, zero.cols());
            alpha_deriv.push_back(temp_mat); // u

            temp_mat.setZero(1, quNodes.cols());
            beta.push_back(temp_mat); // u

            temp_mat.setZero(1, zero.cols());
            beta_0.push_back(temp_mat); // u

            temp_mat.setZero(1, zero.cols());
            beta_deriv.push_back(temp_mat); // u
        }



        if (kindOfEdge[1]) {
            approxGluingData.alphaS(1).eval_into(quNodes.row(1), temp_mat); // 1-dir == PatchID
            alpha.push_back(temp_mat); // v

            approxGluingData.alphaS(1).eval_into(zero.row(0), temp_mat); // 1-dir == PatchID
            alpha_0.push_back(temp_mat); // v

            approxGluingData.alphaS(1).deriv_into(zero.row(0), temp_mat); // 1-dir == PatchID
            alpha_deriv.push_back(temp_mat); // v

            approxGluingData.betaS(1).eval_into(quNodes.row(1), temp_mat); // 1-dir == PatchID
            beta.push_back(temp_mat); // v

            approxGluingData.betaS(1).eval_into(zero.row(0), temp_mat); // 1-dir == PatchID
            beta_0.push_back(temp_mat); // v

            approxGluingData.betaS(1).deriv_into(zero.row(0), temp_mat); // 1-dir == PatchID
            beta_deriv.push_back(temp_mat); // v
        }
        else
        {
            temp_mat.setOnes(1, quNodes.cols());
            alpha.push_back(temp_mat); // u

            temp_mat.setOnes(1, zero.cols());
            alpha_0.push_back(temp_mat); // u

            temp_mat.setZero(1, zero.cols());
            alpha_deriv.push_back(temp_mat); // u

            temp_mat.setZero(1, quNodes.cols());
            beta.push_back(temp_mat); // u

            temp_mat.setZero(1, zero.cols());
            beta_0.push_back(temp_mat); // u

            temp_mat.setZero(1, zero.cols());
            beta_deriv.push_back(temp_mat); // u
        }

        // Geo data:
        gsMatrix<> geo_jac = geo.jacobian(zero);
        gsMatrix<T> geo_der2 = geo.deriv2(zero);

        // Compute dd^^(i_k) and dd^^(i_k-1)
        gsMatrix<> dd_ik_plus, dd_ik_minus;
        gsMatrix<> dd_ik_minus_deriv, dd_ik_plus_deriv;
        dd_ik_minus = -1/(alpha_0[0](0,0)) * (geo_jac.col(1) +
                                              beta_0[0](0,0) * geo_jac.col(0));

        dd_ik_plus = 1/(alpha_0[1](0,0)) * (geo_jac.col(0) +
                                            beta_0[1](0,0) * geo_jac.col(1));

        gsMatrix<> geo_deriv2_12(2,1), geo_deriv2_11(2,1), geo_deriv2_22(2,1);
        geo_deriv2_12.row(0) = geo_der2.row(2);
        geo_deriv2_12.row(1) = geo_der2.row(5);
        geo_deriv2_11.row(0) = geo_der2.row(0);
        geo_deriv2_11.row(1) = geo_der2.row(3);
        geo_deriv2_22.row(0) = geo_der2.row(1);
        geo_deriv2_22.row(1) = geo_der2.row(4);
        gsMatrix<> alpha_squared_u = alpha_0[0]*alpha_0[0];
        gsMatrix<> alpha_squared_v = alpha_0[1]*alpha_0[1];

        dd_ik_minus_deriv = -1/(alpha_squared_u(0,0)) * // N^2
                            ((geo_deriv2_12 + (beta_deriv[0](0,0) * geo_jac.col(0) +
                                               beta_0[0](0,0) * geo_deriv2_11))*alpha_0[0](0,0) -
                             (geo_jac.col(1) + beta_0[0](0,0) * geo_jac.col(0)) *
                             alpha_deriv[0](0,0));

        dd_ik_plus_deriv = 1/(alpha_squared_v(0,0)) *
                           ((geo_deriv2_12 + (beta_deriv[1](0,0) * geo_jac.col(1) +
                                              beta_0[1](0,0) * geo_deriv2_22))*alpha_0[1](0,0) -
                            (geo_jac.col(0) + beta_0[1](0,0) * geo_jac.col(1)) *
                            alpha_deriv[1](0,0));

/*
        gsInfo << "Transversal\n";
        if (kindOfEdge[1])
            gsInfo << "dd_ik_minus_deriv " << dd_ik_minus_deriv << "\n";
        if (kindOfEdge[0])
            gsInfo << "dd_ik_minus_deriv " << dd_ik_plus_deriv << "\n";

        gsInfo << geo_jac.col(0) << " : " << geo_jac.col(1) << "\n";
*/
        // Comupute d_(0,0)^(i_k), d_(1,0)^(i_k), d_(0,1)^(i_k), d_(1,1)^(i_k) ; i_k == 2
        std::vector<gsMatrix<>> d_ik;
        d_ik.push_back(Phi.col(0));
        d_ik.push_back(Phi.block(0,1,6,2) * geo_jac.col(0) ); // deriv into u
        d_ik.push_back(Phi.block(0,1,6,2) * geo_jac.col(1) ); // deriv into v
        d_ik.push_back((geo_jac(0,0) * Phi.col(3) + geo_jac(1,0) * Phi.col(4))*geo_jac(0,1) +
                       (geo_jac(0,0) * Phi.col(4) + geo_jac(1,0) * Phi.col(5))*geo_jac(1,1) +
                       Phi.block(0,1,6,1) * geo_der2.row(2) +
                       Phi.block(0,2,6,1) * geo_der2.row(5)); // Hessian

        // Compute d_(*,*)^(il,ik)
        std::vector<gsMatrix<>> d_ilik_minus, d_ilik_plus;
        d_ilik_minus.push_back(Phi.col(0));
        d_ilik_minus.push_back(Phi.block(0,1,6,2) * geo_jac.col(0));
        d_ilik_minus.push_back((geo_jac(0,0) * Phi.col(3) + geo_jac(1,0) * Phi.col(4))*geo_jac(0,0) +
                               (geo_jac(0,0) * Phi.col(4) + geo_jac(1,0) * Phi.col(5))*geo_jac(1,0) +
                               Phi.block(0,1,6,1) * geo_der2.row(0) +
                               Phi.block(0,2,6,1) * geo_der2.row(3));
        d_ilik_minus.push_back(Phi.block(0,1,6,2) * dd_ik_minus);
        d_ilik_minus.push_back((geo_jac(0,0) * Phi.col(3) + geo_jac(1,0) * Phi.col(4))*dd_ik_minus(0,0) +
                               (geo_jac(0,0) * Phi.col(4) + geo_jac(1,0) * Phi.col(5))*dd_ik_minus(1,0) +
                               Phi.block(0,1,6,1) * dd_ik_minus_deriv.row(0) +
                               Phi.block(0,2,6,1) * dd_ik_minus_deriv.row(1));

        d_ilik_plus.push_back(Phi.col(0));
        d_ilik_plus.push_back(Phi.block(0,1,6,2) * geo_jac.col(1));
        d_ilik_plus.push_back((geo_jac(0,1) * Phi.col(3) + geo_jac(1,1) * Phi.col(4))*geo_jac(0,1) +
                              (geo_jac(0,1) * Phi.col(4) + geo_jac(1,1) * Phi.col(5))*geo_jac(1,1) +
                              Phi.block(0,1,6,1) * geo_der2.row(1) +
                              Phi.block(0,2,6,1) * geo_der2.row(4));
        d_ilik_plus.push_back(Phi.block(0,1,6,2) * dd_ik_plus);
        d_ilik_plus.push_back((geo_jac(0,1) * Phi.col(3) + geo_jac(1,1) * Phi.col(4))*dd_ik_plus(0,0) +
                              (geo_jac(0,1) * Phi.col(4) + geo_jac(1,1) * Phi.col(5))*dd_ik_plus(1,0) +
                              Phi.block(0,1,6,1) * dd_ik_plus_deriv.row(0) +
                              Phi.block(0,2,6,1) * dd_ik_plus_deriv.row(1));

        for (index_t i = 0; i < 6; i++)
        {

            rhsVals.at(i) = d_ilik_minus.at(0)(i,0) * (c_0_plus.at(0).cwiseProduct(c_0.at(1)) -
                                                       beta[0].cwiseProduct(c_0_plus_deriv.at(0).cwiseProduct(c_1.at(1)))) +
                            d_ilik_minus.at(1)(i,0) * (c_1_plus.at(0).cwiseProduct(c_0.at(1)) -
                                                       beta[0].cwiseProduct(c_1_plus_deriv.at(0).cwiseProduct(c_1.at(1)))) +
                            d_ilik_minus.at(2)(i,0) * (c_2_plus.at(0).cwiseProduct(c_0.at(1)) -
                                                       beta[0].cwiseProduct(c_2_plus_deriv.at(0).cwiseProduct(c_1.at(1)))) -
                            d_ilik_minus.at(3)(i,0) * alpha[0].cwiseProduct(c_0_minus.at(0).cwiseProduct(c_1.at(1))) -
                            d_ilik_minus.at(4)(i,0) * alpha[0].cwiseProduct(c_1_minus.at(0).cwiseProduct(c_1.at(1))); // f*_(ik-1,ik)

            //if (kindOfEdge[0])
                //rhsVals.at(i).setZero();

            //if (!kindOfEdge[1])
                rhsVals.at(i) += d_ilik_plus.at(0)(i,0) * (c_0_plus.at(1).cwiseProduct(c_0.at(0)) -
                                                       beta[1].cwiseProduct(c_0_plus_deriv.at(1).cwiseProduct(c_1.at(0)))) +
                             d_ilik_plus.at(1)(i,0) * (c_1_plus.at(1).cwiseProduct(c_0.at(0)) -
                                                       beta[1].cwiseProduct(c_1_plus_deriv.at(1).cwiseProduct(c_1.at(0)))) +
                             d_ilik_plus.at(2)(i,0) * (c_2_plus.at(1).cwiseProduct(c_0.at(0)) -
                                                       beta[1].cwiseProduct(c_2_plus_deriv.at(1).cwiseProduct(c_1.at(0)))) +
                             d_ilik_plus.at(3)(i,0) * alpha[1].cwiseProduct(c_0_minus.at(1).cwiseProduct(c_1.at(0))) +
                             d_ilik_plus.at(4)(i,0) * alpha[1].cwiseProduct(c_1_minus.at(1).cwiseProduct(c_1.at(0))); // f*_(ik+1,ik)

            rhsVals.at(i) -= d_ik.at(0)(i,0) * c_0.at(0).cwiseProduct(c_0.at(1)) + d_ik.at(2)(i,0) * c_0.at(0).cwiseProduct(c_1.at(1)) +
                             d_ik.at(1)(i,0) * c_1.at(0).cwiseProduct(c_0.at(1)) + d_ik.at(3)(i,0) * c_1.at(0).cwiseProduct(c_1.at(1)); // f*_(ik)

            localMat.at(i).setZero(numActive, numActive);
            localRhs.at(i).setZero(numActive, rhsVals.at(i).rows());//multiple right-hand sides
        }



    } // evaluate

    inline void assemble(gsDomainIterator<T>    & element,
                         const gsVector<T>      & quWeights)
    {
        gsMatrix<T> & basisVals  = basisData;
        for (index_t i = 0; i < 6; i++)
        {
            // ( u, v)
            localMat.at(i).noalias() =
                    basisData * quWeights.asDiagonal() * basisData.transpose();

            for (index_t k = 0; k < quWeights.rows(); ++k) // loop over quadrature nodes
            {
                // Multiply weight by the geometry measure
                const T weight = quWeights[k];

                localRhs.at(i).noalias() += weight * (basisVals.col(k) * rhsVals.at(i).col(k).transpose());
            }
        }
    }

    inline void localToGlobal(const int patchIndex,
                              const std::vector<gsMatrix<T> >    & eliminatedDofs,
                              std::vector< gsSparseSystem<T> >     & system)
    {
        gsMatrix<index_t> actives_temp;
        for (unsigned i = 0; i < system.size(); i++) // 6
        {
            // Map patch-local DoFs to global DoFs
            system.at(i).mapColIndices(actives, patchIndex, actives_temp);
            // Add contributions to the system matrix and right-hand side
            system.at(i).push(localMat.at(i), localRhs.at(i), actives_temp, eliminatedDofs[i], 0, 0);
        }
    }

    protected:
        gsMatrix<index_t> actives;
        gsMatrix<T> basisData;
        index_t numActive;

    protected:
        // Local values of the right hand side
        std::vector< gsMatrix<T> >  rhsVals;

    protected:
        // Local matrices
        std::vector< gsMatrix<T> > localMat;
        std::vector< gsMatrix<T> > localRhs;

    }; // class gsVisitorG1BasisVertex
} // namespace gismo