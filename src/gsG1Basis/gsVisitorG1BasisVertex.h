/** @file gsVisitorG1Basis.h

    @brief Visitor for the G1 Basis.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): P. Weinmüller
*/

#pragma once

#include "gsG1Basis/gsGluingData.h"

namespace gismo
{
template <class T>
class gsVisitorG1BasisVertex
{
public:

    gsVisitorG1BasisVertex()
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

        // Set Geometry evaluation flags
        md.flags = NEED_MEASURE ;

        localMat.resize(6);
        localRhs.resize(6);

        rhsVals.resize(6);
    }

    // Evaluate on element.
    inline void evaluate(gsBasis<T>       & basis, //
                         gsBasis<T>       & basis_geo,
                         std::vector<gsBSplineBasis<T>>       & basis_plus,
                         std::vector<gsBSplineBasis<T>>      & basis_minus,
                         const gsGeometry<T>    & geo, // patch
                         gsMatrix<T>            & quNodes,
                         std::vector<gsGluingData<T>>  & gluingData,
                         std::vector<bool> & isBoundary,
                         real_t sigma,
                         gsOptionList optionList)
    {
        md.points = quNodes;

        // Compute the active basis functions
        // Assumes actives are the same for all quadrature points on the elements
        basis.active_into(md.points.col(0), actives);

        // Evaluate basis functions on element
        basis.eval_into(md.points, basisData);

        // Compute geometry related values
        geo.computeMap(md);

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

            gsBSplineBasis<T> bsp_temp = dynamic_cast<gsBSplineBasis<> & >(basis_geo.component(i));
            real_t p = bsp_temp.maxDegree();
            real_t h_geo = bsp_temp.knots().at(p + 2);

            basis_geo.component(i).evalSingle_into(0, md.points.row(i),b_0); // first
            basis_geo.component(i).evalSingle_into(1, md.points.row(i),b_1); // second

            basis_plus[i].evalSingle_into(0, md.points.row(i),b_0_plus);
            basis_plus[i].evalSingle_into(1, md.points.row(i),b_1_plus);
            basis_plus[i].evalSingle_into(2, md.points.row(i),b_2_plus);

            basis_plus[i].derivSingle_into(0, md.points.row(i),b_0_plus_deriv);
            basis_plus[i].derivSingle_into(1, md.points.row(i),b_1_plus_deriv);
            basis_plus[i].derivSingle_into(2, md.points.row(i),b_2_plus_deriv);

            basis_minus[i].evalSingle_into(0, md.points.row(i),b_0_minus);
            basis_minus[i].evalSingle_into(1, md.points.row(i),b_1_minus);

            c_0.push_back(b_0 + b_1);
            c_1.push_back((h_geo / p) * b_1);

            c_0_minus.push_back(b_0_minus + b_1_minus);
            c_1_minus.push_back(h_geo/ (p-1) * b_1_minus);

            // TODO IF CASE
            // WORKS ONLY FOR p=3 AND r=1
            c_0_plus.push_back(b_0_plus + b_1_plus + b_2_plus);
            c_1_plus.push_back((h_geo / p) * (b_1_plus + 3 * b_2_plus));
            c_2_plus.push_back((h_geo * h_geo / (p * (p-1))) * 2 * b_2_plus);

            c_0_plus_deriv.push_back(b_0_plus_deriv + b_1_plus_deriv + b_2_plus_deriv);
            c_1_plus_deriv.push_back((h_geo / p) * (b_1_plus_deriv + 3 * b_2_plus_deriv));
            c_2_plus_deriv.push_back((h_geo * h_geo / (p * (p-1))) * 2 * b_2_plus_deriv);
        }

        // Compute dd^^(i_k) and dd^^(i_k-1)
        gsMatrix<> zero;
        zero.setZero(2,1);

        gsMatrix<> dd_ik_plus, dd_ik_minus, dd_ik_minus_deriv, dd_ik_plus_deriv;
        dd_ik_minus = -1/(gluingData.at(0).get_alpha_tilde().eval(zero.row(0))(0,0)) * (geo.jacobian(zero).col(1) +
            gluingData.at(0).get_beta_tilde().eval(zero.row(0))(0,0) * geo.jacobian(zero).col(0));

        dd_ik_plus = 1/(gluingData.at(1).get_alpha_tilde().eval(zero.row(0))(0,0)) * (geo.jacobian(zero).col(0) +
            gluingData.at(1).get_beta_tilde().eval(zero.row(0))(0,0) * geo.jacobian(zero).col(1));

        gsMatrix<> geo_deriv2_12(2,1), geo_deriv2_11(2,1), geo_deriv2_22(2,1);
        geo_deriv2_12.row(0) = geo.deriv2(zero).row(2);
        geo_deriv2_12.row(1) = geo.deriv2(zero).row(5);
        geo_deriv2_11.row(0) = geo.deriv2(zero).row(0);
        geo_deriv2_11.row(1) = geo.deriv2(zero).row(3);
        geo_deriv2_22.row(0) = geo.deriv2(zero).row(1);
        geo_deriv2_22.row(1) = geo.deriv2(zero).row(4);
        gsMatrix<> alpha_squared_u = gluingData.at(0).get_alpha_tilde().eval(zero.row(0))*gluingData.at(0).get_alpha_tilde().eval(zero.row(0));
        gsMatrix<> alpha_squared_v = gluingData.at(1).get_alpha_tilde().eval(zero.row(0))*gluingData.at(1).get_alpha_tilde().eval(zero.row(0));

        dd_ik_minus_deriv = -1/(alpha_squared_u(0,0)) * // N^2
            ((geo_deriv2_12 + (gluingData.at(0).get_beta_tilde().deriv(zero.row(0))(0,0) * geo.jacobian(zero).col(0) +
             gluingData.at(0).get_beta_tilde().eval(zero.row(0))(0,0) * geo_deriv2_11))*gluingData.at(0).get_alpha_tilde().eval(zero.row(0))(0,0) -
            (geo.jacobian(zero).col(1) + gluingData.at(0).get_beta_tilde().eval(zero.row(0))(0,0) * geo.jacobian(zero).col(0)) *
            gluingData.at(0).get_alpha_tilde().deriv(zero.row(0))(0,0));

        dd_ik_plus_deriv = 1/(alpha_squared_v(0,0)) *
            ((geo_deriv2_12 + (gluingData.at(1).get_beta_tilde().deriv(zero.row(0))(0,0) * geo.jacobian(zero).col(1) +
            gluingData.at(1).get_beta_tilde().eval(zero.row(0))(0,0) * geo_deriv2_22))*gluingData.at(1).get_alpha_tilde().eval(zero.row(0))(0,0) -
            (geo.jacobian(zero).col(0) + gluingData.at(1).get_beta_tilde().eval(zero.row(0))(0,0) * geo.jacobian(zero).col(1)) *
            gluingData.at(1).get_alpha_tilde().deriv(zero.row(0))(0,0));

        // Comupute d_(0,0)^(i_k), d_(1,0)^(i_k), d_(0,1)^(i_k), d_(1,1)^(i_k) ; i_k == 2
        std::vector<gsMatrix<>> d_ik;
        d_ik.push_back(Phi.col(0));
        d_ik.push_back(Phi.block(0,1,6,2) * geo.jacobian(zero).col(0) ); // deriv into u
        d_ik.push_back(Phi.block(0,1,6,2) * geo.jacobian(zero).col(1) ); // deriv into v
        d_ik.push_back((geo.jacobian(zero)(0,0) * Phi.col(3) + geo.jacobian(zero)(1,0) * Phi.col(4))*geo.jacobian(zero)(0,1) +
                       (geo.jacobian(zero)(0,0) * Phi.col(4) + geo.jacobian(zero)(1,0) * Phi.col(5))*geo.jacobian(zero)(1,1) +
                        Phi.block(0,1,6,1) * geo.deriv2(zero).row(2) +
                        Phi.block(0,2,6,1) * geo.deriv2(zero).row(5)); // Hessian

        // Compute d_(*,*)^(il,ik)
        std::vector<gsMatrix<>> d_ilik_minus, d_ilik_plus;
        d_ilik_minus.push_back(Phi.col(0));
        d_ilik_minus.push_back(Phi.block(0,1,6,2) * geo.jacobian(zero).col(0));
        d_ilik_minus.push_back((geo.jacobian(zero)(0,0) * Phi.col(3) + geo.jacobian(zero)(1,0) * Phi.col(4))*geo.jacobian(zero)(0,0) +
            (geo.jacobian(zero)(0,0) * Phi.col(4) + geo.jacobian(zero)(1,0) * Phi.col(5))*geo.jacobian(zero)(1,0) +
            Phi.block(0,1,6,1) * geo.deriv2(zero).row(0) +
            Phi.block(0,2,6,1) * geo.deriv2(zero).row(3));
        d_ilik_minus.push_back(Phi.block(0,1,6,2) * dd_ik_minus);
        d_ilik_minus.push_back((geo.jacobian(zero)(0,0) * Phi.col(3) + geo.jacobian(zero)(1,0) * Phi.col(4))*dd_ik_minus(0,0) +
            (geo.jacobian(zero)(0,0) * Phi.col(4) + geo.jacobian(zero)(1,0) * Phi.col(5))*dd_ik_minus(1,0) +
            Phi.block(0,1,6,1) * dd_ik_minus_deriv.row(0) +
            Phi.block(0,2,6,1) * dd_ik_minus_deriv.row(1));

        d_ilik_plus.push_back(Phi.col(0));
        d_ilik_plus.push_back(Phi.block(0,1,6,2) * geo.jacobian(zero).col(1));
        d_ilik_plus.push_back((geo.jacobian(zero)(0,1) * Phi.col(3) + geo.jacobian(zero)(1,1) * Phi.col(4))*geo.jacobian(zero)(0,1) +
            (geo.jacobian(zero)(0,1) * Phi.col(4) + geo.jacobian(zero)(1,1) * Phi.col(5))*geo.jacobian(zero)(1,1) +
            Phi.block(0,1,6,1) * geo.deriv2(zero).row(1) +
            Phi.block(0,2,6,1) * geo.deriv2(zero).row(4));
        d_ilik_plus.push_back(Phi.block(0,1,6,2) * dd_ik_plus);
        d_ilik_plus.push_back((geo.jacobian(zero)(0,1) * Phi.col(3) + geo.jacobian(zero)(1,1) * Phi.col(4))*dd_ik_plus(0,0) +
            (geo.jacobian(zero)(0,1) * Phi.col(4) + geo.jacobian(zero)(1,1) * Phi.col(5))*dd_ik_plus(1,0) +
            Phi.block(0,1,6,1) * dd_ik_plus_deriv.row(0) +
            Phi.block(0,2,6,1) * dd_ik_plus_deriv.row(1));

        for (index_t i = 0; i < 6; i++)
        {
            rhsVals.at(i) = d_ilik_minus.at(0)(i,0) * (c_0_plus.at(0).cwiseProduct(c_0.at(1)) -
                gluingData.at(0).get_beta_tilde().eval(md.points.row(0)).cwiseProduct(c_0_plus_deriv.at(0).cwiseProduct(c_1.at(1)))) +
                d_ilik_minus.at(1)(i,0) * (c_1_plus.at(0).cwiseProduct(c_0.at(1)) -
                gluingData.at(0).get_beta_tilde().eval(md.points.row(0)).cwiseProduct(c_1_plus_deriv.at(0).cwiseProduct(c_1.at(1)))) +
                d_ilik_minus.at(2)(i,0) * (c_2_plus.at(0).cwiseProduct(c_0.at(1)) -
                gluingData.at(0).get_beta_tilde().eval(md.points.row(0)).cwiseProduct(c_2_plus_deriv.at(0).cwiseProduct(c_1.at(1)))) -
                d_ilik_minus.at(3)(i,0) * gluingData.at(0).get_alpha_tilde().eval(md.points.row(0)).cwiseProduct(c_0_minus.at(0).cwiseProduct(c_1.at(1))) -
                d_ilik_minus.at(4)(i,0) * gluingData.at(0).get_alpha_tilde().eval(md.points.row(0)).cwiseProduct(c_1_minus.at(0).cwiseProduct(c_1.at(1))); // f*_(ik-1,ik)

            rhsVals.at(i) += d_ilik_plus.at(0)(i,0) * (c_0_plus.at(1).cwiseProduct(c_0.at(0)) -
                gluingData.at(1).get_beta_tilde().eval(md.points.row(1)).cwiseProduct(c_0_plus_deriv.at(1).cwiseProduct(c_1.at(0)))) +
                d_ilik_plus.at(1)(i,0) * (c_1_plus.at(1).cwiseProduct(c_0.at(0)) -
                gluingData.at(1).get_beta_tilde().eval(md.points.row(1)).cwiseProduct(c_1_plus_deriv.at(1).cwiseProduct(c_1.at(0)))) +
                d_ilik_plus.at(2)(i,0) * (c_2_plus.at(1).cwiseProduct(c_0.at(0)) -
                gluingData.at(1).get_beta_tilde().eval(md.points.row(1)).cwiseProduct(c_2_plus_deriv.at(1).cwiseProduct(c_1.at(0)))) +
                d_ilik_plus.at(3)(i,0) * gluingData.at(1).get_alpha_tilde().eval(md.points.row(1)).cwiseProduct(c_0_minus.at(1).cwiseProduct(c_1.at(0))) +
                d_ilik_plus.at(4)(i,0) * gluingData.at(1).get_alpha_tilde().eval(md.points.row(1)).cwiseProduct(c_1_minus.at(1).cwiseProduct(c_1.at(0))); // f*_(ik+1,ik)

            rhsVals.at(i) -= d_ik.at(0)(i,0) * c_0.at(0).cwiseProduct(c_0.at(1)) + d_ik.at(2)(i,0) * c_0.at(0).cwiseProduct(c_1.at(1)) +
                d_ik.at(1)(i,0) * c_1.at(0).cwiseProduct(c_0.at(1)) + d_ik.at(3)(i,0) * c_1.at(0).cwiseProduct(c_1.at(1)); // f*_(ik)

            localMat.at(i).setZero(numActive, numActive);
            localRhs.at(i).setZero(numActive, rhsVals.at(i).rows());//multiple right-hand sides

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
                basisData * quWeights.asDiagonal() *
                    md.measures.asDiagonal() * basisData.transpose();

            for (index_t k = 0; k < quWeights.rows(); ++k) // loop over quadrature nodes
            {
                // Multiply weight by the geometry measure
                const T weight = quWeights[k] * md.measure(k);

                localRhs.at(i).noalias() += weight * (basisVals.col(k) * rhsVals.at(i).col(k).transpose());
            }
        }
    }

    inline void localToGlobal(const int patchIndex,
                              const std::vector<gsMatrix<T> >    & eliminatedDofs,
                              std::vector< gsSparseSystem<T> >     & system)
    {
        gsMatrix<unsigned> actives_temp;
        for (unsigned i = 0; i < system.size(); i++) // 6
        {
            // Map patch-local DoFs to global DoFs
            system.at(i).mapColIndices(actives, patchIndex, actives_temp);
            // Add contributions to the system matrix and right-hand side
            system.at(i).push(localMat.at(i), localRhs.at(i), actives_temp, eliminatedDofs[0], 0, 0);
        }
    }

protected:
    gsMatrix<unsigned> actives;
    gsMatrix<T> basisData;
    index_t numActive;

protected:
    // Local values of the right hand side
    std::vector< gsMatrix<T> >  rhsVals;

protected:
    // Local matrices
    std::vector< gsMatrix<T> > localMat;
    std::vector< gsMatrix<T> > localRhs;

    gsMapData<T> md;

}; // class gsVisitorG1BasisVertex
} // namespace gismo