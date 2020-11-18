//
// Created by afarahat on 4/16/20.
//


#pragma once

#include "gsG1Basis/gsG1ASGluingData.h"


namespace gismo
{
template <class T>
class gsG1ASVisitorBasisVertex
{
public:

    gsG1ASVisitorBasisVertex()
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
                         gsMatrix<>  & gluingData,
                         std::vector<bool> & isBoundary,
                         gsMatrix<>  & Phi,
                         gsG1OptionList g1OptionList)
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

//        if (g1OptionList.getInt("gluingData") == gluingData::global)
//        {
//            alpha.push_back(gluingData[0].get_alpha_tilde().eval(md.points.row(0))); // u
//            alpha.push_back(gluingData[1].get_alpha_tilde().eval(md.points.row(1))); // v
//            alpha_0.push_back(gluingData[0].get_alpha_tilde().eval(zero.row(0))); // u
//            alpha_0.push_back(gluingData[1].get_alpha_tilde().eval(zero.row(0))); // v
//            alpha_deriv.push_back(gluingData[0].get_alpha_tilde().deriv(zero.row(0))); // u
//            alpha_deriv.push_back(gluingData[1].get_alpha_tilde().deriv(zero.row(0))); // v
//
//            beta.push_back(gluingData[0].get_beta_tilde().eval(md.points.row(0))); // u
//            beta.push_back(gluingData[1].get_beta_tilde().eval(md.points.row(1))); // v
//            beta_0.push_back(gluingData[0].get_beta_tilde().eval(zero.row(0))); // u
//            beta_0.push_back(gluingData[1].get_beta_tilde().eval(zero.row(0))); // v
//            beta_deriv.push_back(gluingData[0].get_beta_tilde().deriv(zero.row(0))); // u
//            beta_deriv.push_back(gluingData[1].get_beta_tilde().deriv(zero.row(0))); // v
//
//        }
//        else if (g1OptionList.getInt("gluingData") == gluingData::local)
//        {
//            alpha.push_back(gluingData[0].get_local_alpha_tilde(0).eval(md.points.row(0))); // u
//            alpha.push_back(gluingData[1].get_local_alpha_tilde(0).eval(md.points.row(1))); // v
//            alpha_0.push_back(gluingData[0].get_local_alpha_tilde(0).eval(zero.row(0))); // u
//            alpha_0.push_back(gluingData[1].get_local_alpha_tilde(0).eval(zero.row(0))); // v
//            alpha_deriv.push_back(gluingData[0].get_local_alpha_tilde(0).deriv(zero.row(0))); // u
//            alpha_deriv.push_back(gluingData[1].get_local_alpha_tilde(0).deriv(zero.row(0))); // v
//
//            beta.push_back(gluingData[0].get_local_beta_tilde(0).eval(md.points.row(0))); // u
//            beta.push_back(gluingData[1].get_local_beta_tilde(0).eval(md.points.row(1))); // v
//            beta_0.push_back(gluingData[0].get_local_beta_tilde(0).eval(zero.row(0))); // u
//            beta_0.push_back(gluingData[1].get_local_beta_tilde(0).eval(zero.row(0))); // v
//            beta_deriv.push_back(gluingData[0].get_local_beta_tilde(0).deriv(zero.row(0))); // u
//            beta_deriv.push_back(gluingData[1].get_local_beta_tilde(0).deriv(zero.row(0))); // v
//        }

        std::vector<gsMatrix<>> alpha, beta, alpha_0, beta_0, alpha_deriv, beta_deriv;

        // Point zero
        gsMatrix<> zero;
        zero.setZero(2,1);

        gsMatrix<> zeros(1, md.points.cols());
        zeros.setZero();

        // Point One
        gsMatrix<> one;
        one.setOnes(2,1);

        gsMatrix<> ones(1, md.points.cols());
        ones.setOnes();

        alpha.push_back( gluingData(0, 0) * ( ones - md.points.row(0) ) + gluingData(1, 0) * md.points.row(0) ); // u
        alpha.push_back( gluingData(0, 1) * ( ones - md.points.row(1) ) + gluingData(1, 1) * md.points.row(1) ); // v


        alpha_0.push_back( gluingData(0, 0) * ( one.row(0) - zero.row(0) ) + gluingData(1, 0) * zero.row(0) ); // u
        alpha_0.push_back( gluingData(0, 1) * ( one.row(0) - zero.row(0) ) + gluingData(1, 1) * zero.row(0) ); // v
        alpha_deriv.push_back( ( gluingData(1, 0) - gluingData(0, 0) ) * ones.col(0) ); // u
        alpha_deriv.push_back( ( gluingData(1, 1) - gluingData(0, 1) ) * ones.col(0) ); // v

        beta.push_back( gluingData(2, 0) * ( ones - md.points.row(0) ) + gluingData(3, 0) * md.points.row(0) ); // u
        beta.push_back( gluingData(2, 1) * ( ones - md.points.row(1) ) + gluingData(3, 1) * md.points.row(1) ); // v
        beta_0.push_back( gluingData(2, 0) * ( one.row(0) - zero.row(0) ) + gluingData(3, 0) * zero.row(0) ); // u
        beta_0.push_back( gluingData(2, 1) * ( one.row(0) - zero.row(0) ) + gluingData(3, 1) * zero.row(0) ); // v
        beta_deriv.push_back( ( gluingData(3, 0) - gluingData(2, 0) ) * ones.col(0) ); // u
        beta_deriv.push_back( ( gluingData(3, 1) - gluingData(2, 1) ) * ones.col(0) ); // v

        // Compute dd^^(i_k) and dd^^(i_k-1)
        gsMatrix<> dd_ik_plus, dd_ik_minus;
        gsMatrix<> dd_ik_minus_deriv, dd_ik_plus_deriv;

        dd_ik_minus = ( -1 / alpha_0[0](0,0) ) * ( geo.jacobian(zero).col(1) +
                                                            beta_0[0](0,0) * geo.jacobian(zero).col(0) );

        dd_ik_plus = ( 1 / alpha_0[1](0,0) ) * ( geo.jacobian(zero).col(0) +
                                                          beta_0[1](0,0) * geo.jacobian(zero).col(1) );

        gsMatrix<> geo_deriv2_12(geo.targetDim(), 1), geo_deriv2_11(geo.targetDim(), 1), geo_deriv2_22(geo.targetDim(), 1);

                   geo_deriv2_12.row(0) = geo.deriv2(zero).row(2);
                   geo_deriv2_12.row(1) = geo.deriv2(zero).row(5);

                   geo_deriv2_11.row(0) = geo.deriv2(zero).row(0);
                   geo_deriv2_11.row(1) = geo.deriv2(zero).row(3);

                   geo_deriv2_22.row(0) = geo.deriv2(zero).row(1);
                   geo_deriv2_22.row(1) = geo.deriv2(zero).row(4);

        if(geo.parDim() + 1 == geo.targetDim()) // In the surface case the dimension of the second derivative vector is 9x1
        {
            geo_deriv2_12.row(2) = geo.deriv2(zero).row(8);

            geo_deriv2_11.row(2) = geo.deriv2(zero).row(6);

            geo_deriv2_22.row(2) = geo.deriv2(zero).row(7);
        }

        gsMatrix<> alpha_squared_u = alpha_0[0]*alpha_0[0];
        gsMatrix<> alpha_squared_v = alpha_0[1]*alpha_0[1];

        dd_ik_minus_deriv = -1/(alpha_squared_u(0,0)) * // N^2
            ( ( geo_deriv2_12 + (beta_deriv[0](0,0) * geo.jacobian(zero).col(0) + beta_0[0](0,0) * geo_deriv2_11) ) * alpha_0[0](0,0) -
              ( geo.jacobian(zero).col(1) + beta_0[0](0,0) * geo.jacobian(zero).col(0) ) * alpha_deriv[0](0,0) );


        dd_ik_plus_deriv = 1/(alpha_squared_v(0,0)) *
            ( ( geo_deriv2_12 + (beta_deriv[1](0,0) * geo.jacobian(zero).col(1) + beta_0[1](0,0) * geo_deriv2_22) ) * alpha_0[1](0,0) -
              ( geo.jacobian(zero).col(0) + beta_0[1](0,0) * geo.jacobian(zero).col(1) ) * alpha_deriv[1](0,0) );

        //if (isBoundary[0] == false)
        //    gsInfo << dd_ik_minus_deriv << "\n";
        //if (isBoundary[1] == false)
        //    gsInfo << dd_ik_plus_deriv << "\n";

        // Comupute d_(0,0)^(i_k), d_(1,0)^(i_k), d_(0,1)^(i_k), d_(1,1)^(i_k) ; i_k == 2
        std::vector<gsMatrix<>> d_ik;
        d_ik.push_back(Phi.row(0).transpose());

        d_ik.push_back(Phi.block(1, 0, geo.targetDim(), 6).transpose() * geo.jacobian(zero).col(0) ); // deriv into u


//        gsInfo << "d_ik: " << d_ik.back()  << "\n";
//        gsInfo << "======================================= \n";
//        gsInfo << "geo.jacobian(zero).col(0): " << geo.jacobian(zero).col(0)  << "\n";
//        gsInfo << "======================================= \n";


        d_ik.push_back(Phi.block(1, 0, geo.targetDim(), 6).transpose() * geo.jacobian(zero).col(1) ); // deriv into v

        // Hessian
        if(geo.parDim() + 1 == geo.targetDim()) // In the surface case the dimension of the second derivative vector is 9x1
        {
            d_ik.push_back( (geo.jacobian(zero)(0,0) * Phi.row(4).transpose() + geo.jacobian(zero)(1,0) * Phi.row(7).transpose() + geo.jacobian(zero)(2,0) * Phi.row(10).transpose()) * geo.jacobian(zero)(0,1) +
                            (geo.jacobian(zero)(0,0) * Phi.row(5).transpose() + geo.jacobian(zero)(1,0) * Phi.row(8).transpose() + geo.jacobian(zero)(2,0) * Phi.row(11).transpose()) * geo.jacobian(zero)(1,1) +
                            (geo.jacobian(zero)(0,0) * Phi.row(6).transpose() + geo.jacobian(zero)(1,0) * Phi.row(9).transpose() + geo.jacobian(zero)(2,0) * Phi.row(12).transpose()) * geo.jacobian(zero)(2,1) +
                            Phi.block(1, 0, 1, 6).transpose() * geo.deriv2(zero).row(2) +
                            Phi.block(2, 0, 1, 6).transpose() * geo.deriv2(zero).row(5) +
                            Phi.block(3, 0, 1, 6).transpose() * geo.deriv2(zero).row(8) );

        }
        else
        {
            d_ik.push_back( (geo.jacobian(zero)(0,0) * Phi.col(3) + geo.jacobian(zero)(1,0) * Phi.col(4)) * geo.jacobian(zero)(0,1) +
                            (geo.jacobian(zero)(0,0) * Phi.col(4) + geo.jacobian(zero)(1,0) * Phi.col(5)) * geo.jacobian(zero)(1,1) +
                            Phi.block(0,1, 6,1) * geo.deriv2(zero).row(2) +
                            Phi.block(0,2, 6,1) * geo.deriv2(zero).row(5)); // Hessian
        }

//        gsInfo << "d_ik: " << d_ik.back() << " : " << Phi << "\n";
        // Compute d_(*,*)^(il,ik)
        std::vector<gsMatrix<>> d_ilik_minus, d_ilik_plus;

//      d_(*,*)^(ik-1,ik)
        d_ilik_minus.push_back(Phi.row(0).transpose());

        d_ilik_minus.push_back(Phi.block(1, 0, geo.targetDim(), 6).transpose() * geo.jacobian(zero).col(0));

        if(geo.parDim() + 1 == geo.targetDim()) // In the surface case the dimension of the second derivative vector is 9x1
        {
            d_ilik_minus.push_back( (geo.jacobian(zero)(0,0) * Phi.row(4).transpose() + geo.jacobian(zero)(1,0) * Phi.row(7).transpose() + geo.jacobian(zero)(2,0) * Phi.row(10).transpose()) * geo.jacobian(zero)(0,0) +
                                    (geo.jacobian(zero)(0,0) * Phi.row(5).transpose() + geo.jacobian(zero)(1,0) * Phi.row(8).transpose() + geo.jacobian(zero)(2,0) * Phi.row(11).transpose()) * geo.jacobian(zero)(1,0) +
                                    (geo.jacobian(zero)(0,0) * Phi.row(6).transpose() + geo.jacobian(zero)(1,0) * Phi.row(9).transpose() + geo.jacobian(zero)(2,0) * Phi.row(12).transpose()) * geo.jacobian(zero)(2,0) +
                                    Phi.block(1, 0, 1, 6).transpose() * geo.deriv2(zero).row(0) +
                                    Phi.block(2, 0, 1, 6).transpose() * geo.deriv2(zero).row(3) +
                                    Phi.block(3, 0, 1, 6).transpose() * geo.deriv2(zero).row(6) );
        }
        else
        {
            d_ilik_minus.push_back( (geo.jacobian(zero)(0,0) * Phi.col(3) + geo.jacobian(zero)(1,0) * Phi.col(4))*geo.jacobian(zero)(0,0) +
                                    (geo.jacobian(zero)(0,0) * Phi.col(4) + geo.jacobian(zero)(1,0) * Phi.col(5))*geo.jacobian(zero)(1,0) +
                                    Phi.block(0,1,6,1) * geo.deriv2(zero).row(0) +
                                    Phi.block(0,2,6,1) * geo.deriv2(zero).row(3));
        }

        d_ilik_minus.push_back(Phi.block(1, 0, geo.targetDim(), 6).transpose() * dd_ik_minus);

        if(geo.parDim() + 1 == geo.targetDim()) // In the surface case the dimension of the second derivative vector is 9x1
        {
            d_ilik_minus.push_back( (geo.jacobian(zero)(0,0) * Phi.row(4).transpose() + geo.jacobian(zero)(1,0) * Phi.row(7).transpose() + geo.jacobian(zero)(2,0) * Phi.row(10).transpose()) * dd_ik_minus(0,0) +
                                    (geo.jacobian(zero)(0,0) * Phi.row(5).transpose() + geo.jacobian(zero)(1,0) * Phi.row(8).transpose() + geo.jacobian(zero)(2,0) * Phi.row(11).transpose()) * dd_ik_minus(1,0) +
                                    (geo.jacobian(zero)(0,0) * Phi.row(6).transpose() + geo.jacobian(zero)(1,0) * Phi.row(9).transpose() + geo.jacobian(zero)(2,0) * Phi.row(12).transpose()) * dd_ik_minus(2,0) +
                                    Phi.block(1, 0, 1, 6).transpose() * dd_ik_minus_deriv.row(0) +
                                    Phi.block(2, 0, 1, 6).transpose() * dd_ik_minus_deriv.row(1) +
                                    Phi.block(3, 0, 1, 6).transpose() * dd_ik_minus_deriv.row(2) );
        }
        else
        {
            d_ilik_minus.push_back( (geo.jacobian(zero)(0,0) * Phi.col(3) + geo.jacobian(zero)(1,0) * Phi.col(4)) * dd_ik_minus(0,0) +
                                    (geo.jacobian(zero)(0,0) * Phi.col(4) + geo.jacobian(zero)(1,0) * Phi.col(5)) * dd_ik_minus(1,0) +
                                    Phi.block(0,1,6,1) * dd_ik_minus_deriv.row(0) +
                                    Phi.block(0,2,6,1) * dd_ik_minus_deriv.row(1));
        }


//      d_(*,*)^(ik+1,ik)
        d_ilik_plus.push_back(Phi.row(0).transpose());

        d_ilik_plus.push_back(Phi.block(1, 0, geo.targetDim(), 6).transpose() * geo.jacobian(zero).col(1));

        if(geo.parDim() + 1 == geo.targetDim()) // In the surface case the dimension of the second derivative vector is 9x1
        {
            d_ilik_plus.push_back( (geo.jacobian(zero)(0,1) * Phi.row(4).transpose() + geo.jacobian(zero)(1,1) * Phi.row(7).transpose() + geo.jacobian(zero)(2,1) * Phi.row(10).transpose()) * geo.jacobian(zero)(0,1) +
                                   (geo.jacobian(zero)(0,1) * Phi.row(5).transpose() + geo.jacobian(zero)(1,1) * Phi.row(8).transpose() + geo.jacobian(zero)(2,1) * Phi.row(11).transpose()) * geo.jacobian(zero)(1,1) +
                                   (geo.jacobian(zero)(0,1) * Phi.row(6).transpose() + geo.jacobian(zero)(1,1) * Phi.row(9).transpose() + geo.jacobian(zero)(2,1) * Phi.row(12).transpose()) * geo.jacobian(zero)(2,1) +
                                    Phi.block(1, 0, 1, 6).transpose() * geo.deriv2(zero).row(1) +
                                    Phi.block(2, 0, 1, 6).transpose() * geo.deriv2(zero).row(4) +
                                    Phi.block(3, 0, 1, 6).transpose() * geo.deriv2(zero).row(7) );
        }
        else
        {
            d_ilik_plus.push_back( (geo.jacobian(zero)(0,1) * Phi.col(3) + geo.jacobian(zero)(1,1) * Phi.col(4)) * geo.jacobian(zero)(0,1) +
                                   (geo.jacobian(zero)(0,1) * Phi.col(4) + geo.jacobian(zero)(1,1) * Phi.col(5)) * geo.jacobian(zero)(1,1) +
                                   Phi.block(0,1,6,1) * geo.deriv2(zero).row(1) +
                                   Phi.block(0,2,6,1) * geo.deriv2(zero).row(4) );
        }

        d_ilik_plus.push_back(Phi.block(1, 0, geo.targetDim(), 6).transpose() * dd_ik_plus);

        if(geo.parDim() + 1 == geo.targetDim()) // In the surface case the dimension of the second derivative vector is 9x1
        {
            d_ilik_plus.push_back( (geo.jacobian(zero)(0,1) * Phi.row(4).transpose() + geo.jacobian(zero)(1,1) * Phi.row(7).transpose() + geo.jacobian(zero)(2,1) * Phi.row(10).transpose()) * dd_ik_plus(0,0) +
                                   (geo.jacobian(zero)(0,1) * Phi.row(5).transpose() + geo.jacobian(zero)(1,1) * Phi.row(8).transpose() + geo.jacobian(zero)(2,1) * Phi.row(11).transpose()) * dd_ik_plus(1,0) +
                                   (geo.jacobian(zero)(0,1) * Phi.row(6).transpose() + geo.jacobian(zero)(1,1) * Phi.row(9).transpose() + geo.jacobian(zero)(2,1) * Phi.row(12).transpose()) * dd_ik_plus(2,0) +
                                    Phi.block(1, 0, 1, 6).transpose() * dd_ik_plus_deriv.row(0) +
                                    Phi.block(2, 0, 1, 6).transpose() * dd_ik_plus_deriv.row(1) +
                                    Phi.block(3, 0, 1, 6).transpose() * dd_ik_plus_deriv.row(2) );
        }
        else
        {
            d_ilik_plus.push_back( (geo.jacobian(zero)(0,1) * Phi.col(3) + geo.jacobian(zero)(1,1) * Phi.col(4)) * dd_ik_plus(0,0) +
                                   (geo.jacobian(zero)(0,1) * Phi.col(4) + geo.jacobian(zero)(1,1) * Phi.col(5)) * dd_ik_plus(1,0) +
                                   Phi.block(0,1,6,1) * dd_ik_plus_deriv.row(0) +
                                   Phi.block(0,2,6,1) * dd_ik_plus_deriv.row(1) );
        }


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
            localMat.at(i).noalias() = basisData * quWeights.asDiagonal() * md.measures.asDiagonal() * basisData.transpose();

            for (index_t k = 0; k < quWeights.rows(); ++k) // loop over quadrature nodes
            {
                T weight = quWeights[k];
                gsMatrix<> Jk = md.jacobian(k);

                if( Jk.dim().second + 1 == Jk.dim().first)
                {
                    gsMatrix<> G = Jk.transpose() * Jk;
                    real_t detG = G.determinant();
                    weight *= sqrt(detG);
                }
                else
                {
                    weight *=  md.measure(k);
                }
                // Multiply weight by the geometry measure

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
