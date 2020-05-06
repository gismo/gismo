/** @file gsNormL2.h

    @brief Computes the L2 norm, needs for the parallel computing.

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
class gsVisitorGlobalGDNorm
{

public:

    gsVisitorGlobalGDNorm(index_t p = 2)
    {
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
        evFlags = NEED_MEASURE| NEED_VALUE;
    }

    // Evaluate on element.
    void evaluate(const gsMultiPatch<T> & mp,
                  gsBSplineBasis<T>       & basis,
                  const gsMatrix<T> & sol,
                  gsMatrix<T>            & quNodes)
    {
        gsMatrix<unsigned> actives;

        gsMatrix<unsigned> actives_temp;
        basis.active_into(quNodes.col(0), actives_temp);

        actives.resize(actives_temp.rows()*4,1);
        for (index_t j = 0; j < 4; j++)
            for (index_t i = 0; i < actives_temp.rows(); i++)
                actives(i+j*actives_temp.rows(),0) = actives_temp(i,0) + j*basis.size();

        // ++++++++++++++++++++++++++++++++
        // Create alpha^S (x_old) and beta^S (x_old)
        // ++++++++++++++++++++++++++++++++
        gsGeometry<>::uPtr tilde_temp;

        tilde_temp = basis.makeGeometry(sol.block(0,0,basis.size(),1));
        gsBSpline<T> alpha_L = dynamic_cast<gsBSpline<T> &> (*tilde_temp);

        tilde_temp = basis.makeGeometry(sol.block(basis.size(),0,basis.size(),1));
        gsBSpline<T> alpha_R = dynamic_cast<gsBSpline<T> &> (*tilde_temp);

        tilde_temp = basis.makeGeometry(sol.block(2*basis.size(),0,basis.size(),1));
        gsBSpline<T> beta_L = dynamic_cast<gsBSpline<T> &> (*tilde_temp);

        tilde_temp = basis.makeGeometry(sol.block(3*basis.size(),0,basis.size(),1));
        gsBSpline<T> beta_R = dynamic_cast<gsBSpline<T> &> (*tilde_temp);


        gsMatrix<> a_L, a_R, b_L, b_R;
        alpha_L.eval_into(quNodes,a_L);
        alpha_R.eval_into(quNodes,a_R);

        beta_L.eval_into(quNodes,b_L);
        beta_R.eval_into(quNodes,b_R);
        // ++++++++++++++++++++++++++++++++

        // ++++++++++++++++++++++++++++++++
        // Compute the mapping
        // ++++++++++++++++++++++++++++++++
        gsGeometry<> & FR = mp.patch(0);
        gsGeometry<> & FL = mp.patch(1);

        gsMatrix<> pointV(FR.parDim(), quNodes.cols());
        pointV.setZero();
        pointV.row(1) = quNodes;

        gsMatrix<> pointU(FL.parDim(), quNodes.cols());
        pointU.setZero();
        pointU.row(0) = quNodes;

        gsMatrix<> DuFR(FR.targetDim(), quNodes.cols());
        gsMatrix<> DvFR(FR.targetDim(), quNodes.cols());
        gsMatrix<> DvFL(FR.targetDim(), quNodes.cols());

        gsMatrix<> DuFR_DuFR(1, quNodes.cols());
        gsMatrix<> DvFL_DvFL(1, quNodes.cols());
        gsMatrix<> DvFR_DvFR(1, quNodes.cols());


        for(index_t i = 0; i < quNodes.cols(); i++)
        {
            DuFR.col(i) = FR.jacobian(pointV.col(i)).col(0);
            DvFR.col(i) = FR.jacobian(pointV.col(i)).col(1); // Same as DuFL
            DvFL.col(i) = FL.jacobian(pointU.col(i)).col(1);

            // Set scalar product of the jacobian vectors of the geometric mapping
            DuFR_DuFR.col(i) = DuFR.col(i).transpose() * DuFR.col(i);
            DvFL_DvFL.col(i) = DvFL.col(i).transpose() * DvFL.col(i);
            DvFR_DvFR.col(i) = DvFR.col(i).transpose() * DvFR.col(i);

        }
        // ++++++++++++++++++++++++++++++++

        // ++++++++++++++++++++++++++++++++
        // Compute rhs
        // ++++++++++++++++++++++++++++++++
        f1vals.setZero(2,quNodes.cols());
        f1vals.row(0) = a_R.cwiseProduct(DvFL.row(0)) + a_L.cwiseProduct(DuFR.row(0)) + (a_R.cwiseProduct(b_L) + a_L.cwiseProduct(b_R)).cwiseProduct(DvFR.row(0));
        f1vals.row(1) = a_R.cwiseProduct(DvFL.row(1)) + a_L.cwiseProduct(DuFR.row(1)) + (a_R.cwiseProduct(b_L) + a_L.cwiseProduct(b_R)).cwiseProduct(DvFR.row(1));
        // ++++++++++++++++++++++++++++++++

    }

    // assemble on element
    T compute(gsDomainIterator<T>    & ,
              gsVector<T> const      & quWeights,
              T & accumulated)
    {
        T sum(0.0);
        for (index_t k = 0; k < quWeights.rows(); ++k) // loop over quadrature nodes
        {
            const T weight = quWeights[k];
            switch (m_p)
            {
                case 0: // infinity norm
                    // .template lpNorm<Eigen::Infinity>();
                    sum = (f1vals).array().abs().maxCoeff();
                    accumulated = math::max(accumulated, sum);
                    return sum;
                    break;
                case 1:
                    sum += weight * ( f1vals.col(k) ).template lpNorm<1>();
                    break;
                case 2:
                    sum += weight * ( f1vals.col(k)).squaredNorm();
                    break;
                default:
                    sum += weight * ( f1vals.col(k) ).array().abs().pow(m_p).sum();
                    //sum += weight * ( f1vals.col(k) - f2vals.col(k) ).template lpNorm<p>().squared();
            }
        }
        accumulated += sum;
        return sum;
    }

private:

    index_t m_p;

    gsMatrix<T> f1vals;


};






}







