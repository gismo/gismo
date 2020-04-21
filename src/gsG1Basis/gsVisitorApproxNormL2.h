/** @file gsVisitorApproxNormL2.h

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
class gsVisitorApproxNormL2
{

public:

    gsVisitorApproxNormL2(index_t p = 2)
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
    void evaluate(const gsMultiPatch<> * multiPatch,
                  gsGeometryEvaluator<T> & geoEval,
                  const gsBasis<T>       & basis,
                  const gsMultiPatch<T> * approx,
                  gsMatrix<T>            & quNodes,
                  gsOptionList optionList)
    {
        gsMatrix<unsigned> actives;
        gsMatrix<T> basisData;

        basis.active_into(quNodes.col(0), actives);

        // Evaluate basis functions on element
        basis.eval_into(quNodes,basisData);

        f1vals.setZero(1,actives.rows());
        f1vals += approx->patch(0).eval(quNodes);

        // Compute geometry related values
        geoEval.evaluateAt(quNodes);

        // ++++++++++++++++++++++++++++++++
        // Compute alpha^S and beta^S exact
        // ++++++++++++++++++++++++++++++++
        // alpha^S
        gsMatrix<> uv, ev;

        uv.setZero(2,quNodes.cols());
        uv.bottomRows(1) = quNodes; // v

        const index_t d = geoEval.parDim();
        gsVector<> D0(d);

        // ======== Determine bar{beta}^L ========
        const gsGeometry<> & P0 = multiPatch->patch(0); // iFace.second().patch = 0

        for(index_t i = 0; i < uv.cols(); i++)
        {
            P0.jacobian_into(uv.col(i),ev);
            D0 = ev.col(1);
            real_t D1 = 1/ D0.norm();
            uv(0,i) = - D1 * D1 * ev.col(1).transpose() * ev.col(0);

        }

        // alpha:
        /*
        uv.setZero(2,quNodes.cols());
        uv.bottomRows(1) = quNodes; // v
        for (index_t i = 0; i < uv.cols(); i++)
        {
            P0.jacobian_into(uv.col(i), ev);
            uv(0, i) = ev.determinant();
        }
        */


        index_t m_r = optionList.getInt("regularity");

        gsBSplineBasis<> basis_edge = dynamic_cast<gsBSplineBasis<> &>(multiPatch->basis(0).component(1)); // 0 -> v, 1 -> u
        index_t m_p = basis_edge.maxDegree();

        // first,last,interior,mult_ends,mult_interior
        gsKnotVector<T> kv_plus(0,1,0,m_p+1,m_p-1-m_r); // p,r+1 //-1 bc r+1
        gsBSplineBasis<> basis_plus(kv_plus);
        for (size_t i = m_p+1; i < basis_edge.knots().size() - (m_p+1); i = i+(m_p-m_r))
            basis_plus.insertKnot(basis_edge.knot(i),m_p-1-m_r);


        gsKnotVector<T> kv_minus(0,1,0,m_p+1-1,m_p-1-m_r); // p-1,r //-1 bc p-1
        gsBSplineBasis<> basis_minus(kv_minus);
        for (size_t i = m_p+1; i < basis_edge.knots().size() - (m_p+1); i = i+(m_p-m_r))
            basis_minus.insertKnot(basis_edge.knot(i),m_p-1-m_r);


        f2vals.setZero(1,actives.rows());

        gsMatrix<> der_b_plus, b_plus, b_minus;
        basis_plus.evalSingle_into(optionList.getInt("basisID"),quNodes,b_plus);
        basis_plus.derivSingle_into(optionList.getInt("basisID"),quNodes,der_b_plus);
        basis_minus.evalSingle_into(optionList.getInt("basisID"),quNodes,b_minus);

        //f2vals = b_plus + uv.row(0).cwiseProduct(der_b_plus);
        f2vals = b_plus + uv.row(0).cwiseProduct(der_b_plus);
    }

    // assemble on element
    T compute(gsDomainIterator<T>    & ,
              gsGeometryEvaluator<T> & geoEval,
              gsVector<T> const      & quWeights,
              T & accumulated)
    {
        T sum(0.0);
        for (index_t k = 0; k < quWeights.rows(); ++k) // loop over quadrature nodes
        {
            const T weight = quWeights[k] * geoEval.measure(k);
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


};






}







