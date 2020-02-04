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
class gsVisitorNormL2
{
    typedef std::multimap<index_t, std::map<index_t, std::map<index_t, gsMultiPatch<real_t>>>> typedef_g1;

public:

    gsVisitorNormL2(const typedef_g1 & g1,
                    index_t p = 2):
                    m_G1Basis_mp(g1)
    {
        f2param = false;
        g1basis = false;
        g1basis_mp = true;
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
    void evaluate(gsGeometryEvaluator<T> & geoEval,
                  const gsFunction<T>    & _func1,
                  const gsFunction<T>    & _func2,
                  gsMatrix<T>            & quNodes)
    {
        // Evaluate first function
        _func1.eval_into(quNodes, f1vals);

        if (g1basis)
        {
            index_t n = m_G1Basis.at(geoEval.id()).nPatches();

            for (index_t i = 0; i < n; i++)
            {
                f1vals += m_G1Basis.at(geoEval.id()).patch(i).eval(quNodes);
            }
        }
        if (g1basis_mp)
        {
            std::multimap<index_t, std::map<index_t, std::map<index_t, gsMultiPatch<real_t>>>>::iterator i_face;
            for (i_face=m_G1Basis_mp.equal_range(geoEval.id()).first; i_face!=m_G1Basis_mp.equal_range(geoEval.id()).second; ++i_face)
            {
                std::map<index_t, std::map<index_t, gsMultiPatch<real_t>>>::iterator i_side;
                for (i_side = i_face->second.begin(); i_side != i_face->second.end(); ++i_side)
                {
                    for (std::map<index_t, gsMultiPatch<real_t>>::iterator i_mp = i_side->second.begin();
                        i_mp != i_side->second.end(); ++i_mp)
                    {
                        gsMultiPatch<real_t> mp_side = i_mp->second;
                        for (unsigned j = 0; j < mp_side.nPatches(); j++)
                        {
                            f1vals += mp_side.patch(j).eval(quNodes);
                        }
                    }
                }
            }
        }

        // Compute geometry related values
        geoEval.evaluateAt(quNodes);

        // Evaluate second function
        _func2.eval_into( f2param ? quNodes : geoEval.values() , f2vals);
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

    bool f2param;
    bool g1basis;
    bool g1basis_mp;

protected:

    std::vector< gsMultiPatch<>> m_G1Basis;
    typedef_g1 m_G1Basis_mp;
};






}







