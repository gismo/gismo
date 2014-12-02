#pragma once

#include <gsAssembler/gsNorm.h>

namespace gismo
{

template <class T>
class gsErrEstPoissonResidual : public gsNorm<T>
{
    friend class gsNorm<T>;

public:

    // f1 in gsNorm corresponds to discrete Solution
    // f2 in gsNorm corresponds to right-hand-side

    gsErrEstPoissonResidual(const gsField<T> & _discSolution,
             const gsFunction<T> & _rhsFunction,
             bool _rhsFunctionParam = false)
    : gsNorm<T>(_discSolution,_rhsFunction), m_f2param(_rhsFunctionParam)
    {

    }

public:

    T compute(bool storeElWise = true)
    {
        this->apply(*this,storeElWise);
        return this->m_value;
    }


protected:

    void initialize(const gsBasis<T> & basis,
                    gsQuadRule<T> & rule,
                    unsigned      & evFlags) // replace with geoEval ?
    {
        // Setup Quadrature
        m_parDim = basis.dim();

        GISMO_ASSERT(m_parDim == 2 || m_parDim == 3, "Called error estimator with dimension other than 2 or 3.");

        gsVector<index_t> numQuadNodes( m_parDim );
        for (unsigned i = 0; i < m_parDim; ++i)
            numQuadNodes[i] = basis.degree(i) + 1;

        // Setup Quadrature
        rule = gsGaussRule<T>(numQuadNodes);// harmless slicing occurs here

        // Set Geometry evaluation flags
        evFlags = NEED_MEASURE| NEED_VALUE| NEED_JACOBIAN;
    }

    // Evaluate on element.
    inline void evaluate(gsGeometryEvaluator<T> & geoEval,
                         const gsGeometry<T>    & discSolution,
                         const gsFunction<T>    & rhsFunction,
                         gsMatrix<T>            & quNodes)
    {
        // Evaluate first function
        discSolution.deriv2_into(quNodes, m_discSol2ndDer);
        discSolution.eval_into(quNodes, m_discSolVals);

        // Compute geometry related values
        geoEval.evaluateAt(quNodes);

        // Evaluate second function (defined of physical domain)
        rhsFunction.eval_into(geoEval.values(), m_rhsFctVals);
    }

    // assemble on element
    inline T compute(gsDomainIterator<T>    & element,
                     gsGeometryEvaluator<T> & geoEval,
                     gsVector<T> const      & quWeights)
    {
        T sum(0.0);

        for (index_t k = 0; k < quWeights.rows(); ++k) // loop over quadrature nodes
        {
            const T weight = quWeights[k] * geoEval.measure(k);

            const typename gsMatrix<T>::constColumns J = geoEval.jacobian(k);
            gsMatrix<T> sol_der2 = m_discSol2ndDer.col(k);

            // Compute the APPROXIMATION of the
            // transformation of the
            // Laplacian to the physical domain.
            // Note that the term involing the second
            // derivative of the inverse geometry mapping is
            // neglected!

            T sol_Lap(0.0);

            if( m_parDim == 2 )
            {
                gsMatrix<T> Jinv = J.inverse();

                for( int i=0; i < 2; i ++ )
                {
                    sol_Lap += sol_der2(0,0) * Jinv(0,i) * Jinv(0,i) \
                        + sol_der2(2,0) * Jinv(0,i) * Jinv(1,i) \
                        + sol_der2(2,0) * Jinv(1,i) * Jinv(0,i) \
                        + sol_der2(1,0) * Jinv(1,i) * Jinv(1,i);
                }
            }
            else if( m_parDim == 3 )
            {
                gsMatrix<T> Jinv = J.inverse();

                for( int i=0; i < 3; i ++ )
                {
                    sol_Lap += \
                          sol_der2(0,0) * Jinv(0,i) * Jinv(0,i) \
                        + sol_der2(3,0) * Jinv(0,i) * Jinv(1,i) \
                        + sol_der2(4,0) * Jinv(0,i) * Jinv(2,i) \
                        + sol_der2(3,0) * Jinv(1,i) * Jinv(0,i) \
                        + sol_der2(1,0) * Jinv(1,i) * Jinv(1,i) \
                        + sol_der2(5,0) * Jinv(1,i) * Jinv(2,i) \
                        + sol_der2(4,0) * Jinv(2,i) * Jinv(0,i) \
                        + sol_der2(5,0) * Jinv(2,i) * Jinv(1,i) \
                        + sol_der2(2,0) * Jinv(2,i) * Jinv(2,i);
                }
            }

            sum += weight * ( sol_Lap + m_rhsFctVals(0,k) ) \
                    * ( sol_Lap + m_rhsFctVals(0,k) );

        }

        return sum;
    }

private:

    gsMatrix<T> m_discSol2ndDer, m_discSolVals, m_rhsFctVals;

    unsigned m_parDim;

    bool m_f2param;
};


} // namespace gismo

