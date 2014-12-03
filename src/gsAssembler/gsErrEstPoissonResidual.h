/** @file gsErrEstPoissonResidual.h

    @brief Residual-type error estimator for the Poisson problem.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): S. Kleiss
*/

#pragma once

#include <gsAssembler/gsNorm.h>

namespace gismo
{


/** \brief Provides a Residual-type and element-wise error estimator
 * for the Poisson problem.
 *
 * Let the Poission problem on the domain \f$ \Omega \f$ be given by
 * \f[ -\Delta u = f,\quad u = u_D \mathrm{\ on\ } \partial \Omega. \f]
 * The error estimate \f$\eta\f$ for a computed discrete solution
 * \f$ u_h \f$ is then given by
 * \f[ \eta^2 = \sum_K \eta_K^2 = \sum_K \| \Delta u_h + f\|_{L_2(K)}^2, \f]
 * where \f$K\f$ denotes a element/cell of the mesh
 * (and, naturally, \f$ \sum_K\f$ is a sum over all cells).
 *
 * \warning Note that the terms regarding Neumann boundary conditions
 * and jumps of the normal derivative across element interfaces are
 * NEGLECTED in the current version (02.Dec.2014).
 */
template <class T>
class gsErrEstPoissonResidual : public gsNorm<T>
{
    friend class gsNorm<T>;

public:

    // f1 in gsNorm corresponds to discrete Solution
    // f2 in gsNorm corresponds to right-hand-side

    /** \brief Constructor
     * \param _discSolution Discrete solution
     * \param _rhsFunction Right-hand-side-/Source-function \f$ f\f$
     * of the Poisson problem.
     * \param _rhsFunctionParam Flag indicating whether the \em _rhsFunction
     * is parameterized (in this case, the evaluation points must be given
     * on the parameter domain
     */
    gsErrEstPoissonResidual(const gsField<T> & _discSolution,
             const gsFunction<T> & _rhsFunction,
             bool _rhsFunctionParam = false)
    : gsNorm<T>(_discSolution,_rhsFunction), m_f2param(_rhsFunctionParam)
    {

    }

public:

    /** \brief Computes the error estimate.
     *
     * Computes the residual-based error estimate \f$\eta\f$
     * (see class-documentation at the top).
     *
     *
     * \param storeEltWise Bool indicating whether the element-wise
     * errors should be stored also. If <em>storeEletWise = true</em>,
     * the gsVector of element-wise estimates \f$\eta_K\f$
     * can be obtained by
     * calling elementNorms().
     *
     * \returns The total estimated error \f$ \eta \f$.
     *
     */
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
        m_parDim = basis.dim();

        GISMO_ASSERT(m_parDim == 2 || m_parDim == 3, "Called error estimator with dimension other than 2 or 3.");

        // Setup Quadrature
        gsVector<index_t> numQuadNodes( m_parDim );
        for (unsigned i = 0; i < m_parDim; ++i)
            numQuadNodes[i] = basis.degree(i) + 1;

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
        // Evaluate discrete solution
        discSolution.deriv2_into(quNodes, m_discSol2ndDer);
        discSolution.eval_into(quNodes, m_discSolVals);

        // Compute geometry related values
        geoEval.evaluateAt(quNodes);

        // Evaluate right-hand-side function (defined of physical domain)
        rhsFunction.eval_into(geoEval.values(), m_rhsFctVals);
    }

    // assemble on element
    inline T compute(gsDomainIterator<T>    & element,
                     gsGeometryEvaluator<T> & geoEval,
                     gsVector<T> const      & quWeights)
    {
        T hh = T(1.0);
        for( index_t i=0; i < element.upperCorner().size(); i++)
            hh *= ( element.upperCorner()[i] - element.lowerCorner()[i] );

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
            //
            // The transformation is written here exlicitly,
            // because of the special ordering of the second derivatives,
            // and because it should be easier to extend this to the
            // convection-diffusion-reaction-equation starting from this.


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

        return hh * sum;
    }

private:

    gsMatrix<T> m_discSol2ndDer, m_discSolVals, m_rhsFctVals;

    unsigned m_parDim;

    bool m_f2param;
};


} // namespace gismo

