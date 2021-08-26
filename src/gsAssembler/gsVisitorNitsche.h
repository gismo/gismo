/** @file gsVisitorNitsche.h

    @brief Weak (Nitsche-type) imposition of the Dirichlet boundary conditions.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzflaris, S. Moore, S. Takacs
*/

#pragma once

#include <gsAssembler/gsVisitorDg.h>

namespace gismo
{

   /** @brief Visitor for adding the terms associated to weak (Nitsche-type) imposition
     * of the Dirichlet boundary conditions.
     *
     * This visitor adds the following term to the left-hand side (bilinear form):
     * \f[
     *          \alpha \int_{\Gamma_D}  \nabla u \cdot \mathbf{n} v
     *        + \beta  \int_{\Gamma_D}  \nabla v \cdot \mathbf{n} u
     *        + \delta h^{-1} \int_{\Gamma_D}  u  v
     * \f]
     * and the following terms to the right-hand side (linear form):
     * \f[
     *          \beta   \int_{\Gamma_D} g_D \nabla v \cdot \mathbf{n} u
     *        + \delta h^{-1} \int_{\Gamma_D} g_D v \cdot \mathbf{n} u,
     * \f]
     * where \f$g_D\f$ is the Dirichlet value.
     *
     * The default values for Nitsche.Alpha and Nitsche.Beta are \f$\alpha=\beta=1\f$,
     * which corresponds to the standard Nitsche method. The default value for
     * Nitsche.Penalty is -1, which yields \f$\delta=2.5(p+d)(p+1)\f$. If a positive
     * value for Nitsche.Penalty is chosen, that value is taken for \f$\delta\f$.
     *
     * The (normal) grid sizes are estimated based on the geometry function. Use
     * the option Nitsche.ParameterGridSize to just use the grid size on the parameter
     * domain.
     *
     * An analogous visitor for handling the internal smoothness weakly, is the
     * \a gsVisitorDg.
     *
     * @ingroup Assembler
     */


template <class T>
class gsVisitorNitsche
{
public:

    /// @brief Constructor
    ///
    /// @param pde     Reference to \a gsPde object
    /// @param bc      The boundary condition to be realized
    gsVisitorNitsche(const gsPde<T> & pde, const boundary_condition<T> & bc)
    : m_pde(&pde), m_dirdata_ptr( bc.function().get() ), m_penalty(-1), m_side(bc.ps)
    { }

    /// Default options
    static gsOptionList defaultOptions()
    {
        gsOptionList options;
        options.addReal  ("Nitsche.Alpha",
                          "Parameter alpha for Nitsche scheme.",                                                       1);
        options.addReal  ("Nitsche.Beta",
                          "Parameter beta for Nitsche scheme.",                                                        1);
        options.addReal  ("Nitsche.Penalty",
                          "Penalty parameter penalty for Nitsche scheme; if negative, default 2.5(p+d)(p+1) is used.",-1);
        options.addSwitch("Nitsche.ParameterGridSize",
                          "Use grid size on parameter domain for the penalty term.",                               false);
        return options;
    }

    /// Initialize
    void initialize(const gsBasis<T> & basis,
                    const index_t,
                    const gsOptionList & options,
                    gsQuadRule<T> & rule)
    {
        // Setup Quadrature (harmless slicing occurs)
        rule = gsQuadrature::get(basis, options, m_side.direction());

        m_penalty     = options.askReal("Nitsche.Penalty",-1);
        // If not given, use default
        if (m_penalty<0)
        {
            const index_t deg = basis.maxDegree();
            m_penalty = T(2.5) * (deg + basis.dim()) * (deg + 1);
        }

        m_alpha     = options.askReal("Nitsche.Alpha", 1);
        m_beta      = options.askReal("Nitsche.Beta" , 1);

        if (options.getSwitch("DG.ParameterGridSize"))
        {
            m_h     = basis.getMinCellLength();
        }
        else
        {
            GISMO_ENSURE (m_pde, "gsVisitorNitsche::initialize: No PDE given.");
            m_h     = gsVisitorDg<T>::estimateSmallestPerpendicularCellSize(
                          basis,
                          m_pde->patches()[m_side.patch],
                          m_side
                      );
        }

        // Set Geometry evaluation flags
        md.flags = NEED_VALUE | NEED_MEASURE | NEED_GRAD_TRANSFORM;

    }

    /// Evaluate on element
    inline void evaluate(const gsBasis<T>       & basis,
                         const gsGeometry<T>    & geo,
                         // todo: add element here for efficiency
                         const gsMatrix<T>      & quNodes)
    {
        md.points = quNodes;
        // Compute the active basis functions
        // Assumes actives are the same for all quadrature points on the current element
        basis.active_into(md.points.col(0) , actives);
        const index_t numActive = actives.rows();

        // Evaluate basis values and derivatives on element
        basis.evalAllDers_into( md.points, 1, basisData);

        // Compute geometry related values
        geo.computeMap(md);

        // Evaluate the Dirichlet data
        m_dirdata_ptr->eval_into(md.values[0], dirData);

        // Initialize local matrix/rhs
        localMat.setZero(numActive, numActive);
        localRhs.setZero(numActive, m_dirdata_ptr->targetDim() );
    }

    /// Assemble on element
    inline void assemble(gsDomainIterator<T>    & element,
                         const gsVector<T>      & quWeights)
    {
        gsMatrix<T> & bGrads = basisData[1];
        const index_t numActive = actives.rows();

        for (index_t k = 0; k < quWeights.rows(); ++k) // loop over quadrature nodes
        {

            const typename gsMatrix<T>::Block bVals =
                basisData[0].block(0,k,numActive,1);

            // Compute the outer normal vector on the side
            outerNormal(md, k, m_side, unormal);

            // Multiply quadrature weight by the geometry measure
            const T weight = quWeights[k] * unormal.norm();

            // Compute the unit normal vector
            unormal.normalize();

            // Compute physical gradients at k as a Dim x NumActive matrix
            transformGradients(md, k, bGrads, pGrads);

            // Get penalty parameter
            const T mu = m_penalty / m_h;

            // Sum up quadrature point evaluations
            localRhs.noalias() -= weight * ( ( m_beta * pGrads.transpose() * unormal - mu * bVals ) * dirData.col(k).transpose() );

            localMat.noalias() -= weight * (
                                            m_alpha *  bVals * unormal.transpose() * pGrads
                                            + m_beta  * (bVals * unormal.transpose() * pGrads).transpose()
                                            - mu      *  bVals                       * bVals.transpose()
                                        );
        }
    }

    /// Adds the contributions to the sparse system
    inline void localToGlobal(const index_t                     patchIndex,
                              const std::vector<gsMatrix<T> > & ,
                              gsSparseSystem<T>               & system)
    {
        // Map patch-local DoFs to global DoFs
        system.mapColIndices(actives, patchIndex, actives);

        // Add contributions to the system matrix and right-hand side
        system.pushAllFree(localMat, localRhs, actives, 0);
    }

    /// Adds the contributions to the sparse system
    void localToGlobal(const gsDofMapper  & mapper,
                       const gsMatrix<T>  & eliminatedDofs,
                       const index_t        patchIndex,
                       gsSparseMatrix<T>  & sysMatrix,
                       gsMatrix<T>        & rhsMatrix )
    {
        // Local DoFs to global DoFs
        mapper.localToGlobal(actives, patchIndex, actives);
        const index_t numActive = actives.rows();

        // Push element contribution to the global matrix and load vector
        for (index_t j=0; j!=numActive; ++j)
        {
            const unsigned jj = actives(j);
            rhsMatrix.row(jj) += localRhs.row(j);
            for (index_t i=0; i!=numActive; ++i)
            {
                const unsigned ii = actives(i);
//                if ( jj <= ii ) // assuming symmetric problem
                    sysMatrix( ii, jj ) += localMat(i,j);
            }
        }

    }

private:

    /// The underlying PDE
    const gsPde<T> * m_pde;

    /// Dirichlet function
    const gsFunction<T> * m_dirdata_ptr;

    /// Parameter \f$\alpha\f$ for the linear form
    T m_alpha;

    /// Parameter \f$\beta\f$ for the linear form
    T m_beta;

    /// Parameter \f$\delta\f$ for the linear form
    T m_penalty;

    /// Patch side
    patchSide m_side;


private:
    // Basis values
    std::vector<gsMatrix<T> > basisData;
    gsMatrix<T>      pGrads;
    gsMatrix<index_t> actives;

    // Normal and Neumann values
    gsVector<T> unormal;
    gsMatrix<T> dirData;

    // Local  matrix and rhs
    gsMatrix<T> localMat;
    gsMatrix<T> localRhs;

    gsMapData<T> md;

    // Grid size
    T m_h;
};


} // namespace gismo
