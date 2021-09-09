/** @file gsVisitorNeumannBiharmonic.h

    @brief Neumann conditions visitor for 4th order problems.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris, J. Sogn
*/

#pragma once

#include <gsC1BasisDirect/gsG1MultiBasis.h>

namespace gismo
{

/** \brief Visitor for Neumann boundary condition for the biharmonic equation.
 *
 * Visitor for boundary condition term of the form:\n
 * Let g be the function given BVP formulation, typically
 * \f[ g = \Delta u\f].
 * Then this visitor adds the follow term on the right-hand side.
 * \f[ (g, \nabla v \cdot \mathbf{n})_\Gamma \f]
 * Where v is the test function and \f[ \Gamma \f] is the boundary.
 */
    template <class T>
    class gsVisitorLaplaceBoundaryBiharmonic
    {
    public:

        gsVisitorLaplaceBoundaryBiharmonic(const gsPde<T> & , const boundary_condition<T> & s, gsG1MultiBasis<T> g1MultiBasis)
                : neudata_ptr( s.function().get() ), side(s.side()), m_g1MultiBasis(g1MultiBasis)
        { }

        gsVisitorLaplaceBoundaryBiharmonic(const gsFunction<T> & neudata, boxSide s) :
                neudata_ptr(&neudata), side(s)
        { }

        void initialize(const gsBasis<T> & basis,
                        gsQuadRule<T>    & rule)
        {
            const int dir = side.direction();
            gsVector<int> numQuadNodes ( basis.dim() );
            for (int i = 0; i < basis.dim(); ++i)
                numQuadNodes[i] = basis.degree(i) + 1;
            numQuadNodes[dir] = 1;

            // Setup Quadrature
            rule = gsGaussRule<T>(numQuadNodes);// harmless slicing occurs here

            // Set Geometry evaluation flags
            md.flags = NEED_VALUE | NEED_MEASURE | NEED_GRAD_TRANSFORM ;

        }

        void initialize(const gsBasis<T> & basis,
                        const index_t ,
                        const gsOptionList & options,
                        gsQuadRule<T>    & rule)
        {
            // Setup Quadrature (harmless slicing occurs)
            rule = gsGaussRule<T>(basis, options.getReal("quA"),
                                  options.getInt("quB"),
                                  side.direction() );

            // Set Geometry evaluation flags
            md.flags = NEED_VALUE | NEED_MEASURE | NEED_GRAD_TRANSFORM;
        }

        // Evaluate on element.
        inline void evaluate(const gsBasis<T>       & basis, // to do: more unknowns
                             const gsGeometry<T>    & geo,
                // todo: add element here for efficiency
                             gsMatrix<T>            & quNodes)
        {
            md.points = quNodes;
            // Compute the active basis functions
            // Assumes actives are the same for all quadrature points on the elements
            basis.active_into(md.points.col(0), actives);
            numActive = actives.rows();

            // Evaluate basis gradients on element
            basis.deriv_into( md.points, basisGrads);

            // Compute geometry related values
            geo.computeMap(md);

            // Evaluate the Neumann data
            neudata_ptr->eval_into(md.values[0], neuData);

            if (m_g1MultiBasis.nPatches() == 2) // TODO For now: only for two Patch domains
            {
                numg1Active = 0;
                m_g1MultiBasis.active_into(md.points.col(0), g1actives, geo.id());
                numg1Active = g1actives.rows();
                if (g1actives.rows() > 0)
                {
                    m_g1MultiBasis.evalAllDers_into(md.points, 1, g1basisData, geo.id());

                    basisGrads.conservativeResize(basisGrads.rows() + g1basisData[1].rows(), basisGrads.cols());
                    basisGrads.bottomRows(g1basisData[1].rows()) = g1basisData[1];

                    numActive += numg1Active;
                    actives.conservativeResize(actives.rows() + g1actives.rows(), actives.cols());
                    actives.bottomRows(g1actives.rows()) = g1actives;
                }
            }

            // Initialize local matrix/rhs
            localRhs.setZero(numActive, neudata_ptr->targetDim() );
        }

        inline void assemble(gsDomainIterator<T>    & ,
                             gsVector<T> const      & quWeights)
        {
            for (index_t k = 0; k < quWeights.rows(); ++k) // loop over quadrature nodes
            {
                // Compute the outer normal vector on the side
                outerNormal(md, k, side, unormal);

                // Multiply quadrature weight by the measure of normal
                const T weight = quWeights[k] * unormal.norm();
                unormal.normalize();
                //Get gradients of the physical space
                transformGradients(md, k, basisGrads, physBasisGrad);

                localRhs.noalias() += weight *(( physBasisGrad.transpose() * unormal )* neuData.col(k).transpose());
            }
        }

        inline void localToGlobal(const index_t                     patchIndex,
                                  const std::vector<gsMatrix<T> > & ,
                                  gsSparseSystem<T>               & system)
        {
            // Map patch-local DoFs to global DoFs
            system.mapColIndices(actives, patchIndex, actives);

            // Add contributions to the system matrix and right-hand side
            system.pushToRhs(localRhs, actives, 0);
        }

        /*
        void localToGlobal(const gsDofMapper & mapper,
                           const gsMatrix<T> & eliminatedDofs,
                           const index_t       patchIndex,
                           gsSparseMatrix<T> & sysMatrix,
                           gsMatrix<T>       & rhsMatrix )
        {
            // Local DoFs to global DoFs
            mapper.localToGlobal(actives, patchIndex, actives);

            // Push element contribution to the global load vector
            for (index_t j=0; j!=numActive; ++j)
            {
                // convert local DoF index to global DoF index
                const unsigned jj = actives(j);
                if (mapper.is_free_index(jj))
                {
                    rhsMatrix.row(jj) += localRhs.row(j);
                }
            }
        }
        */

    protected:


        // Neumann function
        const gsFunction<T> * neudata_ptr;
        boxSide side;

        // Basis values
        gsMatrix<T> basisGrads;
        gsMatrix<index_t> actives;

        // Normal and Neumann values
        gsMatrix<T>        physBasisGrad;

        gsVector<T> unormal;
        gsMatrix<T> neuData;
        index_t numActive;


        // Local matrix and rhs
        gsMatrix<T> localRhs;

        gsMapData<T> md;

        // Pascal
    protected:
        gsG1MultiBasis<T> m_g1MultiBasis;
        // Basis values
        std::vector<gsMatrix<T> > g1basisData;
        gsMatrix<T>        physG1BasisLaplace;
        gsMatrix<index_t> g1actives;
        index_t numg1Active;
    };


} // namespace gismo
