/** @file gsVisitorNeumannBiharmonic.h

    @brief Neumann conditions visitor for 4th order problems.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris, J. Sogn
*/

#pragma once

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
    class gsVisitorInterfaceNitscheBiharmonicStability
    {
    public:

        gsVisitorInterfaceNitscheBiharmonicStability(const gsPde<T> & )
        { }

        void initialize(const gsBasis<T> & basis1,
                        const gsBasis<T> & basis2,
                        const boundaryInterface & bi,
                        const gsOptionList & options,
                        gsQuadRule<T>    & rule)
        {
            mu = options.getReal("mu");

            side1 = bi.first().side();
            side2 = bi.second().side();
            const int dir = side1.direction();

            gsVector<int> numQuadNodes ( basis1.dim() );
            for (int i = 0; i < basis1.dim(); ++i)
                numQuadNodes[i] = basis1.degree(i) + 1;
            numQuadNodes[dir] = 1;

            // Setup Quadrature
            rule = gsGaussRule<T>(numQuadNodes);// harmless slicing occurs here

            // Set Geometry evaluation flags
            md1.flags = NEED_VALUE | NEED_2ND_DER | NEED_MEASURE | NEED_GRAD_TRANSFORM;
            md2.flags = NEED_VALUE | NEED_2ND_DER | NEED_MEASURE | NEED_GRAD_TRANSFORM;
        }

        // Evaluate on element.
        inline void evaluate(const gsBasis<T>       & basis1,
                             const gsGeometry<T>    & geo1,
                             const gsBasis<T>       & basis2,
                             const gsGeometry<T>    & geo2,
                             gsMatrix<T>            & quNodes1,
                             gsMatrix<T>            & quNodes2)
        {
            md1.points = quNodes1;
            md2.points = quNodes2;

            // Compute the active basis functions
            // Assumes actives are the same for all quadrature points on the elements
            basis1.active_into(md1.points.col(0), actives1);
            numActive1 = actives1.rows();

            basis2.active_into(md2.points.col(0), actives2);
            numActive2 = actives2.rows();

            // Evaluate basis gradients on element
            basis1.evalAllDers_into(md1.points, 2, basisData1);
            basis2.evalAllDers_into(md2.points, 2, basisData2);

            // Compute geometry related values
            geo1.computeMap(md1);
            geo2.computeMap(md2);

            // Initialize local matrix/rhs
            localMatrix1.setZero(numActive1, numActive1);
            localMatrix2.setZero(numActive2, numActive2);

            localMatrix12.setZero(numActive1, numActive2);
            localMatrix21.setZero(numActive2, numActive1);

            localRhs1.setZero(numActive1, 1);
            localRhs2.setZero(numActive2, 1);
        }

        inline void assemble(gsDomainIterator<T>    & element,
                             gsDomainIterator<T>    & ,
                             gsVector<T> const      & quWeights)
        {
            gsMatrix<T> basisVals1 = basisData1[0];
            gsMatrix<T> basisDers1 = basisData1[1];
            gsMatrix<T> basis2Ders1 = basisData1[2];

            gsMatrix<T> basisVals2 = basisData2[0];
            gsMatrix<T> basisDers2 = basisData2[1];
            gsMatrix<T> basis2Ders2 = basisData2[2];

            for (index_t k = 0; k < quWeights.rows(); ++k) // loop over quadrature nodes
            {
                // Multiply quadrature weight by the measure of normal
                const T weight = quWeights[k] * md1.measure(k);

                //Get gradients of the physical space
                transformGradients(md1, k, basisDers1, physBasisGrad1);
                transformGradients(md2, k, basisDers2, physBasisGrad2);

                // Compute physical laplacian at k as a 1 x numActive matrix
                transformLaplaceHgrad(md1, k, basisDers1, basis2Ders1, physBasisLaplace1);
                transformLaplaceHgrad(md2, k, basisDers2, basis2Ders2, physBasisLaplace2);

                const T h = element.getCellSize();
                const T mu_h = mu / (0 != h ? h : 1);

                localMatrix1.noalias() += weight * (physBasisLaplace1.transpose() * physBasisLaplace1);

                localMatrix2.noalias() += weight * (physBasisLaplace2.transpose() * physBasisLaplace2);
            }
        }


        inline void localToGlobal(const index_t                     patchIndex1,
                                  const index_t                     patchIndex2,
                                  const std::vector<gsMatrix<T> > & eliminatedDofs,
                                  gsSparseSystem<T>               & system)
        {
            // Map patch-local DoFs to global DoFs
            system.mapColIndices(actives1, patchIndex1, actives1);

            // Add contributions to the system matrix and right-hand side
            system.push(localMatrix1, localRhs1, actives1, eliminatedDofs[0], 0, 0);

            // Map patch-local DoFs to global DoFs
            system.mapColIndices(actives2, patchIndex2, actives2);

            // Add contributions to the system matrix and right-hand side
            system.push(localMatrix2, localRhs2, actives2, eliminatedDofs[0], 0, 0);
        }

    protected:


        // Neumann function
        boxSide side1, side2;

        // Basis values
        std::vector<gsMatrix<T>> basisData1, basisData2;
        gsMatrix<index_t> actives1, actives2;

        // Normal and Neumann values
        gsMatrix<T> physBasisGrad1, physBasisGrad2, physBasisLaplace1, physBasisLaplace2;

        gsVector<T> unormal1, unormal2;
        gsMatrix<T> neuData;
        index_t numActive1, numActive2;

        // Local matrix and rhs
        gsMatrix<T> localMatrix1, localMatrix2, localMatrix12, localMatrix21;
        gsMatrix<T> localRhs1, localRhs2;

        gsMapData<T> md1, md2;

        real_t mu;
    };


} // namespace gismo
