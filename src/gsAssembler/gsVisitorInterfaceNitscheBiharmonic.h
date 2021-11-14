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
    class gsVisitorInterfaceNitscheBiharmonic
    {
    public:

        gsVisitorInterfaceNitscheBiharmonic(const gsPde<T> & )
        { }

        void initialize(const gsBasis<T> & basis1,
                        const gsBasis<T> & basis2,
                        const boundaryInterface & bi,
                        const gsOptionList & options,
                        gsQuadRule<T>    & rule)
        {
            mu = options.getReal("mu");

            side1 = bi.first().side();
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
            /*
            localMatrix1.setZero(numActive1, numActive1);
            localMatrix2.setZero(numActive2, numActive2);

            localMatrix12.setZero(numActive1, numActive2);
            localMatrix21.setZero(numActive2, numActive1);
*/
            localRhs1.setZero(numActive1, 1);
            localRhs2.setZero(numActive2, 1);

            // Initialize local matrices
            B11.setZero(numActive1, numActive1); B12.setZero(numActive1, numActive2);
            E11.setZero(numActive1, numActive1); E12.setZero(numActive1, numActive2);
            B22.setZero(numActive2, numActive2); B21.setZero(numActive2, numActive1);
            E22.setZero(numActive2, numActive2); E21.setZero(numActive2, numActive1);
        }

        inline void assemble(gsDomainIterator<T>    & element1,
                             gsDomainIterator<T>    & element2,
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
                // Compute the outer normal vector on the side
                outerNormal(md1, k, side1, unormal1);
                //outerNormal(md2, k, side2, unormal2);

                // Multiply quadrature weight by the measure of normal
                const T weight = quWeights[k] * unormal1.norm();
                unormal1.normalize();
                //unormal2.normalize();

                //Get gradients of the physical space
                transformGradients(md1, k, basisDers1, physBasisGrad1);
                transformGradients(md2, k, basisDers2, physBasisGrad2);

                // Compute physical laplacian at k as a 1 x numActive matrix
                transformLaplaceHgrad(md1, k, basisDers1, basis2Ders1, physBasisLaplace1);
                transformLaplaceHgrad(md2, k, basisDers2, basis2Ders2, physBasisLaplace2);

                // Compute element matrices
                const T c1     = weight * T(0.5);
                N1.noalias()   = physBasisGrad1.transpose()*unormal1;
                N2.noalias()   = physBasisGrad2.transpose()*unormal1;

                B11.noalias() += c1 * ( N1 * physBasisLaplace1 );
                B12.noalias() += c1 * ( N1 * physBasisLaplace2 );
                B22.noalias() -= c1 * ( N2 * physBasisLaplace2 );
                B21.noalias() -= c1 * ( N2 * physBasisLaplace1 );

                //const T h = element.getCellSize();
                //const T mu_h = mu / (0 != h ? h : 1);

                const T h1     = element1.getCellSize();
                const T h2     = element2.getCellSize();
                // Maybe, the h should be scaled with the patch diameter, since its the h from the parameterdomain.
                const T c2     = weight * mu * 2*(1./h1 + 1./h2);

                E11.noalias() += c2 * ( N1 * N1.transpose() );
                E12.noalias() += c2 * ( N1 * N2.transpose() );
                E22.noalias() += c2 * ( N2 * N2.transpose() );
                E21.noalias() += c2 * ( N2 * N1.transpose() );

/*
                localMatrix1.noalias() += weight * (physBasisGrad1.transpose() * unormal1 * physBasisLaplace1
                                                    + (physBasisGrad1.transpose() * unormal1 * physBasisLaplace1).transpose()
                                                - mu_h * (physBasisGrad1.transpose() * unormal1) * (physBasisGrad1.transpose() * unormal1).transpose());

                localMatrix2.noalias() += weight * (physBasisGrad2.transpose() * unormal2 * physBasisLaplace2
                                                    + (physBasisGrad2.transpose() * unormal2 * physBasisLaplace2).transpose()
                                                    - mu_h * (physBasisGrad2.transpose() * unormal2) * (physBasisGrad2.transpose() * unormal2).transpose());

                localMatrix12.noalias() += weight * (physBasisGrad1.transpose() * unormal1 * physBasisLaplace2
                                                    + (physBasisGrad2.transpose() * unormal2 * physBasisLaplace1).transpose()
                                                    - mu_h * (physBasisGrad1.transpose() * unormal1) * (physBasisGrad2.transpose() * unormal2).transpose());

                localMatrix21.noalias() += weight * ( physBasisGrad2.transpose() * unormal2 * physBasisLaplace1
                                                     +  (physBasisGrad1.transpose() * unormal1 * physBasisLaplace2).transpose()
                                                     - mu_h * (physBasisGrad2.transpose() * unormal2) * (physBasisGrad1.transpose() * unormal1).transpose());

*/
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
            //system.push(localMatrix1, localRhs1, actives1, eliminatedDofs[0], 0, 0);

            // Map patch-local DoFs to global DoFs
            system.mapColIndices(actives2, patchIndex2, actives2);

            // Add contributions to the system matrix and right-hand side
            //system.push(localMatrix2, localRhs2, actives2, eliminatedDofs[0], 0, 0);

            //system.push(localMatrix12, localRhs1, actives1, actives2, eliminatedDofs[0], 0, 0);
            //system.push(localMatrix21, localRhs2, actives2, actives1, eliminatedDofs[0], 0, 0);

            system.push(-B11 - B11.transpose() + E11, localRhs1,actives1,actives1,eliminatedDofs.front(),0,0);
            system.push(-B21 - B12.transpose() - E21, localRhs2,actives2,actives1,eliminatedDofs.front(),0,0);
            system.push(-B12 - B21.transpose() - E12, localRhs1,actives1,actives2,eliminatedDofs.front(),0,0);
            system.push(-B22 - B22.transpose() + E22, localRhs2,actives2,actives2,eliminatedDofs.front(),0,0);
        }

    protected:

        gsVector<> m_valuePenalty;

        // Neumann function
        boxSide side1;

        // Basis values
        std::vector<gsMatrix<T>> basisData1, basisData2;
        gsMatrix<index_t> actives1, actives2;

        // Normal and Neumann values
        gsMatrix<T> physBasisGrad1, physBasisGrad2, physBasisLaplace1, physBasisLaplace2;

        gsVector<T> unormal1;
        index_t numActive1, numActive2;

        // Local matrix and rhs
        gsMatrix<T> localRhs1, localRhs2;

        // Auxiliary element matrices
        gsMatrix<T> B11, B12, E11, E12, N1,
                B22, B21, E22, E21, N2;

        gsMapData<T> md1, md2;

        real_t mu;
    };


} // namespace gismo
