/** @file gsVisitorNitscheBiharmonic.h

    @brief First-type Nitsche BC imposition visitor for 
           the biharmonic problem.

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): S. Moore
*/

#pragma once

namespace gismo
{
/** @brief
    Visitor for the first-type Nitsche BC of the biharmonic problem.

    It sets up an assembler and assembles the system patch wise and
    combines the patch-local stiffness matrices into a global system.
    Dirichlet boundary can also be imposed weakly (i.e Nitsche ) 
*/

template <class T>
class gsVisitorNitscheBiharmonic
{
public:
/** \brief Visitor for the weak imposition of the first-type dirichlet 
 *         boundary condition.
 *
 * The visitor adds this term to the bilinear term
 * \f[ (\Delta u, \nabla v)_{\partial \Omega} + (\nabla u, \Delta v )_{\partial \Omega} 
 *       + (\mu* \nabla u, \nabla v)_{\partial \Omega} \f]
 * 
 * and the following term is also added to the linear term
 * \f[ (g_1, \mu* \nabla v + \Delta v)_{\partial \Omega}  \f],
 *
 */
    gsVisitorNitscheBiharmonic(const gsFunction<T> & dirdata, T _penalty, boxSide s) : 
    dirdata_ptr(&dirdata),penalty(_penalty), side(s)
    { }

    void initialize(const gsBasis<T> & basis, 
                    gsQuadRule<T>    & rule)
    {
        const int dir = side.direction();
        gsVector<int> numQuadNodes(basis.dim());
        for (int i = 0; i < basis.dim(); ++i)
            numQuadNodes[i] = basis.degree(i) + 1;
        numQuadNodes[dir] = 1;

        // Setup Quadrature
        rule = gsGaussRule<T>(numQuadNodes);// harmless slicing occurs here

        // Set Geometry evaluation flags
        md.flags = NEED_VALUE | NEED_MEASURE | NEED_GRAD_TRANSFORM | NEED_2ND_DER;
    }

    // Evaluate on element.
    inline void evaluate(const gsBasis<T>    & basis, // to do: more unknowns
                         const gsGeometry<T> & geo,
                         // todo: add element here for efficiency
                         gsMatrix<T>         & quNodes)
    {
        md.points = quNodes;
        // Compute the active basis functions
        // Assumes actives are the same for all quadrature points on the current element
        basis.active_into(md.points.col(0), actives);
        const index_t numActive = actives.rows();

        // Evaluate basis values and derivatives on element
        basis.evalAllDers_into(md.points, 2, basisData);

        // Compute geometry related values
        geo.computeMap(md);

        // Evaluate the Dirichlet data
        dirdata_ptr->eval_into(md.values[0], dirData);

        // Initialize local matrix/rhs
        localMat.setZero(numActive, numActive);
        localRhs.setZero(numActive, dirdata_ptr->targetDim());
    }

    inline void assemble(gsDomainIterator<T> & element,
                         const gsVector<T>   & quWeights)
    {
        const unsigned d = element.dim();

        const index_t numActive = actives.rows();

        for (index_t k = 0; k < quWeights.rows(); ++k) // loop over quadrature nodes
        {
            //const typename gsMatrix<T>::Block basisVals  = basisData.topRows(numActive);
            const typename gsMatrix<T>::Block basisGrads =
                basisData.middleRows(numActive, numActive * d);
            const typename gsMatrix<T>::Block basis2ndDerivs =
                basisData.bottomRows(numActive * (d * (d + 1)) / 2);

            // Compute the outer normal vector on the side
            outerNormal(md, k, side, unormal);

            // Multiply quadrature weight by the geometry measure
            const T weight = quWeights[k] * unormal.norm();

            // Compute the unit normal vector : Dim x 1
            unormal.normalize();

            // Compute physical gradients at k as a Dim x NumActive matrix
            transformGradients(md, k, basisGrads, physBasisGrads);
            // Compute physical laplacian at k as a 1 x numActive matrix
            transformLaplaceHgrad(md, k, basisGrads, basis2ndDerivs, physBasisLaplace);

            // Get penalty parameter
            const T h = element.getCellSize();
            const T mu = penalty / (0 != h ? h : 1);

            // Sum up quadrature point evaluations
            localRhs.noalias() += weight * ((physBasisLaplace.transpose() + mu * physBasisGrads.transpose() * unormal)
                * dirData.col(k).transpose());

            localMat.noalias() += weight * (physBasisGrads.transpose() * unormal * physBasisLaplace
                + (physBasisGrads.transpose() * unormal * physBasisLaplace).transpose()
                - mu * physBasisGrads.transpose() * physBasisGrads);
        }
    }
    
    void localToGlobal(const gsDofMapper & mapper,
                       const gsMatrix<T> & eliminatedDofs,
                       const index_t       patchIndex,
                       gsSparseMatrix<T> & sysMatrix,
                       gsMatrix<T>       & rhsMatrix )
    {
        // Local DoFs to global DoFs
        mapper.localToGlobal(actives, patchIndex, actives);
        const index_t numActive = actives.rows();

        // Push element contribution to the global matrix and load vector
        for (index_t j=0; j!=numActive; ++j)
        {
            const unsigned jj = actives(j);
            rhsMatrix.row(jj) -= localRhs.row(j);
            for (index_t i=0; i!=numActive; ++i)
            {
                const unsigned ii = actives(i);
                if ( jj <= ii ) // assuming symmetric problem (!) probably we must not.
                    sysMatrix( ii, jj ) -= localMat(i,j);
            }
        }

    }

private:
    // Dirichlet function
    const gsFunction<T> * dirdata_ptr;

    // Penalty constant
    T penalty;

    // Side
    boxSide side;

private:
    // Basis values
    gsMatrix<T>       basisData;
    gsMatrix<T>       physBasisGrads;
    gsMatrix<T>       physBasisLaplace;
    gsMatrix<index_t> actives;

    // Normal and Neumann values
    gsVector<T> unormal;
    gsMatrix<T> dirData;

    // Local  matrix and rhs
    gsMatrix<T> localMat;
    gsMatrix<T> localRhs;

    gsMapData<T> md;
};


} // namespace gismo
