/** @file gsVisitorNeumann.h

    @brief Neumann conditions visitor for elliptic problems.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris
*/

#pragma once

#include <gsCore/gsFuncData.h>

namespace gismo
{

   /** @brief Implementation of a Neumann BC for elliptic assemblers.
     *
     *  The visior realizes the Neumann BC
     *  \f[ \nabla u \cdot \mathbf{n} = g_N  \f]
     *  by adding the following term to the linear form:
     *  \f[ ( g_N, v )_{\Gamma_N},  \f]
     *  where \f$ \Gamma_N \f$ is the Neumann boundary.
     *
     *  @ingroup Assembler
     */

template <class T>
class gsVisitorNeumann
{
public:

    /// @brief Constructor
    ///
    /// @param pde     Reference to \a gsPde object (is ignored)
    /// @param bc      The boundary condition to be realized
    gsVisitorNeumann(const gsPde<T> & pde, const boundary_condition<T> & bc)
    : neudata_ptr( bc.function().get() ), side(bc.side())
    { GISMO_UNUSED(pde); }

    /// @brief Constructor
    ///
    /// @param neudata Neumann boundary function
    /// @param s       Side of the geometry where Neumann BC is prescribed
    gsVisitorNeumann(const gsFunction<T> & neudata, boxSide s)
    : neudata_ptr(&neudata), side(s)
    { }

    /// Initialize
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
        md.flags = NEED_VALUE|NEED_JACOBIAN;
    }

    /// Initialize
    void initialize(const gsBasis<T>   & basis,
                    const index_t ,
                    const gsOptionList & options,
                    gsQuadRule<T>      & rule)
    {
        // Setup Quadrature (harmless slicing occurs)
        rule = gsQuadrature::get(basis, options, side.direction());

        // Set Geometry evaluation flags
        md.flags = NEED_VALUE | NEED_MEASURE | NEED_GRAD_TRANSFORM;
    }

    /// Evaluate on element
    inline void evaluate(const gsBasis<T>       & basis, // to do: more unknowns
                         const gsGeometry<T>    & geo,
                         // todo: add element here for efficiency
                         const gsMatrix<T>      & quNodes)
    {
        md.points = quNodes;
        // Compute the active basis functions
        // Assumes actives are the same for all quadrature points on the current element
        basis.active_into(md.points.col(0) , actives);
        const index_t numActive = actives.rows();

        // Evaluate basis functions on element
        basis.eval_into(md.points, basisData);

        // Compute geometry related values
        geo.computeMap(md);

        // Evaluate the Neumann data
        neudata_ptr->eval_into(md.values[0], neuData);

        // Initialize local matrix/rhs
        localRhs.setZero(numActive, neudata_ptr->targetDim() );
    }

    /// Assemble on element
    inline void assemble(gsDomainIterator<T>    & ,
                         const gsVector<T>      & quWeights)
    {
        for (index_t k = 0; k < quWeights.rows(); ++k) // loop over quadrature nodes
        {
            // Compute the outer normal vector on the side
            outerNormal(md, k, side, unormal);

            // Multiply quadrature weight by the measure of normal
            const T weight = quWeights[k] * unormal.norm();

            localRhs.noalias() += weight * basisData.col(k) * neuData.col(k).transpose() ;
        }
    }

    /// Adds the contributions to the sparse system
    inline void localToGlobal(const index_t patchIndex,
                              const std::vector<gsMatrix<T> > & ,
                              gsSparseSystem<T>               & system)
    {
        // Map patch-local DoFs to global DoFs
        system.mapColIndices(actives, patchIndex, actives);

        // Add contributions to the system matrix and right-hand side
        system.pushToRhs(localRhs, actives, 0);
    }

    /// Adds the contributions to the sparse system
    void localToGlobal(const gsDofMapper & mapper,
                       const gsMatrix<T> & eliminatedDofs,
                       const index_t       patchIndex,
                       gsSparseMatrix<T> & sysMatrix,
                       gsMatrix<T>       & rhsMatrix )
    {
        // Local DoFs to global DoFs
        mapper.localToGlobal(actives, patchIndex, actives);
        const index_t numActive = actives.rows();

        // Push element contribution to the global load vector
        for (index_t j=0; j!=numActive; ++j)
        {
            // convert local DoF index to global DoF index
            const unsigned jj = actives(j);
            if (mapper.is_free_index(jj))
                rhsMatrix.row(jj) += localRhs.row(j);
        }
    }

protected:

    // Neumann function
    const gsFunction<T> * neudata_ptr;
    boxSide side;

    // Basis values
    gsMatrix<T>      basisData;
    gsMatrix<index_t> actives;

    // Normal and Neumann values
    gsVector<T> unormal;
    gsMatrix<T> neuData;

    // Local matrix and rhs
    gsMatrix<T> localMat;
    gsMatrix<T> localRhs;

    gsMapData<T> md;
};


} // namespace gismo
