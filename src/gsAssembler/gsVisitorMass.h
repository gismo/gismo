/** @file gsVisitorMass.h

    @brief Mass visitor for assembling element mass matrix

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris
*/

#pragma once

namespace gismo
{
/** 
    @brief The visitor computes element mass integrals

    It sets up an assembler and assembles the mass matrix element-wise and
    combines the patch-local mass matrices into a global matrix.
      
    \ingroup Assembler
*/
template <class T>
class gsVisitorMass
{
public:

    gsVisitorMass()
    { }

    /** \brief Visitor for assembling the mass matrix
     *  
     * \f[ (u, v) \f]  
     */
    gsVisitorMass(const gsPde<T> & pde)
    { GISMO_UNUSED(pde); }

    void initialize(const gsBasis<T> & basis,
                    const index_t ,
                    const gsOptionList & options, 
                    gsQuadRule<T>    & rule)
    {
        // Setup Quadrature (harmless slicing occurs)
        rule = gsQuadrature::get(basis, options); // harmless slicing occurs here

        // Set Geometry evaluation flags
        md.flags = NEED_MEASURE;
    }

    // Evaluate on element.
    inline void evaluate(const gsBasis<T>       & basis, // to do: more unknowns
                         const gsGeometry<T>    & geo,
                         // todo: add element here for efficiency
                         gsMatrix<T>            & quNodes)
    {
        md.points = quNodes;
        // Compute the active basis functions
        // Assumes actives are the same for all quadrature points on the current element
        basis.active_into(md.points.col(0), actives);
        const index_t numActive = actives.rows();

        // Evaluate basis functions on element
        basis.eval_into(md.points, basisData);

        // Compute geometry related values
        geo.computeMap(md);

        // Initialize local matrix/rhs
        localMat.setZero(numActive, numActive);
    }

    inline void assemble(gsDomainIterator<T>    & ,
                         gsVector<T> const      & quWeights)
    {
        localMat.noalias() = 
            basisData * quWeights.asDiagonal() * 
            md.measures.asDiagonal() * basisData.transpose();
    }

    inline void localToGlobal(const index_t                     patchIndex,
                              const std::vector<gsMatrix<T> > & ,
                              gsSparseSystem<T>               & system)
    {
        // Map patch-local DoFs to global DoFs
        system.mapColIndices(actives, patchIndex, actives);

        // Add contributions to the system matrix
        system.pushToMatrix(localMat, actives, 0, 0);
    }

/* -----------------------  to be removed later*/

    void initialize(const gsBasis<T> & basis,
                           gsQuadRule<T> & rule)
    {
        gsVector<short_t> numQuadNodes( basis.dim() );
        for (short_t i = 0; i < basis.dim(); ++i)
            numQuadNodes[i] = basis.degree(i) + 1;
        
        // Setup Quadrature
        rule = gsGaussRule<T>(numQuadNodes);// harmless slicing occurs here

        // Set Geometry evaluation flags
        md.flags = NEED_MEASURE;
    }


    void localToGlobal(const gsDofMapper & mapper,
                       const gsMatrix<T> & eliminatedDofs,
                       const index_t       patchIndex,
                       gsSparseMatrix<T> & sysMatrix,
                       gsMatrix<T>       & rhsMatrix )
    {
        mapper.localToGlobal(actives, patchIndex, actives);

        const index_t numActive = actives.rows();

        for (index_t i = 0; i < numActive; ++i)
        {
            const int ii = actives(i,0); // N_i
            
            if ( mapper.is_free_index(ii) )
            {
                for (index_t j = 0; j < numActive; ++j)
                {
                    const int jj = actives(j,0); // N_j
                    if ( mapper.is_free_index(jj) )
                        //if ( jj <= ii ) // only store lower triangular part
                            sysMatrix.coeffRef(ii, jj) += localMat(i, j); // N_i*N_j
                }
            }
        }
    }


protected:

    // Basis values
    gsMatrix<T>      basisData;
    gsMatrix<index_t> actives;

    // Local matrix
    gsMatrix<T> localMat;

    gsMapData<T> md;
};


} // namespace gismo

