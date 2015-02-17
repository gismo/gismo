/** @file gsAssemblerBase.h

    @brief Provides generic assembler routines

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris
*/

#pragma once

#include <gsCore/gsForwardDeclarations.h>
#include <gsCore/gsBasisRefs.h>
#include <gsCore/gsDofMapper.h>
#include <gsCore/gsStdVectorRef.h>

#include <gsAssembler/gsQuadRule.h>

namespace gismo
{

/** @brief The assembler class provides to generic routines for volume
 * and boundary integrals that are used for for matrix and rhs
 * generation
 *
 * \ingroup Assembler
*/
template <class T>
class gsAssemblerBase
{
private:
    typedef gsStdVectorRef<gsDofMapper> gsDofMappers;

public:

    /// @brief Constructor using a multipatch domain
    /// \note Rest of the data fields should be initialized in a
    /// derived constructor
    gsAssemblerBase(const gsMultiPatch<T> & patches) :
    m_patches(patches)
    { }

    ~gsAssemblerBase()
    { }

    /// @brief Generic assembly routine for volume or boundary integrals
    template<class ElementVisitor>
    void apply(ElementVisitor & visitor, 
               int patchIndex = 0, 
               boxSide side = boundary::none)
    {
        //gsDebug<< "Apply to patch "<< patchIndex <<"("<< side <<")\n";

        const gsBasisRefs<T> bases(m_bases, patchIndex);
        const gsDofMappers mappers(m_dofMappers);
        
        gsMatrix<T> quNodes  ; // Temp variable for mapped nodes
        gsVector<T> quWeights; // Temp variable for mapped weights
        gsQuadRule<T> QuRule; // Reference Quadrature rule
        unsigned evFlags(0);

        // Initialize 
        visitor.initialize(bases, QuRule, evFlags);

        // Initialize geometry evaluator -- TODO: Initialize inside visitor
        typename gsGeometry<T>::Evaluator geoEval( 
            m_patches[patchIndex].evaluator(evFlags));

        // Initialize domain element iterator -- using unknown 0
        typename gsBasis<T>::domainIter domIt = bases[0].makeDomainIterator(side);

        // Start iteration over elements
        for (; domIt->good(); domIt->next() )
        {
            // Map the Quadrature rule to the element
            QuRule.mapTo( domIt->lowerCorner(), domIt->upperCorner(), quNodes, quWeights );

            // Perform required evaluations on the quadrature nodes
            visitor.evaluate(bases, /* *domIt,*/ *geoEval, quNodes);
            
            // Assemble on element
            visitor.assemble(*domIt, *geoEval, quWeights);
            
            // Push to global matrix and right-hand side vector
            visitor.localToGlobal(mappers, m_ddof, patchIndex, m_matrix, m_rhs);
        }
    }


public:

    /// @brief Return the multipatch.
    const gsMultiPatch<T> & patches() const { return m_patches; }

    /// @brief Return the multi-basis
    const gsMultiBasis<T> & multiBasis(index_t k) const { return m_bases[k]; }

    /// @brief Return the DOF mapper for unknown \em i.
    const gsDofMapper& dofMapper(unsigned i = 0) const     { return m_dofMappers[i]; }

    /// @brief Returns the left-hand global matrix
    const gsSparseMatrix<T> & matrix() const { return m_matrix; }

    /// @brief Returns true if only the lower triangular part of the
    /// matrix is computed (for symmetric problems)
    virtual bool isSymmertric() const { return false; }

    /// @brief Returns the left-hand side vector(s)
    /// ( multiple right hand sides possible )
    const gsMatrix<T> & rhs() const { return m_rhs; }

    /// @brief Returns the Dirichlet values (if applicable)
    const gsMatrix<T> & dirValues() const { return m_ddof; }

    /// @brief Returns the number of (free) degrees of freedom
    int numDofs() const { return m_dofs; }

protected:

    /// @brief The multipatch domain
    gsMultiPatch<T> m_patches;

    /// The discretization bases corresponding to \a m_patches and to
    /// the number of solution fields that are to be computed
    /// m_bases[i]: The multi-basis for unknown i
    std::vector< gsMultiBasis<T> > m_bases;
    
    /// The Dof mapper is used to map patch-local DoFs to the global DoFs
    /// One for each unknown, one for each patch
    /// m_dofMappers[i]: DoF Mapper for unknown i
    std::vector<gsDofMapper>  m_dofMappers;
    
    /// @brief Dirichlet DoF fixed values (if applicable)
    gsMatrix<T> m_ddof;

    /// @brief Reference Quadrature rule
    gsQuadRule<T> QuRule;

protected:
    // *** Outputs *** 
    
    /// @brief Global system matrix
    gsSparseMatrix<T> m_matrix;

    /// @brief Right-hand side ( multiple right hand sides possible )
    gsMatrix<T>       m_rhs;

protected:
    // *** Information *** 

    /// @brief Number of degrees of freedom (excluding eliminated,
    /// counting glued dofs only once, etc.)
    int m_dofs;

};


} // namespace gismo

