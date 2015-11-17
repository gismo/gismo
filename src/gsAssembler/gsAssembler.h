/** @file gsAssembler.h

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
#include <gsCore/gsMultiBasis.h>
#include <gsCore/gsDomainIterator.h>

#include <gsPde/gsBoundaryConditions.h>
#include <gsAssembler/gsQuadRule.h>
#include <gsAssembler/gsAssemblerOptions.h>

#include <gsAssembler/gsSparseSystem.h>

#include <gsPde/gsPde.h>

namespace gismo
{

/** @brief The assembler class provides generic routines for volume
  and boundary integrals that are used for for matrix and right-hand
  side generation
  
  \ingroup Assembler
*/
template <class T>
class gsAssembler
{
private:
    typedef gsStdVectorRef<gsDofMapper> gsDofMappers;

    typedef typename gsBoundaryConditions<T>::bcContainer bcContainer;

public: /* Constructors and initializers */

    /// @brief default constructor
    /// \note none of the data fields are inititalized, use
    /// additionally an appropriate initialize function
    gsAssembler()
    { }

    /// @brief Constructor using a multipatch domain, a
    /// vector of multibases and the boundary conditions.
    /// \note Rest of the data fields should be initialized in a
    /// derived constructor
    gsAssembler(const gsPde<T>                         & pde  ,
                const gsStdVectorRef<gsMultiBasis<T> > & bases, 
                const gsAssemblerOptions & opt = gsAssemblerOptions() )
    { 
        initialize(pde, bases, opt);
    }

    gsAssembler(const gsPde<T>        & pde,
                const gsMultiBasis<T> & bases, 
                const gsAssemblerOptions & opt = gsAssemblerOptions() )
    { 
        initialize(pde, bases, opt);
    }

    virtual ~gsAssembler()
    { }

    /// @brief Intitialize function for, sets data fields
    /// using a multi-patch, a vector of multi-basis and boundary conditions.
    /// \note Rest of the data fields should be initialized in the
    /// derived function initializePdeSpecific() .
    void initialize(const gsPde<T>                         & pde,
                    const gsStdVectorRef<gsMultiBasis<T> > & bases, 
                    const gsAssemblerOptions & opt = gsAssemblerOptions() )
    {
        // fixme: add cast op in gsMultiBasis: gsMultiBasis--> gsStdVectorRef
        // add cast operator in gsBasis: gsBasis --> gsMultiBasis -->  gsStdVectorRef
        m_pde_ptr = &pde;
        m_bases.clear();// bug ?
        m_bases = bases;
        m_options = opt;
        initializeDofs();
    }

    void initialize(const gsPde<T>           & pde,
                    const gsMultiBasis<T>    & bases, 
                    const gsAssemblerOptions & opt = gsAssemblerOptions() )
    {
        m_pde_ptr = &pde;
        m_bases.clear();// bug ?
        m_bases.push_back(bases);
        m_options = opt;
        initializeDofs();
    }

public: /* Dof initialization */

    void initializeDofs()
    {
        const gsBoundaryConditions<T> & m_bConditions = m_pde_ptr->bc();

        // Check if boundary conditions are OK
        const int np = m_bases.front().nBases();
        for (typename gsBoundaryConditions<T>::const_iterator it = 
                 m_bConditions.dirichletBegin() ; it != m_bConditions.dirichletEnd(); ++it )
            GISMO_ENSURE( it->ps.patch < np && it->ps.patch > -1,
                          "Problem: a Dirichlet boundary condition is set on a patch id which does not exist.");
        for (typename gsBoundaryConditions<T>::const_iterator it = 
                 m_bConditions.neumannBegin() ; it != m_bConditions.neumannEnd(); ++it )
            GISMO_ENSURE( it->ps.patch < np && it->ps.patch > -1,
                          "Problem: a Neumann boundary condition is set on a patch id which does not exist.");
        
        
        // option: same basis for all coordinates of unknowns
        // poisson: 1 unknown, 1 coordinate (many right-hand sides)
        // stokes: 2 unknowns, coords: d + 1, (possibly many right-hand sides)

        // for all unknowns
        // 

        // Initialize DoF mapper
        const bool conforming = ( m_options.intStrategy == iFace::glue );
        std::vector<gsDofMapper> mappers(1); // number of unknowns 
        if ( m_options.dirStrategy == dirichlet::elimination )
            m_bases.front().getMapper(conforming, m_bConditions, mappers.front() );
        else
            m_bases.front().getMapper(conforming, mappers.front() );

        m_dofs = mappers.front().freeSize();

        m_system = gsSparseSystem<T>(mappers, 
                                     m_pde_ptr->numUnknowns() , 
                                     m_pde_ptr->numUnknowns() );
    }

    void finalize()
    {
        m_system.matrix().makeCompressed();
    }

public: /* Element visitors */

    /// @brief Iterates over all elements of the domain and applies
    /// the \a ElementVisitor
    template<class ElementVisitor>
    void push()
    {
        for (unsigned np=0; np < m_pde_ptr->domain().nPatches(); ++np )
        {
            ElementVisitor visitor(*m_pde_ptr);
            //Assemble (fill m_matrix and m_rhs) on patch np
            apply(visitor, np);
        }
    }

    /// @brief Iterates over all elements of the boundaries \a BCs and
    /// applies the \a ElementVisitor
    template<class BElementVisitor>
    void push(const bcContainer & BCs) // ---------------- push()
    {
        for (typename bcContainer::const_iterator it 
                 = BCs.begin(); it!= BCs.end(); ++it)
        {
            BElementVisitor visitor(*m_pde_ptr, *it);
            //Assemble (fill m_matrix and m_rhs) contribution from this BC
            apply(visitor, it->patch(), it->side());
        }
    }

    /// @brief Iterates over all elements of the domain and applies
    /// the \a ElementVisitor
    template<class ElementVisitor>
    void push(const ElementVisitor & visitor)
    {
        for (unsigned np=0; np < m_pde_ptr->domain().nPatches(); ++np )
        {
            ElementVisitor curVisitor = visitor;
            //Assemble (fill m_matrix and m_rhs) on patch np
            apply(curVisitor, np);
        }
    }

public:  /* Virtual assembly routines */

    /// @brief Main assemble routine
    virtual void assemble();

    /// @brief Main non-linear assemble routine with input from
    /// current solution
    virtual void assemble(const gsMultiPatch<T> & curSolution);

public:  /* Dirichlet degrees of freedom computation */

    /// @brief forces the Assembler to calculete the Dirichlet dofs again.
    void computeDirichletDofs();

    void setDirichletDofs(const gsMatrix<T> & vals, int unk = 0);

    void computeDirichletDofsIntpl();
    
    void computeDirichletDofsL2Proj();

    /// Enforce Dirichlet boundary conditions by diagonal penalization
    void penalizeDirichlet();

    // index_t numFixedDofs(int unk = 0) {return m_dofMappers[unk].boundarySize();}

public:  /* Dirichlet degrees of freedom computation */

    /// @brief Reconstruct solution from computed solution vector
    virtual void constructSolution(const gsMatrix<T>& solVector, 
                                   gsMultiPatch<T>& result, int unk = 0) const;

    // fixme: remove (currently for compatibility)
    GISMO_DEPRECATED
    gsField<T> *  constructSolution(const gsMatrix<T>& solVector, int unk = 0) const;

    /// @brief Update solution by adding the computed solution vector
    /// to the current solution
    virtual void updateSolution(const gsMatrix<T>& solVector, 
                                gsMultiPatch<T>& result) const
    {GISMO_NO_IMPLEMENTATION}


public:

    /// @brief Return the Pde
    const gsPde<T> & pde() const { return *m_pde_ptr; }

    /// @brief Return the multipatch
    const gsMultiPatch<T> & patches() const { return m_pde_ptr->patches(); }

    /// @brief Return the multi-basis
    const gsMultiBasis<T> & multiBasis(index_t k) const { return m_bases[k]; }

    /// @brief Returns the number of multi-bases
    std::size_t numMultiBasis() const {return m_bases.size(); }

    /// @brief Returns the left-hand global matrix
    const gsSparseMatrix<T> & matrix() const { return m_system.matrix(); }

    /// @brief Returns the left-hand side vector(s)
    /// ( multiple right hand sides possible )
    const gsMatrix<T> & rhs() const { return m_system.rhs(); }

    /// @brief Returns the left-hand global matrix
    const gsSparseSystem<T> & system() const { return m_system; }

    /// @brief Allocates sparse matrix and right-hand side
    void reserveSparseSystem() 
    {
        int nz = 1;
        for (int i = 0; i != m_bases.front().dim(); ++i) // to do: improve
            nz *= m_options.bdA * m_bases.front().maxDegree(i) + m_options.bdB; // "3 *" for dg
        
        m_system.reserve(nz, m_pde_ptr->numRhs());
    }

    /// @brief Allocates sparse matrix and right-hand side
    void setSparseSystem(gsSparseSystem<T> & sys) 
    {
        m_system.swap(sys);
    }

    // Returns true if only the lower triangular part of the
    // matrix is computed (for symmetric problems)
    //virtual bool isSymmertric() const { return false; }

    /// @brief Returns the Dirichlet values (if applicable)
    const gsMatrix<T> & dirValues() const { return m_ddof; }

    /// @brief Sets the Dirichlet values to given \a values (if applicable)
    void setDirichletValues(gsMatrix<T> values) { m_ddof.swap(values); }

    /// @brief Sets any Dirichlet values to homogeneous (if applicable)
    void homogenizeDirichlet() { m_ddof.setZero(); }

    /// @brief Returns the number of (free) degrees of freedom
    int numDofs() const { return m_dofs; }

    /// @brief Returns the options of the assembler
    const gsAssemblerOptions & options() const { return m_options; }


protected:

    /// @brief Generic assembly routine for volume or boundary integrals
    template<class ElementVisitor>
    void apply(ElementVisitor & visitor, 
               int patchIndex = 0, 
               boxSide side = boundary::none);
    
protected:

    /// The PDE: contains multi-patch domain, boundary conditions and
    /// coeffcient functions
    const gsPde<T> * m_pde_ptr;

    /// The discretization bases corresponding to \a m_patches and to
    /// the number of solution fields that are to be computed
    /// m_bases[i]: The multi-basis for unknown i
    std::vector< gsMultiBasis<T> > m_bases;
    
protected:

    /// Fixed DoF values (if applicable, for instance eliminated Dirichlet DoFs)
    gsMatrix<T> m_ddof;  // fixme: std::vector<gsMatrix<T> > m_fixedDofs;

    /// Options
    gsAssemblerOptions m_options;

protected: // *** Outputs *** 
    
    /// Global sparse linear system
    gsSparseSystem<T> m_system;

protected: // *** Information *** 

    /// Number of degrees of freedom (excluding eliminated,
    /// counting glued dofs only once, etc.)
    int m_dofs;
};


template <class T>
template<class ElementVisitor>
void gsAssembler<T>::apply(ElementVisitor & visitor, 
                           int patchIndex, 
                           boxSide side)
{
    //gsDebug<< "Apply to patch "<< patchIndex <<"("<< side <<")\n";
    
    const gsBasisRefs<T> bases(m_bases, patchIndex   );
    
    gsQuadRule<T> QuRule ; // Quadrature rule
    gsMatrix<T> quNodes  ; // Temp variable for mapped nodes
    gsVector<T> quWeights; // Temp variable for mapped weights
    unsigned evFlags(0);
    
    // Initialize reference quadrature rule and visitor data
    visitor.initialize(bases, patchIndex, m_options, QuRule, evFlags);
    
    //fixme: gsMapData<T> mapData;
    // Initialize geometry evaluator 
    typename gsGeometry<T>::Evaluator geoEval( 
        m_pde_ptr->patches()[patchIndex].evaluator(evFlags));
    
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
        //visitor.localToGlobal(mappers, m_ddof, patchIndex, m_system.matrix(), m_system.rhs());
        visitor.localToGlobal(patchIndex, m_ddof, m_system);
    }
}

} // namespace gismo

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsAssembler.hpp)
#endif
