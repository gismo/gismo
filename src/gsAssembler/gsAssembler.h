/** @file gsAssembler.h

    @brief Provides generic assembler routines

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): S. Kleiss, A. Mantzaflaris, J. Sogn
*/

#pragma once

#include <gsCore/gsForwardDeclarations.h>
#include <gsCore/gsBasisRefs.h>
#include <gsCore/gsDofMapper.h>
#include <gsCore/gsStdVectorRef.h>
#include <gsCore/gsMultiBasis.h>
#include <gsCore/gsDomainIterator.h>
#include <gsCore/gsAffineFunction.h>

#include <gsIO/gsOptionList.h>

#include <gsPde/gsPde.h>
#include <gsPde/gsBoundaryConditions.h>

#include <gsAssembler/gsQuadRule.h>
#include <gsAssembler/gsSparseSystem.h>



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

protected: // *** Input data members *** 

    /// The PDE: contains multi-patch domain, boundary conditions and
    /// coeffcient functions
    memory::shared_ptr<gsPde<T> > m_pde_ptr;

    /// The discretization bases corresponding to the patches and to
    /// the number of solution fields that are to be computed. One
    /// entry of the vector corresponds to one or more unknown (in case
    /// of a system of PDEs)
    std::vector< gsMultiBasis<T> > m_bases;
    
    /// Options
    gsOptionList m_options;

protected: // *** Output data members *** 
    
    /// Global sparse linear system
    gsSparseSystem<T> m_system;

    /// Fixed DoF values (if applicable, for instance eliminated Dirichlet DoFs). One
    /// entry of the vector corresponds to exactly one unknown, i.e. the size of m_ddof
    /// must fit m_system.colBlocks().
    std::vector<gsMatrix<T> > m_ddof;

public:

    gsAssembler() : m_options(defaultOptions())
    { }

    virtual ~gsAssembler()
    { }

    /// @brief Returns the list of default options for assembly
    static gsOptionList defaultOptions();
    
    /// @brief Create an empty Assembler of the derived type and return a
    /// pointer to it. Call the initialize functions to set the members
    virtual gsAssembler * create() const;

    /// @brief Clone this Assembler, making a deep copy.
    virtual gsAssembler * clone() const;

    /// @brief Intitialize function for, sets data fields
    /// using the pde, a vector of multi-basis and assembler options.
    void initialize(const gsPde<T>                         & pde,
                    const gsStdVectorRef<gsMultiBasis<T> > & bases,
                    //note: merge with initialize(.., const gsMultiBasis<T> ,..) ?
                    const gsOptionList & opt = defaultOptions() )
    {
        typename gsPde<T>::Ptr _pde = memory::make_shared_not_owned(&pde);
        initialize(_pde,bases,opt);
    }

    /// @brief Intitialize function for, sets data fields
    /// using the pde, a vector of multi-basis and assembler options.
    void initialize(typename gsPde<T>::Ptr pde,
                    const gsStdVectorRef<gsMultiBasis<T> > & bases,
                    //note: merge with initialize(.., const gsMultiBasis<T> ,..) ?
                    const gsOptionList & opt = defaultOptions() )
    {
        m_pde_ptr = pde;
        m_bases = bases;
        m_options = opt;
        refresh(); // virtual call to derived
        GISMO_ASSERT( check(), "Something went wrong in assembler initialization");
    }

    /// @brief Intitialize function for, sets data fields
    /// using the pde, a multi-basis and assembler options.
    void initialize(const gsPde<T>           & pde,
                    const gsMultiBasis<T>    & bases,
                    const gsOptionList & opt = defaultOptions() )
    {
        typename gsPde<T>::Ptr _pde = memory::make_shared_not_owned(&pde);
        initialize(_pde,bases,opt);
    }

    void initialize(typename gsPde<T>::Ptr pde,
                    const gsMultiBasis<T>    & bases,
                    const gsOptionList & opt = defaultOptions() )
    {
        m_pde_ptr = pde;
        m_bases.clear();
        m_bases.push_back(bases);
        m_options = opt;
        refresh(); // virtual call to derived
        GISMO_ASSERT( check(), "Something went wrong in assembler initialization");
    }


    /// @brief Intitialize function for, sets data fields
    /// using the pde, a vector of bases and assembler options.
    void initialize(const gsPde<T>           & pde,
                    const gsBasisRefs<T>     & basis,
                    const gsOptionList & opt = defaultOptions() )
    {
        GISMO_ASSERT(pde.domain().nPatches() ==1,
                     "You cannot initialize a multipatch domain with bases on a single patch");
        m_pde_ptr = memory::make_shared_not_owned(&pde);

        m_bases.clear();
        m_bases.reserve(basis.size());
        for(std::size_t c=0;c<basis.size();c++)
            m_bases.push_back(gsMultiBasis<T>(basis[c]));

        m_options = opt;
        refresh(); // virtual call to derived
        GISMO_ASSERT( check(), "Something went wrong in assembler initialization");
    }

    /// @brief checks for consistency and legal values of the stored members.
    bool check();

    /// @brief finishes the assembling of the system matrix,
    /// i.e. calls its .makeCompressed() method.
    void finalize() { m_system.matrix().makeCompressed(); }

    /** @brief Penalty constant for patch \a k, used for Nitsche and
    / Discontinuous Galerkin methods
    */
    T penalty(int k) const
    {
        const int deg = m_bases[0][k].maxDegree();
        return (deg + m_bases[0][k].dim()) * (deg + 1) * T(2.0);
    }
    
    /// @brief Provides an estimation of the number of non-zero matrix
    /// entries per column. This value can be used for sparse matrix
    /// memory allocation
    index_t numColNz() const
    {
        return m_system.numColNz(m_bases.front()[0], m_options);
    }

public:  /* Virtual assembly routines*/

    /// @brief Creates the mappers and setups the sparse system.
    /// to be implemented in derived classes, see scalarProblemGalerkinRefresh()
    /// for a possible implementation
    virtual void refresh();

    /// @brief Main assemble routine, to be implemented in derived classes
    virtual void assemble();

    /// @brief Main non-linear assemble routine with input from
    /// current solution
    virtual void assemble(const gsMultiPatch<T> & curSolution);

    gsOptionList & options() {return m_options;}

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
    /// applies the \a BElementVisitor
    template<class BElementVisitor>
    void push(const bcContainer & BCs)
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

    /// @brief Iterates over all elements of interfaces and
    /// applies the \a InterfaceVisitor
    template<class InterfaceVisitor>
    void pushInterface()
    {
        InterfaceVisitor visitor(*m_pde_ptr);
            
        const gsMultiPatch<T> & mp = m_pde_ptr->domain();
        for ( typename gsMultiPatch<T>::const_iiterator
                  it = mp.iBegin(); it != mp.iEnd(); ++it )
        {
            const boundaryInterface & iFace = //recover master elemen
                ( m_bases[0][it->first() .patch].numElements(it->first() .side() ) <
                  m_bases[0][it->second().patch].numElements(it->second().side() ) ?
                  it->getInverse() : *it );
            
            this->apply(visitor, iFace);
        }
    }


public:  /* Dirichlet degrees of freedom computation */

    /// @brief Triggers computation of the Dirichlet dofs
    /// \param[in] unk the considered unknown
    void computeDirichletDofs(int unk = 0);

    /// @brief the user can manually set the dirichlet Dofs for a given patch and
    /// unknown, based on the Basis coefficients
    /// \param[in] coefMatrix the coefficients of the function
    /// \param[in] unk the consideren unknown
    /// \param[in] patch the patch index
    void setFixedDofs(const gsMatrix<T> & coefMatrix, int unk = 0, int patch = 0);

    /// @brief the user can manually set the dirichlet Dofs for a given patch and
    /// unknown.
    /// \param[in] vals the values of the eliminated dofs.
    /// \param[in] unk the considered unknown
    void setFixedDofVector(gsMatrix<T> & vals, int unk = 0);

    /// Enforce Dirichlet boundary conditions by diagonal penalization
    /// \param[in] unk the considered unknown
    void penalizeDirichletDofs(int unk = 0);

    /// @brief Sets any Dirichlet values to homogeneous (if applicable)
    /// \param[in] unk the considered unknown
    void homogenizeFixedDofs(int unk = 0)
    {
        if(unk==-1)
        {
            for(std::size_t i=0;i<m_ddof.size();++i)
                m_ddof[i].setZero();
        }
        else
            m_ddof[unk].setZero();
    }

    // index_t numFixedDofs(int unk = 0) {return m_dofMappers[unk].boundarySize();}

    /// @brief Returns all the Dirichlet values (if applicable)
    const std::vector<gsMatrix<T> > & allFixedDofs() const { return m_ddof; }

    /// @brief Returns the Dirichlet values for a unknown (if applicable)
    /// \param[in] unk the considered unknown
    const gsMatrix<T> & fixedDofs(int unk=0) const { return m_ddof[unk]; }

    GISMO_DEPRECATED
    const gsMatrix<T> & dirValues(int unk=0) const { return m_ddof[unk]; }//remove

protected:  /* Helpers for Dirichlet degrees of freedom computation */

    /// @brief calculates the values of the eliminated dofs based on Interpolation.
    /// \param[in] mapper the dofMapper for the considered unknown
    /// \param[in] mbasis the multipabasis for the considered unknown
    /// \param[in] unk_ the considered unknown
    void computeDirichletDofsIntpl(const gsDofMapper     & mapper,
                                   const gsMultiBasis<T> & mbasis,
                                   const int unk_ = 0);
    
    /// @brief calculates the values of the eliminated dofs based on L2 Projection.
    /// \param[in] mapper the dofMapper for the considered unknown
    /// \param[in] mbasis the multipabasis for the considered unknown
    /// \param[in] unk_ the considered unknown
    void computeDirichletDofsL2Proj(const gsDofMapper     & mapper,
                                    const gsMultiBasis<T> & mbasis,
                                    const int unk_ = 0);

public:  /* Solution reconstruction */

    /// @brief Construct solution from computed solution vector for a single unknows
    /// \param[in] solVector the solution vector obtained from the linear system
    /// \param[out] result the solution in form of a gsMultiBasis
    /// \param[in] unk the considered unknown
    virtual void constructSolution(const gsMatrix<T>& solVector,
                                   gsMultiPatch<T>& result, int unk = 0) const;



    /// @brief Construct solution from computed solution vector for a set of unknowns.
    /// The result is a vectorfield, where each component is given the corresponding
    /// entry of \par unknowns. This method assumes that the specified unknowns have the
    /// same basis.
    /// \param[in] solVector the solution vector obtained from the linear system
    /// \param[out] result the solution seen as vectorfield in form of a gsMultiBasis
    /// \param[in] unknowns the considered vector of unknowns
    virtual void constructSolution(const gsMatrix<T>& solVector,
                                   gsMultiPatch<T>& result,
                                   const gsVector<index_t>  & unknowns) const;

    gsField<T> constructSolution(const gsMatrix<T>& solVector, int unk = 0) const;

    /// @brief Update solution by adding the computed solution vector
    /// to the current solution specified by \par result. This method assumes that all
    /// unknowns have the same basis.
    /// \param[in] solVector the solution vector obtained from the linear system
    /// \param[out] result the solution in form of a gsMultiBasis, \par solVector is added to the
    ///                   coefficients of result.
    /// \param[in] theta damping factor for update, theta = 1 corresponds to a full step.
    virtual void updateSolution(const gsMatrix<T>& solVector,
                                gsMultiPatch<T>& result, T theta = T(1)) const;

public: // *** Accessors *** 

    /// @brief Return the Pde
    const gsPde<T> & pde() const { return *m_pde_ptr; }

    /// @brief Return the multipatch
    const gsMultiPatch<T> & patches() const { return m_pde_ptr->patches(); }

    /// @brief Return the multi-basis
    const gsMultiBasis<T> & multiBasis(index_t k = 0) const { return m_bases[k]; }

    /// @brief Return the multi-basis. Note: if the basis is altered,
    /// then refresh() should be called
    gsMultiBasis<T> & multiBasis(index_t k = 0) { return m_bases[k]; }

    /// @brief Returns the number of multi-bases
    std::size_t numMultiBasis() const {return m_bases.size(); }

    /// @brief Returns the left-hand global matrix
    const gsSparseMatrix<T> & matrix() const { return m_system.matrix(); }

    /// @brief Returns the left-hand side vector(s)
    /// ( multiple right hand sides possible )
    const gsMatrix<T> & rhs() const { return m_system.rhs(); }
    gsMatrix<T> & rhs() { return m_system.rhs(); }
    
    /// @brief Returns the left-hand global matrix
    const gsSparseSystem<T> & system() const { return m_system; }
    gsSparseSystem<T> & system() { return m_system; }
    
    /// @brief Swaps the actual sparse system with the given one
    void setSparseSystem(gsSparseSystem<T> & sys)
    {
        GISMO_ASSERT(sys.initialized(), "Sparse system must be initialized first");
        m_system.swap(sys);
    }

    /// @brief Returns the number of (free) degrees of freedom
    int numDofs() const
    {
        index_t sum = 0;
        for (index_t c = 0; c!= m_system.numColBlocks(); ++c)
            sum += m_system.colMapper(c).freeSize();
        return sum;
    }

protected:

    /// @brief A prototype of the refresh function for a "standard"
    /// scalar problem.  Creats one mapper based on the set options
    /// and initializes the sparse system (without allocating memory.
    void scalarProblemGalerkinRefresh();

protected:

    /// @brief Generic assembly routine for volume or boundary integrals
    /// \param[in] visitor The visitor for the boundary or volume integral
    /// \param[in] patchIndex The considered patch
    /// \param[in] side The considered boundary side, only necessary for boundary
    /// integrals.
    template<class ElementVisitor>
    void apply(ElementVisitor & visitor,
               int patchIndex = 0,
               boxSide side = boundary::none);

    /// @brief Generic assembly routine for patch-interface integrals
    template<class InterfaceVisitor>
    void apply(InterfaceVisitor & visitor,
               const boundaryInterface & bi);
};

template <class T>
template<class ElementVisitor>
void gsAssembler<T>::apply(ElementVisitor & visitor, 
                           int patchIndex,
                           boxSide side)
{
    //gsDebug<< "Apply to patch "<< patchIndex <<"("<< side <<")\n";
    
    const gsBasisRefs<T> bases(m_bases, patchIndex);
    
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
        visitor.localToGlobal(patchIndex, m_ddof, m_system);
    }
}


template <class T>
template<class InterfaceVisitor>
void gsAssembler<T>::apply(InterfaceVisitor & visitor,
                           const boundaryInterface & bi)
{
    //gsDebug<<"Apply DG on "<< bi <<".\n";
    
    const gsAffineFunction<T> interfaceMap(m_pde_ptr->patches().getMapForInterface(bi));
    
    const int patch1      = bi.first().patch;
    const int patch2      = bi.second().patch;
    const gsBasis<T> & B1 = m_bases[0][patch1];// (!) unknown 0
    const gsBasis<T> & B2 = m_bases[0][patch2];

    gsQuadRule<T> QuRule ; // Quadrature rule
    gsMatrix<T> quNodes1, quNodes2;// Mapped nodes
    gsVector<T> quWeights;         // Mapped weights
    // Evaluation flags for the Geometry map
    unsigned evFlags(0);
    
    const int bSize1      = B1.numElements( bi.first() .side() );
    const int bSize2      = B2.numElements( bi.second().side() );
    const int ratio = bSize1 / bSize2;
    GISMO_ASSERT(bSize1 >= bSize2 && bSize1%bSize2==0,
                 "DG assumes nested interfaces. Got bSize1="<<
                 bSize1<<", bSize2="<<bSize2<<"." );
    
    // Initialize
    visitor.initialize(B1, B2, bi, m_options, QuRule, evFlags);
    
    // Initialize geometry evaluators
    typename gsGeometry<T>::Evaluator geoEval1(
        m_pde_ptr->patches()[patch1].evaluator(evFlags));
    typename gsGeometry<T>::Evaluator geoEval2(
        m_pde_ptr->patches()[patch2].evaluator(evFlags));
    
    // Initialize domain element iterators
    typename gsBasis<T>::domainIter domIt1 = B1.makeDomainIterator( bi.first() .side() );
    typename gsBasis<T>::domainIter domIt2 = B2.makeDomainIterator( bi.second().side() );
    
    //typename gsBasis<T>::domainIter domIt = B2.makeDomainIterator(B2, bi);
    
    int count = 0;
    // iterate over all boundary grid cells on the "left"
    for (; domIt1->good(); domIt1->next() )
    {
        count++;
        // Get the element of the other side in domIter2
        //domIter1->adjacent( bi.orient, *domIter2 );
        
        // Compute the quadrature rule on both sides
        QuRule.mapTo( domIt1->lowerCorner(), domIt1->upperCorner(), quNodes1, quWeights);
        interfaceMap.eval_into(quNodes1,quNodes2);
        
        // Perform required evaluations on the quadrature nodes
        visitor.evaluate(B1, *geoEval1, B2, *geoEval2, quNodes1, quNodes2);
        
        // Assemble on element
        visitor.assemble(*domIt1,*domIt2, *geoEval1, *geoEval2, quWeights);
        
        // Push to global patch matrix (m_rhs is filled in place)
        visitor.localToGlobal(patch1, patch2, m_ddof, m_system);
        
        if ( count % ratio == 0 ) // next master element ?
        {
            domIt2->next();
        }
        
    }
}


} // namespace gismo

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsAssembler.hpp)
#endif
