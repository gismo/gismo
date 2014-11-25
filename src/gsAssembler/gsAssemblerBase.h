/** @file gsAssemblerBase.h

    @brief Provides generic assembler routines

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris
*/

#pragma once

#include <gsCore/gsBasisRefs.h>
#include <gsCore/gsStdVectorRef.h>

namespace gismo
{

/** @brief The assembler class provides to generic routines for volume
 * and boundary integrals that are used for for matrix and rhs
 * generation
*/
template <class T>
class gsAssemblerBase
{
private:
    typedef gsStdVectorRef<gsDofMapper> gsDofMappers;

public:

    /// Constructor using a multipatch domain
    gsAssemblerBase(const gsMultiPatch<T> & patches) :
    m_patches(patches)
    { }

    ~gsAssemblerBase()
    { }

    /// Generic assembly routine for volume or boundary integrals
    template<class ElementVisitor>
    void apply(ElementVisitor & visitor, 
               int patchIndex = 0, 
               boundary::side side = boundary::none)
    {
        //gsDebug<< "Apply to patch "<< patchIndex <<"("<< side <<")\n";

        gsSparseMatrix<T> & patchMatrix = m_matrix;// put as argument?
        
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
            visitor.evaluate(bases, /* *domIt,*/ *geoEval, quNodes);//
            
            // Assemble on element
            visitor.assemble(*domIt, *geoEval, quWeights);
            
            // Push to global patch matrix
            // Note: m_rhs is filled in place
            visitor.localToGlobal(mappers, m_ddof, patchIndex, patchMatrix, m_rhs);
        }
    }


    /// Generic assembly routine for patch-interface integrals
    template<class InterfaceVisitor>
    void apply(InterfaceVisitor & visitor, 
               const boundaryInterface & bi)
    {
        //gsDebug<<"Apply DG on "<< bi <<".\n";

        const gsDofMappers mappers(m_dofMappers);
        const int patch1      = bi[0].patch;
        const int patch2      = bi[1].patch;
        const gsBasis<T> & B1 = m_bases[0][patch1];// (!) unknown 0
        const gsBasis<T> & B2 = m_bases[0][patch2];
        const int bSize1      = B1.numElements();
        const int bSize2      = B2.numElements();
        GISMO_ASSERT(bSize1 >= bSize2 && bSize1%bSize2==0,
                     "DG assumes nested interfaces.");
        const boundary::side & side1 = bi[0].side;
        const boundary::side & side2 = bi[1].side;
        
        gsQuadRule<T> QuRule;         // Reference Quadrature rule
        gsMatrix<T> quNodes1, quNodes2;// Mapped nodes
        gsVector<T> quWeights;         // Mapped weights
        // Evaluation flags for the Geometry map
        unsigned evFlags(0);
        
        // Initialize 
        visitor.initialize(B1, B2, QuRule, evFlags);
        
        // Initialize geometry evaluators
        typename gsGeometry<T>::Evaluator geoEval1( 
            m_patches[patch1].evaluator(evFlags));
        typename gsGeometry<T>::Evaluator geoEval2( 
            m_patches[patch2].evaluator(evFlags));

        // Initialize domain element iterators
        typename gsBasis<T>::domainIter domIt1 = B1.makeDomainIterator(side1);
        typename gsBasis<T>::domainIter domIt2 = B2.makeDomainIterator(side2);
        
        //const index_t numNodes = QuRule.numNodes();

        const int ratio = bSize1 / bSize2;
    
        // const int dir1             = direction(side1);
        // const int dir2             = direction(side2);
        const bool par2               = parameter(side2);
        //GISMO_ASSERT( B1.component(!dir1).size() == B2.component(!dir2).size(), 
        //              "DG method not implemented yet for non matching interfaces");

        int count = 0;
        // iterate over all boundary grid cells on the "left"
        for (; domIt1->good(); domIt1->next() )
        {
            count++;
            // Get the element of the other side in domIter2
            //domIter1->adjacent( bi.orient, *domIter2 );
            
            // Compute the quadrature rule on both sides
            QuRule.mapTo( domIt1->lowerCorner(), domIt1->upperCorner(), quNodes1, quWeights);
            mapGaussNodes( quNodes1, side1, side2, ( par2 ? 1.0 : 0.0 ), quNodes2 );// todo: move to visitor!

            // Perform required evaluations on the quadrature nodes            
            visitor.evaluate(B1, *geoEval1, B2, *geoEval2, quNodes1, quNodes2);

            // Assemble on element
            visitor.assemble(*domIt1, *geoEval1, *geoEval2, quWeights);
            
            // Push to global patch matrix (m_rhs is filled in place)
            visitor.localToGlobal(mappers, patch1, patch2, m_matrix, m_rhs);
            
            if ( count % ratio == 0 ) // next master element ?
                domIt2->next();
        }
    }


public:

    /// Return the multipatch.
    const gsMultiPatch<T> & patches() const { return m_patches; }

    /// Return the multi-basis
    const gsMultiBasis<T> & multiBasis(index_t k) const { return m_bases[k]; }

    /// Return the DOF mapper for unknown \em i.
    const gsDofMapper& dofMapper(unsigned i = 0) const     { return m_dofMappers[i]; }

    /// @brief Returns the left-hand global matrix
    const gsSparseMatrix<T> & matrix() const { return m_matrix; }

    /// @brief Returns the left-hand side vector(s)
    /// ( multiple right hand sides possible )
    const gsMatrix<T> & rhs() const { return m_rhs; }

    /// @brief Returns the number of (free) degrees of freedom
    int numDofs() const { return m_dofs; }
    
protected:
    
    static void mapGaussNodes(const gsMatrix<T> & nodes1, 
                              const  boundary::side & side1,
                              const  boundary::side & side2,
                              T fixedParam,
                              gsMatrix<T> & nodes2 );    
protected:

    /// The multipatch domain
    gsMultiPatch<T> m_patches;

    /// The discretization bases corresponding to \a m_patches and to
    /// the number of solution fields that are to be computed
    /// m_bases[i]: The multi-basis for unknown i
    std::vector< gsMultiBasis<T> > m_bases;
    
    /// The Dof mapper is used to map patch-local DoFs to the global DoFs
    /// One for each unknown, one for each patch
    /// m_dofMappers[i]: DoF Mapper for unknown i
    std::vector<gsDofMapper>  m_dofMappers;
    
    // Dirichlet DoF fixed values (if applicable)
    gsMatrix<T> m_ddof;

    // *** Outputs *** 
    
    /// Global matrix
    gsSparseMatrix<T> m_matrix;

    /// Right-hand side ( multiple right hand sides possible )
    gsMatrix<T>       m_rhs;

    // *** Information *** 

    /// number of degrees of freedom (excluding eliminated etc)
    int m_dofs;

};



template <class T>
void gsAssemblerBase<T>::mapGaussNodes(const gsMatrix<T> & nodes1, 
                                    const  boundary::side & side1,
                                    const  boundary::side & side2,
                                    T fixedParam,
                                    gsMatrix<T> & nodes2 )
{
    const int dir1 = direction(side1);
    const int dir2 = direction(side2);
    const int d    = nodes1.rows();
    nodes2.resize( d, nodes1.cols() );
    
    if ( dir1 == dir2 )
    {
        nodes2.row(dir2).setConstant(fixedParam);
        nodes2.topRows(dir2)        = nodes1.topRows( dir2 );
        nodes2.bottomRows(d-dir2-1) = nodes1.bottomRows( d-dir2-1 );
    }
    else
    {
        GISMO_ASSERT( nodes1.rows() == 2, "Implemented for 2D");
            nodes2.row(dir2).setConstant(fixedParam);
            nodes2.row( !dir2 )  = nodes1.row( dir2 );
    }
}



} // namespace gismo

