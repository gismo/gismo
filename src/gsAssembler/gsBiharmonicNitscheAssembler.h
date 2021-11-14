/** @file gsBiharmonicNitscheAssembler.h

    @brief Provides assembler for a homogenius Biharmonic equation.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): J. Sogn
*/


#pragma once

#include <gsAssembler/gsAssembler.h>

#include <gsPde/gsBiharmonicPde.h>

#include <gsAssembler/gsVisitorBiharmonic.h>
#include <gsAssembler/gsVisitorNeumann.h>
#include <gsAssembler/gsVisitorNeumannBiharmonic.h>
#include <gsAssembler/gsVisitorInterfaceNitscheBiharmonic.h>
#include <gsAssembler/gsVisitorInterfaceNitscheBiharmonicStability.h>

namespace gismo
{

/** @brief
    Implementation of a homogeneous Biharmonic Assembler.

    It sets up an assembler and assembles the system patch wise and
    combines the patch-local stiffness matrices into a global system.
    Dirichlet boundary can only be enforced strongly (i.e Nitsche is
    not implemented).
*/
    template <class T, class bhVisitor = gsVisitorBiharmonic<T> >
    class gsBiharmonicNitscheAssembler : public gsAssembler<T>
    {
    public:
        typedef gsAssembler<T> Base;

    public:
/** @brief
    Constructor of the assembler object.

    \param[in] patches is a gsMultiPatch object describing the geometry.
    \param[in] bases a multi-basis that contains patch-wise bases
    \param[in] bconditions  is a gsBoundaryConditions object that holds boundary conditions on the form:
    \f[ \text{Dirichlet: } u = g \text{ on } \Gamma, \text{ and Neumann: } \nabla \Delta u \cdot \mathbf{n} = h \text{ on } \Gamma\f]
    \param[in] bconditions2 is a gsBoundaryConditions object that holds Neumann boundary conditions on the form:
    \f[\text{Neumann: } \nabla \Delta u \cdot \mathbf{n} = g\, \rightarrow \,(g,\nabla v \cdot \mathbf{n})_\Gamma, \f] where \f$ g \f$ is the Neumann data,
    \f$ v \f$ is the test function and \f$ \Gamma\f$ is the boundary side.
    \param[in] rhs is the right-hand side of the Biharmonic equation, \f$\mathbf{f}\f$.
    \param[in] dirStrategy option for the treatment of Dirichlet boundary in the \em bconditions object.
    \param[in] intStrategy option for the treatment of patch interfaces
*/
        gsBiharmonicNitscheAssembler( gsMultiPatch<T> const         & patches,
                               gsMultiBasis<T> const         & bases,
                               gsBoundaryConditions<T> const & bconditions,
                               gsBoundaryConditions<T> const & bconditions2,
                               const gsFunction<T>           & rhs,
                               gsOptionList const & optionList)
                : m_ppde(patches,bconditions,bconditions2,rhs)
        {
            m_options.setInt("DirichletStrategy", dirichlet::elimination);
            m_options.setInt("InterfaceStrategy",  iFace::glue);

            m_options.addReal("mu", "Mu", optionList.getReal("mu"));

            valuePenalty.setZero(patches.nInterfaces());

            Base::initialize(m_ppde, bases, m_options);
        }

        void refresh();

        void assemble();

        //stability
        void pushInterface();
        void apply(gsVisitorInterfaceNitscheBiharmonicStability<T> & visitor,
                   const boundaryInterface & bi);

        gsVector<> get_valuePenalty() { return valuePenalty; }

    protected:

        // fixme: add constructor and remove this
        gsBiharmonicPde<T> m_ppde;

        gsVector<> valuePenalty;

        // Members from gsAssembler
        using Base::m_pde_ptr;
        using Base::m_bases;
        using Base::m_ddof;
        using Base::m_options;
        using Base::m_system;

        // for stabilization
        gsSparseSystem<> m_system_stab;
    };

    template <class T, class bhVisitor>
    void gsBiharmonicNitscheAssembler<T,bhVisitor>::refresh()
    {
        /*
        // We use predefined helper which initializes the system matrix
        // rows and columns using the same test and trial space
        gsDofMapper map(m_bases[0]);

        gsMatrix<index_t> act;
        for (typename gsBoundaryConditions<T>::bcContainer::const_iterator it
                = m_ppde.bcFirstKind().dirichletSides().begin(); it!= m_ppde.bcFirstKind().dirichletSides().end(); ++it)
        {
            act = m_bases[0].basis(it->patch()).boundaryOffset(it->side().index(), 0); // First
            map.markBoundary(it->patch(), act);
        }

        for (typename gsBoundaryConditions<T>::bcContainer::const_iterator it
                = m_ppde.bcSecondKind().neumannSides().begin(); it!= m_ppde.bcSecondKind().neumannSides().end(); ++it)
        {
            act = m_bases[0].basis(it->patch()).boundaryOffset(it->side().index(), 1); // Second
            // without the first and the last (already marked from dirichlet boundary)
            map.markBoundary(it->patch(), act.block(1,0,act.rows()-2,1));
        }

        map.finalize();

        // 2. Create the sparse system
        m_system = gsSparseSystem<T>(map);
         */
        Base::scalarProblemGalerkinRefresh();
    }

    template <class T, class bhVisitor>
    void gsBiharmonicNitscheAssembler<T,bhVisitor>::assemble()
    {
        GISMO_ASSERT(m_system.initialized(), "Sparse system is not initialized, call refresh()");

        // Reserve sparse system
        const index_t nz = gsAssemblerOptions::numColNz(m_bases[0][0],2,1,0.333333);
        m_system.reserve(nz, this->pde().numRhs());

        // Compute the Dirichlet Degrees of freedom (if needed by m_options)
        Base::computeDirichletDofs();

        // Assemble volume integrals
        Base::template push<bhVisitor >();

        // Neumann conditions of first kind
        //Base::template push<gsVisitorNeumann<T> >(
        //        m_ppde.bcFirstKind().neumannSides() );

        // For the stability parameter
        //gsSparseMatrix<> matrix_B = m_system.matrix();
        //m_system_stab = m_system;
        //m_system_stab.setZero();
        //pushInterface();

        // Neumann conditions of second kind // TODO Rename to Laplace
        Base::template push<gsVisitorNeumannBiharmonic<T> >(
                m_ppde.bcSecondKind().laplaceSides() );

        // Add interface integrals
        Base::template pushInterface<gsVisitorInterfaceNitscheBiharmonic<T>>( valuePenalty );

        if ( m_options.getInt("InterfaceStrategy") == iFace::dg )
            gsWarn <<"DG option ignored.\n";

        /*
        // If requested, force Dirichlet boundary conditions by Nitsche's method
        this->template push<gsVisitorNitscheBiharmonic<T> >(
        m_ppde.bcSecondKind().dirichletSides() );
        */

        // Assembly is done, compress the matrix
        Base::finalize();

        //gsInfo << Base::m_system.matrix().toDense() << "\n";
        //gsInfo << Base::m_system.rhs() << "\n";
    }

    template <class T, class bhVisitor>
    void gsBiharmonicNitscheAssembler<T,bhVisitor>::pushInterface()
    {
        gsVisitorInterfaceNitscheBiharmonicStability<T> visitor(*m_pde_ptr);

        const gsMultiPatch<T> & mp = m_pde_ptr->domain();
        for ( typename gsMultiPatch<T>::const_iiterator
                      it = mp.iBegin(); it != mp.iEnd(); ++it )
        {
            const boundaryInterface & iFace = //recover master elemen
                    ( m_bases[0][it->first() .patch].numElements(it->first() .side() ) <
                      m_bases[0][it->second().patch].numElements(it->second().side() ) ?
                      it->getInverse() : *it );

            apply(visitor, iFace);
        }
    }

    template <class T, class bhVisitor>
    void gsBiharmonicNitscheAssembler<T,bhVisitor>::apply(gsVisitorInterfaceNitscheBiharmonicStability<T> & visitor,
                               const boundaryInterface & bi)
    {
        gsRemapInterface<T> interfaceMap(m_pde_ptr->patches(), m_bases[0], bi);

        const index_t patchIndex1      = bi.first().patch;
        const index_t patchIndex2      = bi.second().patch;
        const gsBasis<T> & B1 = m_bases[0][patchIndex1];// (!) unknown 0
        const gsBasis<T> & B2 = m_bases[0][patchIndex2];

        gsQuadRule<T> quRule ; // Quadrature rule
        gsMatrix<T> quNodes1, quNodes2;// Mapped nodes
        gsVector<T> quWeights;         // Mapped weights

        // Initialize
        visitor.initialize(B1, B2, bi, m_options, quRule);

        const gsGeometry<T> & patch1 = m_pde_ptr->patches()[patchIndex1];
        const gsGeometry<T> & patch2 = m_pde_ptr->patches()[patchIndex2];

        // Initialize domain element iterators
        typename gsBasis<T>::domainIter domIt = interfaceMap.makeDomainIterator();
        //int count = 0;

        // iterate over all boundary grid cells on the "left"
        for (; domIt->good(); domIt->next() )
        {
            //count++;

            // Compute the quadrature rule on both sides
            quRule.mapTo( domIt->lowerCorner(), domIt->upperCorner(), quNodes1, quWeights);
            interfaceMap.eval_into(quNodes1,quNodes2);

            // Perform required evaluations on the quadrature nodes
            visitor.evaluate(B1, patch1, B2, patch2, quNodes1, quNodes2);

            // Assemble on element
            visitor.assemble(*domIt,*domIt, quWeights);

            // Push to global patch matrix (m_rhs is filled in place)
            visitor.localToGlobal(patchIndex1, patchIndex2, m_ddof, m_system_stab);
        }

    }
} // namespace gismo



