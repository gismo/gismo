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
                               gsMultiBasis<T>         & bases,
                               gsBoundaryConditions<T> const & bconditions,
                               gsBoundaryConditions<T> const & bconditions2,
                               const gsFunction<T>           & rhs,
                               gsOptionList const & optionList)
                : m_ppde(patches,bconditions,bconditions2,rhs)
        {
            m_options.setInt("DirichletStrategy", dirichlet::elimination);
            m_options.setInt("InterfaceStrategy",  iFace::glue);

            m_options.addReal("mu", "Mu", optionList.getReal("mu"));
            m_options.addReal("wu", "Wu", optionList.getReal("wu"));
            m_options.addInt("l", "Level", optionList.getInt("l"));

            valuePenalty.setZero(patches.nInterfaces());

/*
            // Create interior spline space
            std::vector<gsBasis<T> *> basis_container;
            for (size_t np = 0; np < patches.nPatches(); np++)
            {
                gsTensorBSplineBasis<2, T> * basis_inner = dynamic_cast<gsTensorBSplineBasis<2, T>*>(&bases.basis(np));

                // Construct special space for r = p - 1:
                // The first and the last knot (not 0,1) are repeated +1, e.g.
                // deg 3, r = 2: |||| || | [...] | || ||||
                index_t r = optionList.getInt("discreteRegularity"); // Assume same reg for each direction
                for (index_t uv = 0; uv < 2; uv++)
                {
                    if (basis_inner->degree(uv) - r == 1)
                    {
                        T knot_u = basis_inner->knot(uv,basis_inner->degree(uv)+1);
                        if (knot_u != 1)
                            basis_inner->insertKnot(knot_u,uv,1);

                        if (knot_u != 0.5 && knot_u != 1)
                            basis_inner->insertKnot(1-knot_u,uv,1);
                    }
                }

                //gsInfo << "basis u " << basis_inner.knots(0).asMatrix() << "\n";
                gsInfo << "basis v " << basis_inner->knots(1).asMatrix() << "\n";
                //gsBasis<T> * temp_basis = *dynamic_cast<gsBasis<T>*>(basis_inner);
                //basis_container.push_back(temp_basis);
            }

            //gsMultiBasis<T> mb_temp(basis_container, bases.topology());
            //bases.swap(mb_temp);
            gsDebugVar(bases.basis(0));
*/

            Base::initialize(m_ppde, bases, m_options);
        }

        void refresh();

        void assemble();

        //stability
        void solveStabilityParameter();
        void apply(bhVisitor & visitor, size_t patchIndex, size_t local = -1, boxSide side=boundary::none);
        void apply(gsVisitorInterfaceNitscheBiharmonicStability<T> & visitor,
                   const boundaryInterface & bi, bool local = true);

        gsVector<T> get_valuePenalty() { return valuePenalty; }
        void set_valuePenalty(gsVector<T> stab) { valuePenalty = stab; }

    protected:

        // fixme: add constructor and remove this
        gsBiharmonicPde<T> m_ppde;

        gsVector<T> valuePenalty;

        // Members from gsAssembler
        using Base::m_pde_ptr;
        using Base::m_bases;
        using Base::m_ddof;
        using Base::m_options;
        using Base::m_system;

        // for stabilization
        std::vector<gsMatrix<T> > m_ddof_stab;
        gsSparseSystem<> m_system_stabA, m_system_stabB;
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

        if (m_options.getReal("mu") == -2 && m_options.getInt("l") < 2)
        {
            T h1 = m_bases[0][0].getMinCellLength(); // The same for each patch TODO
            solveStabilityParameter();
            valuePenalty = valuePenalty * h1; // Correction Term
        }

        else if (m_options.getReal("mu") == -2 && m_options.getInt("l") > 1)
            valuePenalty = valuePenalty * 1/16; // Correction Term

/*
        if (m_options.getReal("mu") == -2)
        {
            valuePenalty[0] = 4.39;
            valuePenalty[1] = 4.92;
            valuePenalty[2] = 7.68;
            valuePenalty[3] = 4.64;
            valuePenalty[4] = 6.64;
            valuePenalty[5] = 6.102;
            valuePenalty[6] = 6.58;
            valuePenalty[7] = 5.23;
            valuePenalty[8] = 7.58;
            valuePenalty[9] = 7.09;
            valuePenalty[10] = 7.03;
            valuePenalty[11] = 11.5;
            valuePenalty[12] = 10.57;
            valuePenalty[13] = 5.156;
            valuePenalty[14] = 10.606;
            valuePenalty[15] = 3.104;
            valuePenalty[16] = 18.58;
            valuePenalty[17] = 3.718;
            valuePenalty[18] = 5.94;
            valuePenalty[19] = 5.82;
            valuePenalty[20] = 4.298;
            valuePenalty[21] = 6.418;

            valuePenalty = valuePenalty * 4.1;
        }
*/
        if (m_options.getReal("mu") == -1)
            solveStabilityParameter();



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
/*

        // For the stability parameter
        // We use predefined helper which initializes the system matrix
        // rows and columns using the same test and trial space
        gsDofMapper map(m_bases[0]);
        gsDofMapper map2(m_bases[0]);

        gsMatrix<index_t> act;
        for ( typename gsMultiPatch<T>::const_iiterator
                      it = m_ppde.domain().iBegin(); it != m_ppde.domain().iEnd(); ++it )
        {
            index_t side_1 = it->first().side().index();
            index_t side_2 = it->second().side().index();
            index_t dim_uv = m_bases[0].basis(it->first().patch).component(0).size();
            for (index_t i = 4; i < dim_uv; i++)
            {
                act = m_bases[0].basis(it->first().patch).boundaryOffset(side_1, i); // First
                map2.markBoundary(it->first().patch, act);
                map.markBoundary(it->first().patch, act);
                act = m_bases[0].basis(it->second().patch).boundaryOffset(side_2, i); // First
                map.markBoundary(it->second().patch, act);
                map2.markBoundary(it->second().patch, act);
            }
        }
        map.finalize();
        map.print();

        map2.finalize();
        map2.print();

        // 2. Create the sparse system
        m_system_stab = gsSparseSystem<T>(map);
        m_system = gsSparseSystem<T>(map2);

        GISMO_ASSERT(m_system.initialized(), "Sparse system is not initialized, call refresh()");

        // Reserve sparse system
        const index_t nz = gsAssemblerOptions::numColNz(m_bases[0].basis(0).component(0),2,1,0.333333);
        m_system.reserve(nz, this->pde().numRhs());
        m_system_stab.reserve(nz, this->pde().numRhs());

        // Compute the Dirichlet Degrees of freedom (if needed by m_options)
        //Base::computeDirichletDofs();
        m_ddof.resize(m_system.numUnknowns());
        m_ddof[0].setZero(m_system.colMapper(0).boundarySize(), m_system.unkSize(0) * m_pde_ptr->numRhs() );

        // Assemble volume integrals
        Base::template push<bhVisitor >();

        m_system_stab.setZero();
        pushInterface();

        //Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXf> ges(AA, BB);
        Eigen::MatrixXd A = m_system_stab.matrix().toDense();
        Eigen::MatrixXf AA = A.cast <float> ();

        Eigen::MatrixXd B = m_system.matrix().toDense();
        Eigen::MatrixXf BB = B.cast <float> ();
        Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> ges(A, B);
        //gsInfo << "The (complex) numerators of the generalzied eigenvalues are: " << ges.alphas().transpose() << "\n";
        //gsInfo << "The (real) denominatore of the generalzied eigenvalues are: " << ges.betas().transpose() << "\n";
        gsInfo << "The (complex) generalzied eigenvalues are (alphas./beta): " << ges.eigenvalues().transpose() << "\n";
*/

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
    void gsBiharmonicNitscheAssembler<T,bhVisitor>::solveStabilityParameter()
    {


        const gsMultiPatch<T> & mp = m_pde_ptr->domain();

        index_t i = 0;
        for ( typename gsMultiPatch<T>::const_iiterator
                      it = mp.iBegin(); it != mp.iEnd(); ++it, ++i)
        {
            const boundaryInterface & iFace = *it;
            //recover master elemen
/*                    ( m_bases[0][it->first() .patch].numElements(it->first() .side() ) <
                      m_bases[0][it->second().patch].numElements(it->second().side() ) ?
                      it->getInverse() : *it );*/


            gsVector<index_t> sz(2);
            sz[0] = m_bases[0][iFace.first().patch].size();
            sz[1] = m_bases[0][iFace.second().patch].size();

            gsDofMapper map(sz);
            gsDofMapper map2(sz);

            index_t side_11 = iFace.first().side().index();
            index_t side_22 = iFace.second().side().index();

            gsMatrix<index_t> matched_dofs1, matched_dofs2;
            matched_dofs1 = m_bases[0][iFace.first().patch].boundaryOffset(side_11, 0); // First
            matched_dofs2 = m_bases[0][iFace.second().patch].boundaryOffset(side_22, 0); // Second

            map.matchDofs(0, matched_dofs1, 1, matched_dofs2);
            map2.matchDofs(0, matched_dofs1, 1, matched_dofs2);

            gsMatrix<index_t> act;
            index_t dim_uv = m_bases[0][iFace.first().patch].component(0).size();
            for (index_t u_i = 2; u_i < dim_uv; u_i++)
            {
                act = m_bases[0][iFace.first().patch].boundaryOffset(side_11, u_i); // First
                map.markBoundary(0, act);
                map2.markBoundary(0, act);
                act = m_bases[0][iFace.second().patch].boundaryOffset(side_22, u_i); // Second
                map.markBoundary(1, act);
                map2.markBoundary(1, act);
            }



            map2.finalize();
            map.finalize();

            m_system_stabA = gsSparseSystem<T>(map);
            m_system_stabB = gsSparseSystem<T>(map2);

            const index_t nz = dim_uv;
            m_system_stabA.reserve(nz, this->pde().numRhs());
            m_system_stabB.reserve(nz, this->pde().numRhs());
            // Compute the Dirichlet Degrees of freedom (if needed by m_options)
            m_ddof_stab.resize(m_system_stabB.numUnknowns());
            m_ddof_stab[0].setZero(m_system_stabB.colMapper(0).boundarySize(), m_system_stabB.unkSize(0) * m_pde_ptr->numRhs() );

            // Assemble volume integrals
            m_system_stabB.setZero();
            for (size_t np = 0; np < 2; ++np)
            {
                bhVisitor visitor2(*m_pde_ptr);
                //Assemble (fill m_matrix and m_rhs) on patch np
                apply(visitor2, np == 0 ? iFace.first().patch : iFace.second().patch, np);
            }

            //Base::template push<bhVisitor >();

            m_system_stabA.setZero();
            gsVisitorInterfaceNitscheBiharmonicStability<T> visitor(*m_pde_ptr);
            apply(visitor, iFace);

            //Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXf> ges(AA, BB);
            Eigen::MatrixXd A = m_system_stabA.matrix().toDense().cast<double>();
            //Eigen::MatrixXf AA = A.cast <float> ();

            Eigen::MatrixXd B = m_system_stabB.matrix().toDense().cast<double>();
            //Eigen::MatrixXf BB = B.cast <float> ();
            Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> ges(A, B);
            //gsInfo << "The (complex) numerators of the generalzied eigenvalues are: " << ges.alphas().transpose() << "\n";
            //gsInfo << "The (real) denominatore of the generalzied eigenvalues are: " << ges.betas().transpose() << "\n";
            //gsInfo << "Max: " << ges.eigenvalues().transpose() << "\n";
            //gsInfo << "Max 2: " << ges.eigenvalues().array().abs().maxCoeff() << "\n";

            //Eigen::GeneralizedEigenSolver<Eigen::MatrixXd> ges2;
            //ges2.compute(A, B);
            //gsInfo << "The (complex) generalzied eigenvalues are (alphas./beta): " << ges.eigenvalues().transpose() << "\n";

            valuePenalty[i] = ges.eigenvalues().array().maxCoeff() * m_options.getReal("wu"); // * 8.1;
        }

        //T max_penalty = valuePenalty.maxCoeff();
        //valuePenalty.setConstant(max_penalty);
        gsInfo << "Max value: " << valuePenalty.transpose() << "\n";
/*
        {
            gsDofMapper map(m_bases[0]);
            gsDofMapper map2(m_bases[0]);
            for ( typename gsMultiPatch<T>::const_iiterator
                          it = mp.iBegin(); it != mp.iEnd(); ++it) {
                const boundaryInterface &iFace = //recover master elemen
                        (m_bases[0][it->first().patch].numElements(it->first().side()) <
                         m_bases[0][it->second().patch].numElements(it->second().side()) ?
                         it->getInverse() : *it);

                gsMatrix<index_t> act;
                index_t side_11 = iFace.first().side().index();
                index_t side_22 = iFace.second().side().index();
                index_t dim_uv = m_bases[0][iFace.first().patch].component(0).size();
                for (index_t u_i = 5; u_i < dim_uv; u_i++) {
                    act = m_bases[0][iFace.first().patch].boundaryOffset(side_11, u_i); // First
                    map.markBoundary(iFace.first().patch, act);
                    map2.markBoundary(iFace.first().patch, act);
                    act = m_bases[0][iFace.second().patch].boundaryOffset(side_22, u_i); // Second
                    map.markBoundary(iFace.second().patch, act);
                    map2.markBoundary(iFace.second().patch, act);
                }

                gsMatrix<index_t> matched_dofs1, matched_dofs2;
                matched_dofs1 = m_bases[0][iFace.first().patch].boundaryOffset(side_11, 0); // First
                matched_dofs2 = m_bases[0][iFace.second().patch].boundaryOffset(side_22, 0); // Second

                map.matchDofs(iFace.first().patch, matched_dofs1, iFace.second().patch, matched_dofs2);
                map2.matchDofs(iFace.first().patch, matched_dofs1, iFace.second().patch, matched_dofs2);
            }

            map2.finalize();
            map.finalize();

            m_system_stabA = gsSparseSystem<T>(map);
            m_system_stabB = gsSparseSystem<T>(map2);

            const index_t nz = m_bases[0][0].component(0).size();;
            m_system_stabA.reserve(nz, this->pde().numRhs());
            m_system_stabB.reserve(nz, this->pde().numRhs());
            // Compute the Dirichlet Degrees of freedom (if needed by m_options)
            m_ddof_stab.resize(m_system_stabB.numUnknowns());
            m_ddof_stab[0].setZero(m_system_stabB.colMapper(0).boundarySize(), m_system_stabB.unkSize(0) * m_pde_ptr->numRhs() );

            // Assemble volume integrals
            m_system_stabB.setZero();
            for (size_t np = 0; np < mp.nPatches(); ++np)
            {
                bhVisitor visitor2(*m_pde_ptr);
                //Assemble (fill m_matrix and m_rhs) on patch np
                apply(visitor2, np);
            }

            //Base::template push<bhVisitor >();

            m_system_stabA.setZero();
            typename gsMultiPatch<T>::const_iiterator
                    it = mp.iBegin();
            //for ( typename gsMultiPatch<T>::const_iiterator
            //              it = mp.iBegin(); it != mp.iEnd(); ++it, ++i) {
                const boundaryInterface &iFace = //recover master elemen
                        (m_bases[0][it->first().patch].numElements(it->first().side()) <
                         m_bases[0][it->second().patch].numElements(it->second().side()) ?
                         it->getInverse() : *it);
                gsVisitorInterfaceNitscheBiharmonicStability<T> visitor(*m_pde_ptr);
                apply(visitor, iFace, false);
           // }

            //Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXf> ges(AA, BB);
            Eigen::MatrixXd A = m_system_stabA.matrix().toDense();
            //Eigen::MatrixXf AA = A.cast <float> ();

            Eigen::MatrixXd B = m_system_stabB.matrix().toDense();
            //Eigen::MatrixXf BB = B.cast <float> ();
            Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> ges(A, B);

            //Eigen::GeneralizedEigenSolver<Eigen::MatrixXd> ges2;
            //ges2.compute(A, B);
            gsInfo << "The (complex) generalzied eigenvalues are (alphas./beta): " << ges.eigenvalues().array().maxCoeff()*8 << "\n";

            //valuePenalty[i] = ges.eigenvalues().array().maxCoeff()*8;


        }
*/
    }

    template <class T, class bhVisitor>
    void gsBiharmonicNitscheAssembler<T,bhVisitor>::apply(bhVisitor & visitor,
                               size_t patchIndex,
                               size_t local,
                               boxSide side)
    {
        //gsDebug<< "Apply to patch "<< patchIndex <<"("<< side <<")\n";

        const gsBasis<T> & bases = m_bases[0][patchIndex]; // only for same basis for each patch

#pragma omp parallel
        {
            gsQuadRule<T> quRule ; // Quadrature rule
            gsMatrix<T> quNodes  ; // Temp variable for mapped nodes
            gsVector<T> quWeights; // Temp variable for mapped weights

            bhVisitor
#ifdef _OPENMP
            // Create thread-private visitor
    visitor_(visitor);
    const int tid = omp_get_thread_num();
    const int nt  = omp_get_num_threads();
#else
                    &visitor_ = visitor;
#endif

            // Initialize reference quadrature rule and visitor data
            visitor_.initialize(bases, patchIndex, m_options, quRule);

            const gsGeometry<T> & patch = m_pde_ptr->domain()[patchIndex];

            // Initialize domain element iterator -- using unknown 0
            typename gsBasis<T>::domainIter domIt = bases.makeDomainIterator(side);

            // Start iteration over elements
#ifdef _OPENMP
            for ( domIt->next(tid); domIt->good(); domIt->next(nt) )
#else
            for (; domIt->good(); domIt->next() )
#endif
            {
                // Map the Quadrature rule to the element
                quRule.mapTo( domIt->lowerCorner(), domIt->upperCorner(), quNodes, quWeights );

                // Perform required evaluations on the quadrature nodes
                visitor_.evaluate(bases, patch, quNodes);

                // Assemble on element
                visitor_.assemble(*domIt, quWeights);

                // Push to global matrix and right-hand side vector
#pragma omp critical(localToGlobal)
                visitor_.localToGlobal(local == -1 ? patchIndex : local, m_ddof_stab, m_system_stabB); // omp_locks inside
            }
        }//omp parallel

    }

    template <class T, class bhVisitor>
    void gsBiharmonicNitscheAssembler<T,bhVisitor>::apply(gsVisitorInterfaceNitscheBiharmonicStability<T> & visitor,
                               const boundaryInterface & bi, bool local)
    {
        gsRemapInterface<T> interfaceMap(m_pde_ptr->domain(), m_bases[0], bi);
        const index_t patchIndex1      = local ? 0 : bi.first().patch;
        const index_t patchIndex2      = local ? 1 : bi.second().patch;
        const gsBasis<T> & B1 = m_bases[0][bi.first().patch];// (!) unknown 0
        const gsBasis<T> & B2 = m_bases[0][bi.second().patch];

        gsQuadRule<T> quRule ; // Quadrature rule
        gsMatrix<T> quNodes1, quNodes2;// Mapped nodes
        gsVector<T> quWeights;         // Mapped weights

        // Initialize
        visitor.initialize(B1, B2, bi, m_options, quRule);

        const gsGeometry<T> & patch1 = m_pde_ptr->domain()[bi.first().patch];
        const gsGeometry<T> & patch2 = m_pde_ptr->domain()[bi.second().patch];

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
            visitor.localToGlobal(patchIndex1, patchIndex2, m_ddof_stab, m_system_stabA);
        }

    }
} // namespace gismo



