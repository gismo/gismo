/** @file gsBiharmonicAssembler.h

    @brief Provides assembler for a homogenius Biharmonic equation.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): J. Sogn & P. Weinmüller
*/


#pragma once

#include <gsAssembler/gsAssembler.h>

#include <gsPde/gsBiharmonicPde.h>

#include <gsAssembler/gsVisitorBiharmonic.h>
#include <gsAssembler/gsVisitorNeumann.h>
#include <gsAssembler/gsVisitorLaplaceBoundaryBiharmonic.h>
//#include <gsAssembler/gsVisitorNitscheBiharmonic.h>

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
class gsBiharmonicAssembler : public gsAssembler<T>
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
    gsBiharmonicAssembler( gsMultiPatch<T> const         & patches,
                           gsMultiBasis<T> const         & bases,
                           gsBoundaryConditions<T> const & bconditions,
                           gsBoundaryConditions<T> const & bconditions2,
                           const gsFunction<T>           & rhs,
                           dirichlet::strategy           dirStrategy,
                           iFace::strategy               intStrategy = iFace::glue)
    : m_ppde(patches,bconditions,bconditions2,rhs)
    {
        m_options.setInt("DirichletStrategy", dirStrategy);
        m_options.setInt("InterfaceStrategy", intStrategy);

        Base::initialize(m_ppde, bases, m_options);
    }

    void refresh();
    
    void assemble();

    void computeDirichletAndNeumannDofs();

protected:

    // fixme: add constructor and remove this
    gsBiharmonicPde<T> m_ppde;

    // Members from gsAssembler
    using Base::m_pde_ptr;
    using Base::m_bases;
    using Base::m_ddof;
    using Base::m_options;
    using Base::m_system;
};

template <class T, class bhVisitor>
void gsBiharmonicAssembler<T,bhVisitor>::refresh()
{
    // We use predefined helper which initializes the system matrix
    // rows and columns using the same test and trial space

    // Pascal
    gsDofMapper map(m_bases[0]);

    gsMatrix<unsigned> act;
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
    // END
}

template <class T, class bhVisitor>
void gsBiharmonicAssembler<T,bhVisitor>::assemble()
{
    GISMO_ASSERT(m_system.initialized(), "Sparse system is not initialized, call refresh()");

    // Reserve sparse system
    const index_t nz = gsAssemblerOptions::numColNz(m_bases[0][0],2,1,0.333333);
    m_system.reserve(nz, this->pde().numRhs());

    // Compute the Dirichlet Degrees of freedom (if needed by m_options)
    if (m_ppde.bcSecondKind().laplaceSides().size() != 0)
        Base::computeDirichletDofs();
    else
        computeDirichletAndNeumannDofs();

    // Assemble volume integrals
    Base::template push<bhVisitor >();

    // Neuman conditions of first kind
    Base::template push<gsVisitorNeumann<T> >(
        m_ppde.bcFirstKind().neumannSides() );

    // Laplace conditions of second kind
    //Base::template push<gsVisitorLaplaceBoundaryBiharmonic<T> >(
    //    m_ppde.bcSecondKind().laplaceSides() );

    if ( m_options.getInt("InterfaceStrategy") == iFace::dg )
        gsWarn <<"DG option ignored.\n";

    /*
    // If requested, force Dirichlet boundary conditions by Nitsche's method
    this->template push<gsVisitorNitscheBiharmonic<T> >(
    m_ppde.bcSecondKind().dirichletSides() );
    */
    
    // Assembly is done, compress the matrix
    Base::finalize();
}


// Pascal
template <class T, class bhVisitor>
void gsBiharmonicAssembler<T,bhVisitor>::computeDirichletAndNeumannDofs()
{

    const gsMultiBasis<T> & mbasis = m_bases[m_system.colBasis(0)];
    const gsDofMapper & mapper =
        dirichlet::elimination == m_options.getInt("DirichletStrategy") ?
        m_system.colMapper(0) :
        mbasis.getMapper(dirichlet::elimination,
                         static_cast<iFace::strategy>(m_options.getInt("InterfaceStrategy")),
                         m_pde_ptr->bc(), 0);

    if(m_ddof.size()==0)
        m_ddof.resize(m_system.numUnknowns());

    m_ddof[0].resize( mapper.boundarySize(), m_system.unkSize(0)*m_system.rhs().cols());  //m_pde_ptr->numRhs() );

    // Set up matrix, right-hand-side and solution vector/matrix for
    // the L2-projection
    gsSparseEntries<T> projMatEntries;
    gsMatrix<T>        globProjRhs;
    globProjRhs.setZero( mapper.boundarySize(), m_system.unkSize(0)*m_system.rhs().cols() );

    // Temporaries
    gsVector<T> quWeights;

    gsMatrix<T> rhsVals, rhsVals2;
    gsMatrix<unsigned> globIdxAct;
    gsMatrix<T> basisVals, basisGrads, physBasisGrad;

    gsVector<T> unormal;

    real_t lambda = 1e-2;

    gsMapData<T> md(NEED_VALUE | NEED_MEASURE | NEED_GRAD_TRANSFORM | NEED_JACOBIAN);


    typename gsBoundaryConditions<T>::const_iterator
        iter_dir = m_pde_ptr->bc().dirichletBegin();

    for ( typename gsBoundaryConditions<T>::const_iterator
              iter = m_ppde.bcSecondKind().neumannBegin();
          iter != m_ppde.bcSecondKind().neumannEnd(); ++iter )
    {
        if (iter->isHomogeneous() )
            continue;

        GISMO_ASSERT(iter_dir->function()->targetDim() == m_system.unkSize(0)*m_system.rhs().cols(),
                     "Given Dirichlet boundary function does not match problem dimension. "
                         << iter_dir->function()->targetDim()<<" != "<<m_system.unkSize(0)<<"\n");

        const int unk = iter_dir->unknown();

        const int patchIdx   = iter_dir->patch();
        const gsBasis<T> & basis = (m_bases[unk])[patchIdx];

        GISMO_ASSERT(iter_dir->patch() == patchIdx && iter_dir->side().index() == iter->side().index(),
                     "Given Dirichlet boundary edge does not match the neumann edge."
                         <<iter_dir->patch()<<" != "<<patchIdx<<" and "
                         <<iter_dir->side().index()<<" != "<<iter->side().index()<<"\n");

        const gsGeometry<T> & patch = m_pde_ptr->patches()[patchIdx];

        // Set up quadrature to degree+1 Gauss points per direction,
        // all lying on iter->side() except from the direction which
        // is NOT along the element

        gsGaussRule<T> bdQuRule(basis, 1.0, 1, iter->side().direction());

        // Create the iterator along the given part boundary.
        typename gsBasis<T>::domainIter bdryIter = basis.makeDomainIterator(iter->side());

        for(; bdryIter->good(); bdryIter->next() )
        {
            bdQuRule.mapTo( bdryIter->lowerCorner(), bdryIter->upperCorner(),
                            md.points, quWeights);

            //geoEval->evaluateAt( md.points );
            patch.computeMap(md);

            // the values of the boundary condition are stored
            // to rhsVals. Here, "rhs" refers to the right-hand-side
            // of the L2-projection, not of the PDE.
            rhsVals = iter_dir->function()->eval( m_pde_ptr->domain()[patchIdx].eval( md.points ) );
            rhsVals2 = iter->function()->eval( m_pde_ptr->domain()[patchIdx].eval( md.points ) );

            basis.eval_into( md.points, basisVals);
            basis.deriv_into( md.points, basisGrads);

            // Indices involved here:
            // --- Local index:
            // Index of the basis function/DOF on the patch.
            // Does not take into account any boundary or interface conditions.
            // --- Global Index:
            // Each DOF has a unique global index that runs over all patches.
            // This global index includes a re-ordering such that all eliminated
            // DOFs come at the end.
            // The global index also takes care of glued interface, i.e., corresponding
            // DOFs on different patches will have the same global index, if they are
            // glued together.
            // --- Boundary Index (actually, it's a "Dirichlet Boundary Index"):
            // The eliminated DOFs, which come last in the global indexing,
            // have their own numbering starting from zero.

            // Get the global indices (second line) of the local
            // active basis (first line) functions/DOFs:
            basis.active_into(md.points.col(0), globIdxAct );
            gsMatrix<unsigned > localIdx = globIdxAct;
            mapper.localToGlobal( globIdxAct, patchIdx, globIdxAct);

            // Out of the active functions/DOFs on this element, collect all those
            // which correspond to a boundary DOF.
            // This is checked by calling mapper.is_boundary_index( global Index )

            // eltBdryFcts stores the row in basisVals/globIdxAct, i.e.,
            // something like a "element-wise index"
            std::vector<index_t> eltBdryFcts;
            eltBdryFcts.reserve(mapper.boundarySize());
            for( index_t i=0; i < globIdxAct.rows(); i++)
                if( mapper.is_boundary_index(globIdxAct(i,0)) )
                    eltBdryFcts.push_back( i );

            //          Surface case
            if(md.dim.first + 1 == md.dim.second)
            {
                gsMatrix<T> geoMapDeriv1 = patch.deriv(md.points); // First derivative of the geometric mapping with respect to the parameter coordinates

                for( index_t k=0; k < md.points.cols(); k++ )
                {
                    // Compute the outer normal vector on the side
                    outerNormal(md, k, iter->side(), unormal);

                    gsMatrix<T> Jk = md.jacobian(k);
                    gsMatrix<T> G = Jk.transpose() * Jk;
                    gsMatrix<T> G_inv = G.cramerInverse();


                    real_t detG = G.determinant();

//                  Multiply quadrature weight by the measure of normal
                    const T weight_k = sqrt(detG) * quWeights[k] ;

                    unormal.normalize();

                    gsInfo << "rhsVals3: " << (Jk * G_inv).transpose() * unormal << "\n";

                    gsInfo << "unormal: " <<  unormal.transpose() << "\n";
                    gsInfo << "rhsVals2: " <<  rhsVals2.col(k).transpose() << "\n";



                    // Only run through the active boundary functions on the element:
                    for( size_t i0=0; i0 < eltBdryFcts.size(); i0++ )
                    {
                        // Each active boundary function/DOF in eltBdryFcts has...
                        // ...the above-mentioned "element-wise index"
                        const unsigned i = eltBdryFcts[i0];
                        // ...the boundary index.
                        const unsigned ii = mapper.global_to_bindex( globIdxAct( i ));

                        for( size_t j0=0; j0 < eltBdryFcts.size(); j0++ )
                        {
                            const unsigned j = eltBdryFcts[j0];
                            const unsigned jj = mapper.global_to_bindex( globIdxAct( j ));

                            // Use the "element-wise index" to get the needed
                            // function value.
                            // Use the boundary index to put the value in the proper
                            // place in the global projection matrix.
                            projMatEntries.add(ii, jj, weight_k * (basisVals(i,k) * basisVals(j,k)  +
                                lambda * ( ( (Jk * G_inv * basisGrads.block(2*i, k, 2, 1)).transpose() * unormal)(0,0) *
                                    ( (Jk * G_inv * basisGrads.block(2*j, k, 2, 1)).transpose() * unormal )(0,0) ) ) );
                        } // for j
                        globProjRhs.row(ii) += weight_k * ( basisVals(i,k) * rhsVals.col(k).transpose() );
                        globProjRhs.row(ii) += weight_k * ( lambda * ( (Jk * G_inv * basisGrads.block(2*i,k,2,1)).transpose() * unormal ) *
                            ( rhsVals2.col(k).transpose() * unormal) ); // unormal is different from planar, everthing else is the same! Does that makes sense?

                    } // for i
                } // for k
            }
            else
            {
                // Do the actual assembly:
                for (index_t k = 0; k < md.points.cols(); k++)
                {
                    // Compute the outer normal vector on the side
                    outerNormal(md, k, iter->side(), unormal);

                    // Multiply quadrature weight by the measure of normal
                    const T weight_k = quWeights[k] * md.measure(k);
                    unormal.normalize();

                    transformGradients(md, k, basisGrads, physBasisGrad);

                    // Only run through the active boundary functions on the element:
                    for (size_t i0 = 0; i0 < eltBdryFcts.size(); i0++)
                    {
                        // Each active boundary function/DOF in eltBdryFcts has...
                        // ...the above-mentioned "element-wise index"
                        const unsigned i = eltBdryFcts[i0];
                        // ...the boundary index.
                        const unsigned ii = mapper.global_to_bindex(globIdxAct(i));

                        for (size_t j0 = 0; j0 < eltBdryFcts.size(); j0++)
                        {
                            const unsigned j = eltBdryFcts[j0];
                            const unsigned jj = mapper.global_to_bindex(globIdxAct(j));

                            // Use the "element-wise index" to get the needed
                            // function value.
                            // Use the boundary index to put the value in the proper
                            // place in the global projection matrix.
                            projMatEntries.add(ii, jj, weight_k * (basisVals(i, k) * basisVals(j, k) + lambda *
                                ((physBasisGrad.col(i).transpose() * unormal)(0, 0)
                                    * (physBasisGrad.col(j).transpose() * unormal)(0, 0))));
                        } // for j

                        globProjRhs.row(ii) += weight_k * (basisVals(i, k) * rhsVals.col(k).transpose() + lambda *
                            (physBasisGrad.col(i).transpose() * unormal) * (rhsVals2.col(k).transpose() * unormal));

                    } // for i
                } // for k
            }
        } // bdryIter
        iter_dir++;
    } // boundaryConditions-Iterator

    gsSparseMatrix<T> globProjMat( mapper.boundarySize(), mapper.boundarySize() );
    globProjMat.setFrom( projMatEntries );
    globProjMat.makeCompressed();

    // Solve the linear system:
    // The position in the solution vector already corresponds to the
    // numbering by the boundary index. Hence, we can simply take them
    // for the values of the eliminated Dirichlet DOFs.
    typename gsSparseSolver<T>::CGDiagonal solver;

    m_ddof[0] = solver.compute( globProjMat ).solve ( globProjRhs );
    //m_ddof[0].setZero();

    gsInfo << "Boundary: " << m_ddof[0] << "\n";
}
// End


} // namespace gismo



