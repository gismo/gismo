
#pragma once

#include <ostream>
#include <gsCore/gsDebug.h>

#include <gsCore/gsDofMapper.h>
#include <gsAssembler/gsGaussAssembler.h>
#include <gsUtils/gsInterpolate.h>

namespace gismo
{

template<class T>
gsGalerkinMethod<T>::~gsGalerkinMethod()
{
  delete m_assembler;
  delete m_dofMapper;
  m_system.free();

  freeAll(m_bases);
}

template<class T>
void gsGalerkinMethod<T>::assemble()
{
    // If the Dirichlet strategy is elimination then precompute
    // Dirichlet dofs (m_dofMapper excludes these from the system
    if ( m_dirStrategy == dirichlet::interpolation )
        computeDirichletDofs();
    
    if (m_idofs == 0 ) // Are there any interior dofs ? 
    {
        gsWarn << " No internal DOFs. Computed dirichlet boundary only.\n" <<"\n" ;
        // Create empty system, for compatibility
        m_system = gsSparseSystem<T>( new gsSparseMatrix<T>, new gsVector<T> );
        return;
    }

    // Assemble the system matrix and right-hand side
    if (m_bvp->nPatches() == 1)
    {
        m_assembler->setGeometry( m_bvp->patch(0) );
        m_system = 
            m_assembler->assemble( *m_bases[0], dofMapper(), m_ddof, m_bvp->pde() );
    }
    else
        assembleMultipatch();
    
    // Enforce Dirichlet boundary conditions by Nitsche's method
    if ( m_dirStrategy == dirichlet::nitsche )
    {
        for ( typename gsBVProblem<T>::const_iterator 
                  it = m_bvp->dirichletBegin(); it != m_bvp->dirichletEnd(); ++it )
        {
            m_assembler->setGeometry( m_bvp->patch(it->patch()) );
            m_assembler->applyBoundary( basis(it->patch()), *it, dofMapper(), m_system );
        }
    }
    
    // Enforce Neumann boundary conditions
    for (typename gsBVProblem<T>::const_iterator it= m_bvp->neumannBegin();
         it != m_bvp->neumannEnd(); ++it)
    {
        m_assembler->setGeometry( m_bvp->patch(it->patch()) );
        m_assembler->applyBoundary( basis(it->patch()), *it,
                                    dofMapper(), m_system );
    }

    // If we are in in dg (Discontinous Galerkin) mode add interface
    // contributions
    if ( m_interfaceStrategy == iFace::dg )
    {
        for ( typename gsMultiPatch<T>::const_iiterator it = 
                    this->m_bvp->patches().iBegin();
              it != this->m_bvp->patches().iEnd(); ++it )
        {
            m_assembler->setGeometry( m_bvp->patch( it->ps1.patch ) );
            m_assembler->applyDG( basis(it->ps1.patch), basis(it->ps2.patch),
                                  m_bvp->patch( it->ps2.patch ),
                                  *it, dofMapper(), m_system );
        }
    }
}


template<class T>
void gsGalerkinMethod<T>::assembleMultipatch()
{
    gsSparseMatrix<T> * K = new gsSparseMatrix<T>(m_idofs,m_idofs) ;
    gsVector<T>       * b = new gsVector<T>(m_idofs);
    b->setZero();

    for (int np=0; np < m_bvp->nPatches(); ++np )
    {
        m_assembler->setGeometry( m_bvp->patch(np) );

        // assemble stiffness matrix and rhs for the local patch np
        gsSparseSystem<T> patchSystem =
          m_assembler->assemble( *m_bases[np], dofMapper(), m_ddof, m_bvp->pde(), np );

        // add result to the global system (K,b)
        *K += *patchSystem.matrix();
        *b += *patchSystem.rhs();

        patchSystem.free();
    }

    K->makeCompressed();

    this->m_system = gsSparseSystem<T>(K, b);
}


template<class T>
gsSparseMatrix<T>* gsGalerkinMethod<T>::assembleMass() const
{
    if (m_idofs == 0 )
    {
        gsWarn << " No internal DOFs.\n";
        return 0;
    }

    if (m_bvp->nPatches() == 1)
      {
        m_assembler->setGeometry( m_bvp->patch(0) );
        return m_assembler->assembleMass( *m_bases[0], dofMapper() );
      }
    else
        GISMO_ERROR( "not implemented" );
        //assembleMultipatch();
}


template<class T>
gsField<T> * gsGalerkinMethod<T>::solve()
{
    // Assemble system
    if (!m_system.matrix())
        assemble();

    // Solve linear system
    gsVector<T> res;

    if ( m_idofs > 0 )
    {
        if (m_bvp->pde().isSymmetric())
        {
            Eigen::ConjugateGradient< gsSparseMatrix<T> > solver;
            //Eigen::ConjugateGradient< gsSparseMatrix<T>, Eigen::Lower, Eigen::IncompleteLUT<T> > solver;
            
            solver.setMaxIterations(4 * m_system.matrix()->rows() );

            res = solver.compute( *m_system.matrix() ).solve ( * m_system.rhs() );
            gsInfo << "residual error: " << solver.error() << "\n";
	    ls_iterations = solver.iterations();
            gsInfo << "    iterations: " << solver.iterations() << "\n";
	    gsInfo << "    Tolerance: " << solver.tolerance() << "\n"; 
        }
        else
        {
            Eigen::BiCGSTAB< gsSparseMatrix<T>, Eigen::IncompleteLUT<T> > solver;
            //Eigen::BiCGSTAB< gsSparseMatrix<T>, Eigen::IdentityPreconditioner > solver;
            //Eigen::BiCGSTAB< gsSparseMatrix<T>, Eigen::DiagonalPreconditioner > solver;

            //Eigen::SparseQR<gsSparseMatrix<T>, Eigen::COLAMDOrdering<index_t> >  solver;
            //solver.analyzePattern(*m_system.matrix());
            //solver.factorize( *m_system.matrix() );
            //res = solver.solve ( * m_system.rhs() );

            res = solver.compute( * m_system.matrix() ).solve ( * m_system.rhs() );
            gsInfo << "residual error: " << solver.error() << "\n";
	    ls_iterations = solver.iterations();
            gsInfo << "    iterations: " << solver.iterations() << "\n";
	    gsInfo << "    Tolerance: " << solver.tolerance() << "\n"; 
        }
        //Eigen::SimplicialLDLT< gsSparseMatrix<T> > solver;
        //Eigen::SparseLU<gsSparseMatrix<T>, Eigen::COLAMDOrdering<index_t> >  solver;
    }
    
    return reconstructSolution(res);
}



template<class T>
gsField<T> * gsGalerkinMethod<T>::reconstructSolution(const gsVector<T>& data) const
{
    if (m_bvp->nPatches() == 1)     // single-patch case
    {
        return new gsField<T>( this->m_bvp->patches(),
                               reconstructPatchSolution(data, 0) );
    }
    else        // multipatch case
    {
        std::vector<gsFunction<T> * > sols ;

        const gsMultiPatch<T> & mp = this->m_bvp->patches();

        for (size_t np=0; np < mp.nPatches(); ++np )
        {    
            sols.push_back( reconstructPatchSolution( data, np ) );
        }

        return new gsField<T>( mp, sols );
    }
}


template<class T>
gsGeometry<T> * gsGalerkinMethod<T>::reconstructPatchSolution(const gsVector<T>& data, int p) const
{
    // Reconstruct solution coefficients on patch p
    const int sz  = m_bases[p]->size();
    const int fsz = dofMapper().freeSize();
    const int dim = m_bvp->pde().fieldDim();

    gsMatrix<T> coeffs(sz, dim);
    
    for (index_t i = 0; i < sz; ++i)
    {
        if ( dofMapper().is_free(i, p) ) // internal or interface
        {
            for (int k = 0; k < dim; ++k)
                coeffs(i,k) = data[ k * fsz + dofMapper().index(i, p) ];
        }
        else // eliminated Dof: fill with Dirichlet data
        {
            coeffs.row(i) = m_ddof.row( dofMapper().bindex(i, p) );
        }
    }
    return m_bases[p]->makeGeometry( give(coeffs) );
}


template<class T>
void gsGalerkinMethod<T>::initAssembler()
{
    if (m_assembler)
        delete m_assembler;

    m_assembler = new gsGaussAssembler<T>( m_bvp->patch(0) );
}


template<class T>
void gsGalerkinMethod<T>::initDofMapper()
{
    // Setup the degrees of freedom, interface matching etc
    if (m_dofMapper)
        delete m_dofMapper;
    
    if ( m_dirStrategy == dirichlet::elimination || 
         m_dirStrategy == dirichlet::homogeneous)
        m_dofMapper = new gsDofMapper( gsMultiBasis<T>(m_bases, this->m_bvp->patches()), 
                                       m_bvp->boundaryConditions() );
    else
        m_dofMapper = new gsDofMapper( gsMultiBasis<T>(m_bases, this->m_bvp->patches()) );
    if ( m_interfaceStrategy == iFace::glue )
    {
        for ( gsBoxTopology::const_iiterator it = this->m_bvp->patches().iBegin();
              it != this->m_bvp->patches().iEnd(); ++it )
        {
            gsMatrix<unsigned>
                * b1= m_bases[it->ps1.patch]->boundary( it->ps1.side ),
                * b2= m_bases[it->ps2.patch]->boundary( it->ps2.side );
            
            m_dofMapper->matchInterface( it->ps1.patch, it->ps2.patch, *b1, *b2, it->orient);
            
            delete b1;
            delete b2;
        }
    }
    m_dofMapper->finalize();

    m_idofs = m_dofMapper->freeSize();
    m_dofs  = m_dofMapper->size();
}



template<class T>
void gsGalerkinMethod<T>::computeDirichletDofs()
{
    m_ddof.resize( dofMapper().boundarySize(), m_bvp->pde().fieldDim() );

    for ( typename gsBVProblem<T>::const_iterator 
              it = m_bvp->dirichletBegin(); it != m_bvp->dirichletEnd(); ++it )
    {
        // Get dofs on this boundary
        gsMatrix<unsigned> * boundary = m_bases[it->patch()]->boundary(it->side()) ;

        // If the condition is homogeneous then fill with zeros
        if ( it->isHomogeneous() )
        {
            for (index_t k=0; k!= boundary->size(); ++k)
            {
                const int ii= dofMapper().bindex( (*boundary)(k) , it->patch() );
                m_ddof.row(ii).setZero();
            }
            delete boundary;
            continue;
        }

        // Get the side information
        int dir = direction( it->side() );
        index_t param = (parameter( it->side() ) ? 1 : 0);

        // Compute grid of points on the face ("face anchors")
        std::vector< gsVector<T> > rr;
        rr.reserve( this->patches().parDim() );

        for ( int i=0; i < this->patches().parDim(); ++i)
        {
            if ( i==dir )
            {
                gsVector<T> b(1); 
                b[0] = ( m_bases[it->patch()]->component(i).support() ) (0, param);
                rr.push_back(b);
            }
            else
            {   
                rr.push_back( m_bases[it->patch()]->component(i).anchors()->transpose() );
            }
        }

        // Compute dirichlet values
        gsMatrix<T> fpts = 
            it->function()->eval( m_bvp->patch(it->patch()).eval(  gsPointGrid( rr ) ) );

        // Interpolate dirichlet boundary 
        gsBasis<T> * h = m_bases[it->patch()]->boundaryBasis(it->side());
        gsGeometry<T> * geo = gsInterpolate( *h, *h->anchors(), fpts );
        const gsMatrix<T> & dVals =  geo->coefs();

        // Save corresponding boundary dofs
        for (index_t k=0; k!= boundary->size(); ++k)
        {
            const int ii= dofMapper().bindex( (*boundary)(k) , it->patch() );
            m_ddof.row(ii) = dVals.row(k);
        }
        delete h;
        delete geo;
        delete boundary;
    }
}


} // namespace gismo
