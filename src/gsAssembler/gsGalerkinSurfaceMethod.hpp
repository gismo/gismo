// Assembler for Surface problems R^2 ---> R^3
// Create private classes 
// Galerkin Surface Method to include 
// Surface Convection-Diffusion-Reaction methods and Parabolic methods.

#pragma once

#include <ostream>
#include <gsCore/gsDebug.h>

#include <gsCore/gsDofMapper.h>
#include <gsUtils/gsInterpolate.h>
#include <gsCore/gsMemory.h>


namespace gismo
{

template<class T>
gsGalerkinSurfaceMethod<T>::~gsGalerkinSurfaceMethod()
{
  // TODO: should we delete the m_bases?
  delete m_assembler;
  delete m_dofMapper;
  m_system.free();
}
//-------------------------------------------------------------//

template<class T>
void gsGalerkinSurfaceMethod<T>::assemble()
{
    if ( m_dirStrategy == elimination )
        computeDirichletDofs();
    
    if (m_idofs == 0 )
    {
        gsWarn << " No internal DOFs. Computed dirichlet boundary only.\n" <<"\n" ;
        return;
    }

    if (m_bvp->nPatches() == 1)
      {
        m_assembler->setGeometry( m_bvp->patch(0) );
        m_system = m_assembler->assemble( *m_bases[0], dofMapper(), m_ddof, m_bvp->pde());
      }
    else
        assembleMultipatch();

    // Enforce Dirichlet boundary conditions by Nitsche's method
    if ( m_dirStrategy == nitsche )
    {
        for ( typename gsBVProblem<T>::const_iterator 
                  it = m_bvp->dirichletBegin(); it != m_bvp->dirichletEnd(); ++it )
        {
            m_assembler->setGeometry( m_bvp->patch(it->patch()) );
            m_assembler->applyBoundary( basis(it->patch()), *it,
                                        dofMapper(), m_system );
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
}
//--------------------------------------------------------//

template<class T>
void gsGalerkinSurfaceMethod<T>::assembleMultipatch()
{
    gsSparseMatrix<T> * K = new gsSparseMatrix<T>(m_idofs,m_idofs) ;
    gsVector<T>       * b = new gsVector<T>(m_idofs);
    b->setZero();

    for (int np=0; np < m_bvp->nPatches(); ++np )
    {
        // TO DO: add to the Gauss assembler the DIRICHLET DOFS, so
        // that they are excluded from the DoFs of the linear system,
        // eg. to use .stiffness() method

        m_assembler->setGeometry( m_bvp->patch(np) );

        // assemble stiffness matrix and rhs for the local patch np
        gsSparseSystem<T> patchSystem = m_assembler->assemble(
          basis(np), dofMapper(), m_ddof, m_bvp->pde(), np);

        // add result to the global system (K,b)
        *K += *patchSystem.matrix();
        *b += *patchSystem.rhs();

        patchSystem.free();
    }

    K->makeCompressed();

    this->m_system = gsSparseSystem<T>(K, b);
}

//---------------------------------------------------------------//

template<class T>
gsField<T> * gsGalerkinSurfaceMethod<T>::solve()
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
            res = solver.compute( * m_system.matrix() ).solve ( * m_system.rhs() );
//             gsInfo << "residual error: " << solver.error() << "\n";
            gsInfo << "    iterations: " << solver.iterations() << "\n";
	}
        else
        {
            Eigen::BiCGSTAB< gsSparseMatrix<T>, Eigen::IncompleteLUT<T> > solver;
            res = solver.compute( * m_system.matrix() ).solve ( * m_system.rhs() );
            gsInfo << "residual error: " << solver.error() << "\n";
            gsInfo << "    iterations: " << solver.iterations() << "\n";
        }
        //Eigen::SimplicialLDLT< gsSparseMatrix<T> > solver;
        //Eigen::SparseLU<gsSparseMatrix<T>, Eigen::COLAMDOrdering<index_t> >  solver;
    }
    
    return reconstructSolution(res);
}
    
//------------------------------------------------------------------//

template<class T>
gsField<T> * gsGalerkinSurfaceMethod<T>::reconstructSolution(const gsVector<T>& data) const
{
    if (m_bvp->nPatches() == 1)     // single-patch case
    {
        gsBasis<T> * basis = m_bases[0];

        // Reconstruct solution coefficients
        gsMatrix<T> sol(basis->size(), 1);

        for ( index_t i = 0; i < basis->size(); ++i )
        {
            sol(i) = dofMapper().is_free(i)
                ? data[ dofMapper().index(i) ]
                : m_ddof( dofMapper().bindex(i) );
        }
        
        // energy norm should be of O(h) -- Energy norm!
        //std::cout<<"\nError norm: \n" <<  sol->transpose() * (*K) * (*sol) << std::endl;

        return new gsField<T>( this->m_bvp->patches(), basis->makeGeometry( give(sol) ) );
    }
    else        // multipatch case
    {
        // Reconstruct solution coefficients
        std::vector<gsFunction<T> * > sols ;

        const gsMultiPatch<T> & mp = this->m_bvp->patches();

        for (size_t np=0; np < mp.nPatches(); ++np )
        {    
            gsMatrix<T> sol(m_bases[np]->size(), 1);
            for ( index_t i= 0; i< m_bases[np]->size(); ++i )
            {
                sol(i) = dofMapper().is_free(i, np)    // internal or interface
                    ? data[ dofMapper().index(i,np) ]
                    : m_ddof( dofMapper().bindex(i,np) );
            }
            sols.push_back( m_bases[np]->makeGeometry( give(sol) ) );
        }

        return new gsField<T>( mp, sols );
    }
}

//--------------------------------------------------------------------------//

template<class T>
void gsGalerkinSurfaceMethod<T>::initAssembler()
{
    if (m_assembler)
        delete m_assembler;
    
    gsGaussSurfaceAssembler<T> * ga = new gsGaussSurfaceAssembler<T>( & m_bvp->patch(0) );

    gsVector<int> int_nodes(d);
    for (int i = 0; i != d  ; ++i) 
        int_nodes(i) = m_bases[0]->degree(i) + 1;
    ga->setNumIntNodes( int_nodes );

    m_assembler = ga;
}

//---------------------------------------------------------------------------//

template<class T>
void gsGalerkinSurfaceMethod<T>::initDofMapper()
{
    // Setup the degrees of freedom, interface matching etc
    if (m_dofMapper)
       delete m_dofMapper;
    
    m_dofMapper = new gsDofMapper( m_bases );
    m_dofMapper->setMatchingInterfaces( this->m_bvp->patches() );

    // Mark boundary degrees of freedom in case of elimination of Dirichlet dofs
    if ( m_dirStrategy == elimination )
        for ( typename gsBVProblem<T>::const_iterator 
                  it = m_bvp->dirichletBegin(); it != m_bvp->dirichletEnd(); ++it )
        {
            m_dofMapper->markBoundary( it->ps );   
        }
        
    m_dofMapper->finalize();
    m_idofs = m_dofMapper->freeSize();
    m_dofs  = m_dofMapper->size();
}

//------------------------------------------------------------------------------//

template<class T>
void gsGalerkinSurfaceMethod<T>::computeDirichletDofs()
{
    m_ddof.resize(dofMapper().boundarySize(), 1);

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
		m_ddof(ii, 0) = T(0.0);
	      }
	    delete boundary;	    
	    continue;
	  }

        // Get the side information
        int dir = direction( it->side() );
        index_t param = (parameter( it->side() ) ? 1 : 0);

        // Compute grid of points on the face ("face anchors")
        std::vector< gsVector<T> > rr;
        rr.reserve(this->geometry().parDim());

        for ( int i=0; i!=this->geometry().parDim(); ++i)
            if ( i==dir )
            {
                gsVector<T>  b(1); 
                b[0] = ( m_bases[it->patch()]->component(i).support() ) (0, param);
                rr.push_back(b);
            }
            else
            {   
                rr.push_back(m_bases[it->patch()]->component(i).anchors()->transpose() );
            }

        // Compute dirichlet values
        gsMatrix<T> fpts = m_bvp->patch(it->patch()).eval( gsPointGrid( rr ));

        // Interpolate dirichlet boundary 
        gsBasis<T> * h = m_bases[it->patch()]->boundaryBasis(it->side());
        gsGeometry<T> * geo = gsInterpolate( *h, *h->anchors(), fpts );

        // Save corresponding boundary dofs
        for (index_t k=0; k!= boundary->size(); ++k)
        {
            const int ii= dofMapper().bindex( (*boundary)(k) , it->patch() );
            m_ddof(ii, 0) = geo->coefs()(k,0);
        }
        delete h;
        delete geo;
        delete boundary;
    }
}

//----------------------------------------------------------------//

template<class T>
void gsGalerkinSurfaceMethod<T>::addNitscheBoundary()
{
    for ( typename gsBVProblem<T>::const_iterator 
              it = m_bvp->dirichletBegin(); it != m_bvp->dirichletEnd(); ++it )
    {
        m_assembler->setGeometry( m_bvp->patch(it->patch()) );
        m_assembler->boundaryNeumann( basis(it->patch()), *it,
                                      dofMapper(), m_system );
    }
}

//-----------------------------------------------------------------//

template<class T>
void gsGalerkinSurfaceMethod<T>::addNeumannConditions()
{
    // iterate over all Neumann boundaries
    for (typename gsBVProblem<T>::const_iterator it= m_bvp->neumannBegin();
         it != m_bvp->neumannEnd(); ++it)
    {
        m_assembler->setGeometry( m_bvp->patch(it->patch()) );
        m_assembler->boundaryNeumann( basis(it->patch()), *it,
                                      dofMapper(), m_system );
    }
}
} // namespace gismo
