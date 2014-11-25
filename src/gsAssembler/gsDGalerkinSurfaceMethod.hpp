// Assembler for Discontinuous Galerkin Surface problems R^2 ---> R^3
// Create private classes 
// Galerkin Surface Method to include 
// Surface Convection-Diffusion-Reaction problem as well as Parabolic problem.

#pragma once

#include <ostream>
#include <gsCore/gsDebug.h>

#include <gsCore/gsDofMapper.h>
#include <gsUtils/gsInterpolate.h>
#include <gsCore/gsMemory.h>


namespace gismo
{

template<class T>
gsDGalerkinSurfaceMethod<T>::~gsDGalerkinSurfaceMethod()
{
  // TODO: should we delete the m_bases?
  delete m_assembler;
  delete m_dofMapper;
  m_system.free();
}
//-------------------------------------------------------------//

template<class T>
void gsDGalerkinSurfaceMethod<T>::assemble()
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
        m_assembler->setGeometry( &m_bvp->patch(0) );
        m_system = m_assembler->assembleSurfacePoisson( *m_bases[0], dofMapper(), m_ddof, static_cast<gsPoissonPde<T> &>( m_bvp->pde() ));
      }
    else
        assembleMultipatch();

    if ( m_dirStrategy == nitsche )
        addNitscheBoundary();

}
//--------------------------------------------------------//

template<class T>
void gsDGalerkinSurfaceMethod<T>::assembleMultipatch()
{
    gsSparseMatrix<T> * K = new gsSparseMatrix<T>(m_idofs,m_idofs) ;
    gsVector<T>       * b = new gsVector<T>(m_idofs);
    b->setZero();

    for (int np=0; np < m_bvp->nPatches(); ++np )
    {
        // TO DO: add to the Gauss assembler the DIRICHLET DOFS, so
        // that they are excluded from the DoFs of the linear system,
        // eg. to use .stiffness() method

        m_assembler->setGeometry( &m_bvp->patch(np) );

        // assemble stiffness matrix and rhs for the local patch np
        gsSparseSystem<T> patchSystem = m_assembler->assembleSurfacePoisson(
        *m_bases[np], dofMapper(), m_ddof, static_cast<gsPoissonPde<T>&>(m_bvp->pde() ), np);

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
gsField<T> * gsDGalerkinSurfaceMethod<T>::solve()
{
    // Assemble system
    if (!m_system.matrix())
      assemble();

    // Enforce Neumann boundary conditions
    addNeumannConditions();
    
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
gsField<T> * gsDGalerkinSurfaceMethod<T>::reconstructSolution(const gsVector<T>& data) const
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

        return new gsField<T>( &geometry(), basis->makeGeometry( &sol ) );
    }
    else        // multipatch case
    {
        // Reconstruct solution coefficients
        std::vector<gsFunction<T> * > sols ;
        gsMatrix<T> sol;

        const gsMultiPatch<T> & mp = this->m_bvp->patches();

        for (size_t np=0; np < mp.nPatches(); ++np )
        {    
            sol.resize(m_bases[np]->size(), 1);
            for ( index_t i= 0; i< m_bases[np]->size(); ++i )
            {
                sol(i) = dofMapper().is_free(i, np)    // internal or interface
                    ? data[ dofMapper().index(i,np) ]
                    : m_ddof( dofMapper().bindex(i,np) );
            }
            sols.push_back( m_bases[np]->makeGeometry( &sol ) );
        }

        return new gsField<T>( mp, sols );
    }
}

//------------------------------------------------------------------//

template<class T>
void gsDGalerkinSurfaceMethod<T>::initAssembler()
{
    gsGaussSurfaceAssembler<T> * ga = new gsGaussSurfaceAssembler<T>( & m_bvp->patch(0) );

    gsVector<int> int_nodes(d);
    for (int i = 0; i != d  ; ++i) 
        int_nodes(i) = m_bases[0]->degree(i) + 1;
    ga->setNumIntNodes( int_nodes );

    m_assembler = ga;
}

//----------------------------------------------------------------------//

template<class T>
void gsDGalerkinSurfaceMethod<T>::initDofMapper()
{
    if (m_dofMapper)
        delete m_dofMapper;

    gsBoundaryConditions<T>  empty_bc;
    gsBoxTopology            empty_topology;

    // initialize options to do nothing
    const gsBoundaryConditions<T> *bc      =&empty_bc;
    const gsBoxTopology           *topology=&empty_topology;

     // set interfaces to glue to the patch interfaces
    if ( m_interfaceStrategy == glue )
        topology = &(this->m_bvp->patches());

    // set boundary dofs to eliminate corresponding to Dirichlet sides
    if ( m_dirStrategy == elimination )
        bc = &(m_bvp->boundaryConditions());

    m_dofMapper = new gsDofMapper( m_bases, *topology, *bc );

    m_idofs = m_dofMapper->freeSize();
    m_dofs  = m_dofMapper->size();
}

//-------------------------------------------------------------------//

template<class T>
void gsDGalerkinSurfaceMethod<T>::computeDirichletDofs()
{
    m_ddof.resize(dofMapper().boundarySize(), 1);
    //std::cout<< "# boundary dofs: "<< ddof->rows() <<"\n";

    for ( typename gsBVProblem<T>::const_iterator 
              it = m_bvp->dirichletBegin(); it != m_bvp->dirichletEnd(); ++it )
    {
        // Get the side information
        int dir = direction( it->side() );
        index_t param = (parameter( it->side() ) ? 1 : 0);

        // Compute grid of points on the face ("face anchors")
        std::vector< gsVector<T>* > rr;
        rr.reserve(this->geometry().parDim());

        for ( int i=0; i!=this->geometry().parDim(); ++i)
            if ( i==dir )
            {
                gsVector<T> * b = new gsVector<T>(1); 
                (*b)[0] = ( m_bases[it->patch()]->component(i).support() ) (0, param);
                rr.push_back(b);
            }
            else
            {   
                //gsVector<T> * b = new gsVector<T>(h->anchors()->transpose()); 
                gsVector<T> * b = new gsVector<T>( 
                    safe( m_bases[it->patch()]->component(i).anchors() )->transpose() );
                rr.push_back(b);
            }

        // Compute dirichlet values
        gsMatrix<T> fpts   =  it->function()->eval(
            m_bvp->patch(it->patch()).eval( *gsPointGrid( rr )));

        // free rr
        for (unsigned i = 0; i < rr.size(); ++i)
          delete rr[i];

        // Interpolate dirichlet boundary 
        gsBasis<T> * h = m_bases[it->patch()]->boundaryBasis(it->side());
        gsGeometry<T> * geo = gsInterpolate( *h, *safe(h->anchors()), fpts );
        // Save corresponding boundary dofs
        gsMatrix<unsigned> * boundary = m_bases[it->patch()]->boundary(it->side()) ;
        for (index_t k=0; k!= boundary->size(); ++k)
        {
            int ii= dofMapper().bindex( (*boundary)(k) , it->patch() );
            m_ddof(ii, 0) = geo->coefs()(k,0);
        }
        delete h;
        delete geo;
        delete boundary;
    }
}

//----------------------------------------------------------------//

template<class T>
void gsDGalerkinSurfaceMethod<T>::addNitscheBoundary()
{
    gsSparseMatrix<T> * lhs = m_system.matrix();    
    //gsVector<T> *rhs = m_system.rhs(); // TO DO: Implement Nitsche for non-zero Dir. boundary as well

    for ( typename gsBVProblem<T>::const_iterator 
              it = m_bvp->dirichletBegin(); it != m_bvp->dirichletEnd(); ++it )
    {
        // int dir = direction( it->side() );
        // index_t param = (parameter( it->side() ) ? 1 : 0); // 
        gsBasis<T>  * bbasis = m_bases[it->patch()]->boundaryBasis(it->side() );
        m_assembler->setGeometry( &m_bvp->patch(it->patch()) );
	m_assembler->boundaryNitsche( bbasis,* it->function(), it->side() , lhs );
    

    }
}

//-----------------------------------------------------------------//

template<class T>
void gsDGalerkinSurfaceMethod<T>::addNeumannConditions()
{
    gsVector<T> *rhs = m_system.rhs();
    //gsLog<< "rhs = " << rhs->transpose() << "\n";
    //m_dofMapper->print();

    // Add Neumann boundary conditions
    for (typename gsBVProblem<T>::const_iterator it= m_bvp->neumannBegin();
           it != m_bvp->neumannEnd(); ++it)
    {
        gsBasis<T>  * bb = m_bases[it->patch()]->boundaryBasis(it->side() );
        m_assembler->setGeometry( &m_bvp->patch(it->patch()) );
        gsVector<T> * nmnn = m_assembler->boundaryMoments( bb , * it->function(), it->side() );

        gsWarn << "Neumann values: \n" << nmnn->transpose() << "\n";
        
        gsMatrix<unsigned> *ind = this->basis().boundary(it->side()) ;
        
        for (index_t i = 0; i < ind->size(); ++i)
        {
            unsigned ii = m_dofMapper->index((*ind)(i), it->patch() ); // convert local dof index to global dof index
            gsWarn<<"ii= " << ii<< ",i=" << i << ",ind = " << (*ind)(i,0)<< "\n" ;
            if (m_dofMapper->is_free_index(ii)) // exclude Dirichlet dof from b
                (*rhs)[ii] -= (*nmnn)[i];
        }

        delete ind;
        delete bb;
        delete nmnn;
    }
}
} // namespace gismo
