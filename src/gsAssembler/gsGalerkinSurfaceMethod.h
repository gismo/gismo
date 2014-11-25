#pragma once

#include <vector>

#include <gsAssembler/gsGaussSurfaceAssembler.h>
#include <gsCore/gsField.h>
// #include <gsCore/gsDofMapper.h>
#include <gsPde/gsBVProblem.h>


namespace gismo
{

    enum gsDirichletStrategy
    {
        elimination = 1, nitsche = 2
    };

  /** 
      Galerkin Surface Poisson solver
  */

//-----------------------------------------------------------------------------------//  
template<class T>
class gsGalerkinSurfaceMethod
{
public:
    typedef gismo::gsDirichletStrategy DirichletStrategy;
public:

  /// Default empty constructor
  gsGalerkinSurfaceMethod() { };

    //  Single Patch or multipatch with common basis
    gsGalerkinSurfaceMethod( gsBVProblem<T> * bvp, gsBasis<T> * basis, 
		       DirichletStrategy dstrategy = elimination)
        : m_bvp(bvp)
        {
            const gsMultiPatch<T> & mp = bvp->patches();
            mp.checkConsistency();
            d = mp.patch(0).parDim();
            
            // Collect the given basis, to be used for all patches
	    m_bases.push_back( basis );

	    // Initialize assember and dof mapper
            m_dirStrategy = dstrategy;
            m_assembler = NULL;
            initAssembler();
            m_dofMapper = NULL;
            initDofMapper();
        }

    //  Single Patch or multipatch using the geometry basis
    gsGalerkinSurfaceMethod( gsBVProblem<T> * bvp, DirichletStrategy dstrategy = elimination ) 
        : m_bvp(bvp)
        {
            const gsMultiPatch<T> & mp = bvp->patches();
            mp.checkConsistency();
            d = mp.patch(0).parDim();
            
	    // Collect & clone the basis of the geometry
	    m_bases.resize( mp.nPatches() );
	    for (size_t i=0; i< mp.nPatches(); ++i ) 
		m_bases[i] = mp.patch(i).basis().clone();
            
	    // Initialize assember and dof mapper
            m_dirStrategy = dstrategy;
            m_assembler = NULL;
            initAssembler();
            m_dofMapper = NULL;
            initDofMapper();
        }

    // Multipatch problem
    gsGalerkinSurfaceMethod( gsBVProblem<T> * bvp, 
		       std::vector<gsBasis<T> *> bases, 
		       DirichletStrategy dstrategy = elimination ) 
        : m_bvp(bvp), m_bases(bases)
        {
            assert( (int)bases.size() == bvp->nPatches() );

            const gsMultiPatch<T> & mp = bvp->patches();
            mp.checkConsistency();
            d = mp.patch(0).parDim();

	    // Initialize assember and dof mapper
            m_dirStrategy = dstrategy;
            m_assembler = NULL;
            initAssembler();
            m_dofMapper = NULL;
            initDofMapper();
        }

    //Destructor
  ~gsGalerkinSurfaceMethod();
 
//-----------------------------------------------------------------------------------// 
public:
    
    void assemble();
    
//-----------------------------------------------------------------------------------//    
    gsField<T> * solve();


 //-----------------------------------------------------------------------------------//
    gsGeometry<T>   & geometry() const  { return m_bvp->patches().patch(0); }

    const gsBasis<T>& basis(int i = 0) const    { return *m_bases[i]; }
    const gsDofMapper& dofMapper() const     { return *m_dofMapper; }
    const gsMatrix<T>& boundaryValues() const   { return m_ddof; }

    /// Total number of degrees of freedom, including Dirichlet dofs
    int dofs() const                            { return m_dofs; }

    /// Number of degrees of freedom which we have to solve for
    int freeDofs() const                        { return m_idofs; }

    void set_assembler( gsAssembler<T>* assblr)  
    { 
         if( m_assembler ) 
	    delete m_assembler;
         m_assembler= assblr; 
    }

          gsAssembler<T>& assembler()          { return *m_assembler; }
    const gsAssembler<T>& assembler() const    { return *m_assembler; }

    const gsSparseSystem<T>& linearSystem() const       { return m_system; }

    /// Given a vector of the free dof coefficients, complete it with the Dirichlet data and return a solution field
  gsField<T> * reconstructSolution(const gsVector<T>& data) const;


//-----------------------------------------------------------------------------------//    
    // void setDirichletStrategy( const DirichletStrategy & strategy){ m_dirStrategy = strategy; }
    // DirichletStrategy getDirichletStrategy()                { return m_dirStrategy; }
    
private:

  void initAssembler();
    
  void initDofMapper();

//-----------------------------------------------------------------------------------//
  void computeDirichletDofs();
//-----------------------------------------------------------------------------------//    
  void addNitscheBoundary();
 //-----------------------------------------------------------------------------------//   
  void addNeumannConditions();
//-----------------------------------------------------------------------------------//
  void assembleMultipatch();

 //-----------------------------------------------------------------------------------//
// Problem & solver data
private:

    //BV Problem
    gsBVProblem<T> * m_bvp;
    int d;

    // Discretization data
    std::vector<gsBasis<T> *> m_bases ;

    // Assembly data
    gsGaussSurfaceAssembler<T> * m_assembler;
    gsDofMapper * m_dofMapper;

    // Strategy for dirichlet Dofs
    DirichletStrategy m_dirStrategy;
    
    // data for Dirichlet dofs by Elimination
    gsMatrix<T> m_ddof;

    // Linear solving  data
    // CG
    gsSparseSystem<T> m_system;
    
    // number of degrees of freedom
    int m_dofs;
    int m_idofs;

}; // class gsClass


//////////////////////////////////////////////////
//////////////////////////////////////////////////


} // namespace gismo

// #ifndef GISMO_HEADERS_ONLY

#include "gsGalerkinSurfaceMethod.hpp"

// #endif

