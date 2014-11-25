#pragma once

#include <vector>

#include <gsCore/gsField.h>
#include <gsAssembler/gsAssembler.h>
#include <gsAssembler/gsAssemblerOptions.h>
#include <gsPde/gsBVProblem.h>


namespace gismo
{

/** @brief
    Solves a boundary value problem using a given discretization basis.

    This class obtains a gsBVProblem as well as one or several discretization
    bases and solves the discretized PDE. It sets up an assembler and
    assembles the system patchwise and combines the patch-local stiffness
    matrices into a global system by various methods (see gismo::gsInterfaceStrategy).
    It can also enforce Dirichlet boundary conditions in various ways
    (see gismo::gsDirichletStrategy).

    By calling gsGalerkinMethod::solve(), the linear system is assembled
    and solved, and a gsField containing the computed solution is returned.
*/


/// \todo pass all solver options with a struct gsGalerkinOptions
//struct gsGalerkinOptions / gsGalerkinConfig
//{
//    
//}

template<class T>
class gsGalerkinMethod
{
public:

    /// Default empty constructor
    gsGalerkinMethod() { }
    
    //  Single Patch or multipatch with common basis
    gsGalerkinMethod( const gsBVProblem<T> & bvp, const gsBasis<T> & basis,
                      dirichlet::strategy dstrategy = dirichlet::elimination,
                      iFace::strategy istrategy = iFace::glue        )
        : m_bvp(&bvp), m_assembler(NULL), m_dofMapper(NULL)
        {
            const gsMultiPatch<T> & mp = m_bvp->patches();
            mp.checkConsistency();
            
            // Collect the given basis, to be used for all patches
            // \todo adjust  bases domain to match geometry
            m_bases.push_back( basis.clone() );
            
            // Initialize assember and dof mapper
            m_dirStrategy       = dstrategy;
            m_interfaceStrategy = istrategy;
            initAssembler();
            initDofMapper();
        }

/*
    //  Single Patch or multipatch using the geometry basis
    gsGalerkinMethod( const gsBVProblem<T> & bvp,
                      gsDirichletStrategy dstrategy = elimination,
                      gsInterfaceStrategy istrategy = glue       )
        : m_bvp(&bvp), m_assembler(NULL), m_dofMapper(NULL)
        {
            const gsMultiPatch<T> & mp = m_bvp->patches();
            mp.checkConsistency();
            
            // Collect & clone the basis of the geometry
            // \todo adjust  bases domain to match geometry
            m_bases.resize( mp.nPatches() );
            for (size_t i=0; i< mp.nPatches(); ++i ) 
                m_bases[i] = mp.patch(i).basis().clone();
            
            // Initialize assember and dof mapper
            m_dirStrategy = dstrategy;
            m_interfaceStrategy = istrategy;
            initAssembler();
            initDofMapper();
        }
*/

    // Multipatch problem
    gsGalerkinMethod( const gsBVProblem<T> & bvp,
                      std::vector<gsBasis<T> *> bases,
                      dirichlet::strategy dstrategy = dirichlet::elimination,
                      iFace::strategy istrategy = iFace::glue       )
        : m_bvp(&bvp), m_assembler(NULL), m_dofMapper(NULL)
        {
            m_bases.resize(bases.size() );
            cloneAll( bases.begin(), bases.end(), m_bases.begin() );

            GISMO_ASSERT( (int)m_bases.size() == m_bvp->nPatches(), "Error" );

            const gsMultiPatch<T> & mp = m_bvp->patches();
            mp.checkConsistency();

            // Initialize assember and dof mapper
            m_dirStrategy       = dstrategy;
            m_interfaceStrategy = istrategy;
            initAssembler();
            initDofMapper();
        }

    // Multipatch problem
    gsGalerkinMethod( const gsBVProblem<T> & bvp,
                      const gsMultiBasis<T>& bases,
                      dirichlet::strategy dstrategy = dirichlet::elimination,
                      iFace::strategy istrategy = iFace::glue       )
        : m_bvp(&bvp), m_assembler(NULL), m_dofMapper(NULL)
        {
            m_bases.resize(bases.nBases() );
            cloneAll( bases.begin(), bases.end(), m_bases.begin() );

            GISMO_ASSERT( (int)m_bases.size() == m_bvp->nPatches(), "Error" );

            const gsMultiPatch<T> & mp = m_bvp->patches();
            mp.checkConsistency();

            // Initialize assember and dof mapper
            m_dirStrategy       = dstrategy;
            m_interfaceStrategy = istrategy;
            initAssembler();
            initDofMapper();
        }

    //Destructor
    ~gsGalerkinMethod();
    
public:
    
    /// Assemble the linear system.
    void assemble();

    /// Solve the PDE and return the solution field.
    gsField<T> * solve();

    /// Assemble and return the mass matrix.
    gsSparseMatrix<T>* assembleMass() const;

    /// Return the gsBVProblem describing the boundary value problem.
    const gsBVProblem<T>& bvp() const           { return *m_bvp; }

    /// Return the gsMultiPatch describing the computational domain.
    const gsMultiPatch<T>& patches() const      { return m_bvp->patches(); }

    /// Return the basis for the \a i-th patch.
    const gsBasis<T>& basis(int i = 0) const    { return *m_bases[i]; }

    /// Return the DOF mapper.
    const gsDofMapper& dofMapper() const     { return *m_dofMapper; }

    /// Return the values of the eliminated Dirichlet dofs, if any.
    const gsMatrix<T>& boundaryValues() const   { return m_ddof; }

    /// Total number of degrees of freedom, including Dirichlet dofs.
    int dofs() const                            { return m_dofs; }

    /// Number of degrees of freedom which we have to solve for, i.e., size of the linear system.
    int freeDofs() const                        { return m_idofs; }

    void setAssembler( gsAssembler<T>* assblr)
    {
        if ( m_assembler )
            delete m_assembler;
        m_assembler = assblr;
    }

    /// Return the assembler.
          gsAssembler<T>& assembler()          { return *m_assembler; }
    /// Return the assembler.
    const gsAssembler<T>& assembler() const    { return *m_assembler; }

    /// Return the assembled linear system.
    const gsSparseSystem<T>& linearSystem() const { return m_system; }

    /// @brief Given a vector of the free dof coefficients, complete it with
    /// the Dirichlet data and return a solution field.
    gsField<T> * reconstructSolution(const gsVector<T>& data) const;

    /// @brief Given a vector of the free dof coefficients, complete it with
    /// the Dirichlet data and return the solution on a single patch \a p.
    gsGeometry<T> * reconstructPatchSolution(const gsVector<T>& data, int p) const;

    /// Set the strategy for dealing with Dirichlet dofs.
    void setDirichletStrategy( dirichlet::strategy strategy)
    { 
        if ( m_dirStrategy != strategy )
        {
            m_dirStrategy = strategy; 
            initDofMapper();
        }
    }

    /// Set the strategy for dealing with patch interfaces.
    void setInterfaceStrategy( iFace::strategy strategy)
    { 
        if ( m_interfaceStrategy != strategy )
        {
            m_interfaceStrategy = strategy; 
            initDofMapper();
        }
    }

    // DirichletStrategy getDirichletStrategy()                { return m_dirStrategy; }
    
    int numIterations()
    {
      return ls_iterations;
    }
    
private:
    // disable copying
    gsGalerkinMethod(const gsGalerkinMethod& other);
    gsGalerkinMethod& operator=(const gsGalerkinMethod& other);

    void initAssembler();
    void initDofMapper();
    // Computes the Dirichlet Boundary function by interpolation.
    void computeDirichletDofs();

    void assembleMultipatch();

// Problem & solver data
private:

    //BV Problem
    const gsBVProblem<T> * m_bvp;

    // Discretization data
    std::vector<gsBasis<T> *> m_bases ;

    // Assembly data
    gsAssembler<T> * m_assembler;
    gsDofMapper * m_dofMapper;

    // Strategy for dirichlet Dofs
    dirichlet::strategy m_dirStrategy;

    // Strategy for treating interfaces
    iFace::strategy m_interfaceStrategy;
    
    // data for Dirichlet dofs by Elimination
    gsMatrix<T> m_ddof;

    // Linear solving  data
    // CG
    gsSparseSystem<T> m_system;
    
    // number of degrees of freedom
    int m_dofs;
    int m_idofs;

    // Linear solver log
    int ls_iterations;
    //T ls_error;

}; // class gsClass


//////////////////////////////////////////////////
//////////////////////////////////////////////////


} // namespace gismo


#ifndef GISMO_HEADERS_ONLY
#include GISMO_HPP_HEADER(gsGalerkinMethod.hpp)
#endif
