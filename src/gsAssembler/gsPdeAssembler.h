/** @file gsPdeAssembler.h

    @brief Provides base interface for linear PDE's.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris, J. Sogn
*/

#pragma once

#include <gsPde/gsBoundaryConditions.h>
#include <gsCore/gsMultiBasis.h>

#include <gsCore/gsBasisEvaluator.h>
#include <gsAssembler/gsGaussRule.h>


namespace gismo
{    

/** @brief
    A generic class for Pde solvers.

    \param[in] patches is a gsMultiPatch object describing the geometry.
    \param[in] bconditions is vector of gsBoundaryConditions for each unknown.
    \param[in] bases is a vector of vectors of gsBasis.
    The base vector is for each unknown (example Stokes equation: @b u and p).
    The child vector is for each components of the unknown.

    It sets up an assembler and assembles the system patch wise and combines
    the patch-local stiffness matrices into a global system by various methods
    (see gismo::gsInterfaceStrategy). It can also enforce Dirichlet boundary
    conditions in various ways (see gismo::gsDirichletStrategy).

    By calling solve(), the linear system is assembled and solved, and a
    gsField containing the computed solution is returned.

*/


template<class T>
class gsPdeAssembler
{

public:
    /// Default empty constructors
    gsPdeAssembler() { }

    ///Constructor without basis
    gsPdeAssembler( gsMultiPatch<T> const & patches,
                 gsBoundaryConditions<T> const & bconditions )
        : m_patches(patches)
        {
            m_bconditions.push_back(new gsBoundaryConditions<T>(bconditions));
            m_patches.checkConsistency();  
        }

    ///Constructor without basis.
    gsPdeAssembler( gsMultiPatch<T> const & patches,
                 std::vector<gsBoundaryConditions<T>*> const & bconditions)
        : m_patches(patches)
        {
            for (typename std::vector<gsBoundaryConditions<T>*>::size_type j = 0;
                 j < bconditions.size(); j++)
                {
                if (bconditions[j])
                    m_bconditions.push_back(new gsBoundaryConditions<T>(*bconditions[j]));
                //If NULL point push back NULL pointer
                else
                    m_bconditions.push_back(NULL);
                }
            m_patches.checkConsistency();
        }

    /// Constructor
    gsPdeAssembler( gsMultiPatch<T> const & patches,
                 gsBoundaryConditions<T> const & bconditions,
                 std::vector< gsBasis<T>*> const & bases )
        : m_patches(patches)
        {
            //m_bases.push_back(bases);
            m_bases.push_back( gsMultiBasis<T>(bases,patches));
            m_bconditions.push_back(new gsBoundaryConditions<T>(bconditions));
            m_patches.checkConsistency();
        }

    /// Constructor
    gsPdeAssembler( gsMultiPatch<T> const & patches,
                 gsBoundaryConditions<T> const & bconditions,
                 std::vector< std::vector< gsBasis<T>*> > const & bases )
        : m_patches(patches)
        {
            for (unsigned k = 0; k < bases.size(); ++k)
            {
                m_bases.push_back( gsMultiBasis<T>(bases[k],patches));
            }
            m_bconditions.push_back(new gsBoundaryConditions<T>(bconditions));
            m_patches.checkConsistency();
        }

    /// Constructor
    gsPdeAssembler( gsMultiPatch<T> const & patches,
                 std::vector<gsBoundaryConditions<T>*> const & bconditions,
                 std::vector< gsBasis<T>*> const & bases )
        : m_patches(patches)
        {
            //m_bases.push_back(bases);
            m_bases.push_back( gsMultiBasis<T>(bases,patches));
            for (typename std::vector<gsBoundaryConditions<T>*>::size_type j = 0;
                 j < bconditions.size(); j++)
            {
                m_bconditions.push_back(new gsBoundaryConditions<T>(*bconditions[j]));
            }
            m_patches.checkConsistency();
        }

    /// General constructor
    gsPdeAssembler( gsMultiPatch<T> const & patches,
                 std::vector<gsBoundaryConditions<T>*> const & bconditions,
                 std::vector< std::vector< gsBasis<T>*> > const & bases )
        : m_patches(patches)//, m_bases(bases)
        {
            for (unsigned k = 0; k < bases.size(); ++k)
            {
                m_bases.push_back( gsMultiBasis<T>(bases[k],patches));
            }
            for (typename std::vector<gsBoundaryConditions<T>*>::size_type j = 0; j < bconditions.size(); ++j)
            {
                if (bconditions[j])
                    m_bconditions.push_back(new gsBoundaryConditions<T>(*bconditions[j]));
                //If NULL point push back NULL pointer
                else
                    m_bconditions.push_back(NULL);
            }
            m_patches.checkConsistency();
        }

    /// General constructor with gsMultiBasis vector
    gsPdeAssembler( gsMultiPatch<T> const & patches,
                 std::vector<gsBoundaryConditions<T>*> const & bconditions,
                 std::vector< gsMultiBasis<T> *>  const & bases )
        : m_patches(patches)//, m_bases(bases)
        {
            for (unsigned k = 0; k < bases.size(); ++k)
            {
                m_bases.push_back(gsMultiBasis<T>(*bases[k]));
            }
            for (typename std::vector<gsBoundaryConditions<T>*>::size_type j = 0; j < bconditions.size(); ++j)
            {
                if (bconditions[j])
                    m_bconditions.push_back(new gsBoundaryConditions<T>(*bconditions[j]));
                //If NULL point push back NULL pointer
                else
                    m_bconditions.push_back(NULL);
            }
            m_patches.checkConsistency();
        }

    /// Constructor with gsMultiBasis vector
    gsPdeAssembler( gsMultiPatch<T> const & patches,
                 gsBoundaryConditions<T> const & bconditions,
                 std::vector< gsMultiBasis<T>* > const & bases )
        : m_patches(patches)
        {
            for (unsigned k = 0; k < bases.size(); ++k)
            {
                m_bases.push_back( gsMultiBasis<T>(*bases[k]));
            }
            m_bconditions.push_back(new gsBoundaryConditions<T>(bconditions));
            m_patches.checkConsistency();
        }

    /// General constructor with gsMultiBasis
    gsPdeAssembler( gsMultiPatch<T> const & patches,
                 std::vector<gsBoundaryConditions<T>*> const & bconditions,
                 gsMultiBasis<T> const & bases )
        : m_patches(patches)//, m_bases(bases)
        {
            m_bases.push_back(gsMultiBasis<T>(bases));

            for (typename std::vector<gsBoundaryConditions<T>*>::size_type j = 0; j < bconditions.size(); ++j)
            {
                if (bconditions[j])
                    m_bconditions.push_back(new gsBoundaryConditions<T>(*bconditions[j]));
                //If NULL point push back NULL pointer
                else
                    m_bconditions.push_back(NULL);
            }
            m_patches.checkConsistency();
        }

    /// Constructor with gsMultiBasis
    gsPdeAssembler( gsMultiPatch<T> const & patches,
                 gsBoundaryConditions<T> const & bconditions,
                 gsMultiBasis<T> const & bases )
        : m_patches(patches)
        {
            m_bases.push_back( gsMultiBasis<T>(bases));
            m_bconditions.push_back(new gsBoundaryConditions<T>(bconditions));
            m_patches.checkConsistency();
        }

   ~gsPdeAssembler()
    {
        freeAll( m_dofMapper);
        freeAll( m_solutions);
        freeAll(m_bconditions);
    }

    /// Clear the matrix and right-hand side.
    void cleanSystem()
    {
        m_matrix.resize(0,0);
        m_rhs.resize(0,0);
    }


    /// Prints the object as a string.
    virtual std::ostream &print(std::ostream &os) const
    {
        GISMO_NO_IMPLEMENTATION
    }


public:
    
    /// Initializes: To be overwritten in derived classes.
    /// This will be called in solve()
    ///
    /// Typically:
    /// Computes any fixed DoFs, initializes some internal members,
    /// checks consistency
    virtual void initialize() = 0;

    /// @brief Assembles the final system with all boundary conditions contained.
    /// \a m_matrix and \a m_rhs are filled.
    virtual void assemble() = 0;

    /// @brief Solves the linear system and fills up \a m_sysSolution
    virtual void solveSystem() = 0;

    // Reconstruct the solution to m_system
    /// @brief Creates the solution fields from the given solution vector(s)
    void reconstructSolution(bool hasMatrixRhs = false);

    /// @brief Solve the PDE, e.i initialize, assemble, solve and reconstruct solution.
    virtual void solve()
    {
        initialize ();
        assemble   ();
        solveSystem();
        reconstructSolution();
    }

    /// @brief Solve the PDE, e.i initialize, assemble, solve and reconstruct solution.
    virtual void nextTimeStep()
    {
        solveSystem();
    }

    /// @brief Returns the left-hand side global matrix
    const gsSparseMatrix<T> & systemMatrix() const { return m_matrix; }

    /// @brief Returns the vector right-hand side as a matrix
    /// ( multiple right hand sides possible ).
    const gsMatrix<T> & systemRhs() const { return m_rhs; }

    /// @brief Returns the solution of the linear system
    /// ( multiple columns are possible ).
    const gsMatrix<T> &  systemSolution() const { return m_sysSolution; }

    /// @brief Gives the solution. Usefull when solving is done outside of this class
    void setSolution(const gsMatrix<T> &solution) {m_sysSolution = solution;}

//    /// @brief Accessor for left-hand side global matrix.
//    ///
//    /// As of now, no checks are performed, so the user has to take care that
//    /// size and format of the input match the other problem data.
//    void setSystemMatrix( const gsSparseMatrix<T> & in_matrix ) { m_matrix = in_matrix; }
//    // Not passed by reference on purpose, so that the matrix outside of
//    // gsPdeAssembler will not be modified.

//    /// @brief Accessor for the vector right-hand side as a matrix.
//    ///
//    /// As of now, no checks are performed, so the user has to take care that
//    /// size and format of the input match the other problem data.
//    void setSystemRhs( const gsMatrix<T> in_Rhs ) { m_rhs = in_Rhs; }
//    // Not passed by reference on purpose, so that the matrix outside of
//    // gsPdeAssembler will not be modified.

//Accessors
public:

    /// Return patch number \em i.
    gsGeometry<T> & patch(unsigned i) const {return m_patches.patch(i); }

    /// Return the solution field for component \em i.
    const gsField<T> & solution(unsigned i = 0) const   {return *m_solutions[i]; }

    /// Return the vector of solution fields.
    const std::vector< gsField<T> *> & solutionVector() const {return m_solutions; }

    /// Set current reconstructed solution \todo this member is temporary, to be removed
    void setSolution(gsField<T> * sol) { m_solutions.push_back(sol); }

    /// Return the vector of fixed DoFs.
    const std::vector<gsMatrix<T> > & fixedValues() const {return m_fixedDofs; }

    /// Return the DOF mapper for unknown \em i.
    const gsDofMapper& dofMapper(unsigned i = 0) const     { return *m_dofMapper[i]; }

    /// Return the vector DOF mapper.
    const std::vector< gsDofMapper *> & dofMapperVector() const {return m_dofMapper;}

    /// Return the multipatch.
    const gsMultiPatch<T> & patches() const { return m_patches; }

    /// Returns the dimension of the unknown field for unknown \em i.
    int fieldDim(index_t i = 0) { return this->m_unknownDim[i]; }

    /// Set whether the PDE is (or should be treated as) symmetric or not.
    void setPdeSymmetry( const bool symm ) {m_isSymmetric = symm; }

    /// Return the vector of boundary conditions
    const std::vector<gsBoundaryConditions<T> *> & boundaryConditionsVector() const {return m_bconditions;}

    /// Return the boundary condition of the i-th unknown
    const gsBoundaryConditions<T> & boundaryConditions(unsigned i = 0) const  {return *m_bconditions[i];}


// Helpers
protected:

    /// Setup the degrees of freedom, interface matching etc
    void initDofMapper(iFace::strategy interfaceStrategy,
                       dirichlet::strategy dirStrategy, bool hasMatrixRhs = false);

    /// @brief Computes the Dirichlet Boundary function by interpolation.
    virtual void computeDirichletDofs(dirichlet::strategy dirStrategy);

    /// @brief Computes the Dirichlet DoFs for unknown \a unk, in case
    /// of elimination. NOT IMPLEMENTED
    void eliminateDirichletDofs(int unk)
    {
        // computeDirichletDofs()
    }
    
    //Handling Neumann boundary conditions
    /// \brief Add contribution of Neumann boundary condition to the linear system.
    ///
    /// Add contribution of Neumann boundary condition to the right hand side of
    /// the linear system. More precisely it adds:
    /// \f[
    /// \int_{\Gamma_s} f\, B_i dS
    /// \f]
    /// Where \f$ B_i \f$ are the elements of basis \a B.
    /// \param B is a boundary basis,
    /// \param f is the (Neumann) function
    /// \param s is the side of the boundary
    /// \param mapper       contains the global Dof map for the linear system
    /// \param blockMapper       contains information the block structure of the
    /// right hand side. Where the each element in \a blockMapper is the index
    /// for the first element in a block. Default is no block structure.
    void boundaryNeumann( const gsBasis<T> & B,
                          const int patchIndex,
                          const gsFunction<T> & f,
                          const boxSide s,
                          const gsDofMapper& mapper,
                          gsVector<index_t> blockMapper = gsVector<index_t>::Zero(1));

    /// \brief Add contribution of Neumann boundary condition to the linear system for one basis for each components.
    void boundaryNeumann( std::vector< gsBasis<T> *>  const & B,
                          const int patchIndex,
                          const gsFunction<T> & f,
                          const boxSide s,
                          std::vector< gsDofMapper *> const & mapper,
                          gsVector<index_t> blockMapper = gsVector<index_t>::Zero(1));



    /// \brief Add contribution of Nitsche Dirichlet boundary to global matrix.
    ///
    /// Weakly emposes the Dirichlet boundary condition with the Nitsch method.
    /// It adds the following terms to the global stiffness matrix:
    /// \f[
    /// - \int_{\Gamma_s} (\kappa \nabla B_j \cdot n) \, B_i \,dS - \int_{\Gamma_s} \kappa B_j \, (\nabla B_i \cdot n)\,dS + \int_{\Gamma_s}  \kappa \mu(h) \, B_j \,B_i \,dS
    /// \f]
    /// It adds the following terms to the right hand side:
    /// \f[
    /// - \int_{\Gamma_s} \kappa \, f \, (\nabla B_i \cdot n) \,dS + \int_{\Gamma_s}  \kappa \mu(h) \, f \,B_i \,dS
    /// \f]
    /// \param B is a boundary basis,
    /// \param f is the Dirichlet function
    /// \param s is the side of the boundary
    /// \param mapper contains the global Dof map for the linear system
    /// \param blockMapper       Is a \f$bn \times 2\f$ matrix containing information
    /// the block structure of the global stiffness matrix and right hand side.
    /// Where row \a i in \a blockMapper is the index for the first
    /// element in block \a i. Default is no block structure.
    /// \param kappa is a constant (should be a function)
    void boundaryNitsche( const gsBasis<T> & B,
                          const int patchIndex,
                          const gsFunction<T> & f,
                          const boxSide s,
                          const gsDofMapper& mapper,
                          gsMatrix<index_t> blockMapper = gsMatrix<index_t>::Zero(1,2),
                          const T kappa = 1);

    /// \brief Add contribution of Nitsche Dirichlet boundary condition to the linear system for one basis for each components.
    void boundaryNitsche( std::vector< gsBasis<T> *>  const & B,
                          const int patchIndex,
                          const gsFunction<T> & f,
                          const boxSide s,
                          std::vector< gsDofMapper *> const & mapper,
                          gsMatrix<index_t> blockMapper = gsMatrix<index_t>::Zero(1,2),
                          const T kappa = 1);

    //JS2 should this be private?
    //since reconstructSolution() is the only one calling it?
    /// @brief Creates the solution fields from the given solution
    /// vector for unknown \a unk on patch \a p.
    virtual gsFunction<T> * reconstructPatchSolution(index_t unk, int p, bool hasMatrixRhs = false) const;


    // do not use yet
void remapToFullSystem (std::vector<gsDofMapper*> &input);
void reconstructPatchCoefficients (
    const gsDofMapper &mapper,
          index_t      patch_id,
    const gsMatrix<T> &solution,
    const gsMatrix<T> &eliminated_values,
          gsMatrix<T> &result
    );

protected:
    gsPde<T> * m_pde;

    // The multipatch domain
    gsMultiPatch<T> m_patches;

    // The boundary conditions, one set of conditions for each unknown
    std::vector<gsBoundaryConditions<T> *> m_bconditions;

    // The discretization bases corresponding to \a m_patches and to
    // the number of solution fields that are to be computed
    // m_bases[u][k]: basis for unknown (not for each component) u on patch k
    //std::vector< gsMultiBasis<T> > m_bases;
     std::vector<
         gsMultiBasis<T>
         //std::vector< gsBasis<T>* >  // this should be a "gsMultiBasis" object
         > m_bases;


    // *** Internal *** 

    // To indicate whether the PDE is symmetric or not,
    // or should be treated like a symmetric problem or not.
    bool m_isSymmetric;

    // Global matrix
    gsSparseMatrix<T> m_matrix;

    // Right-hand side ( multiple right hand sides possible )
    gsMatrix<T>       m_rhs;

    // Solution of the linear system
    gsMatrix<T> m_sysSolution;

    // Values for the degrees of freedom that are known in advance and
    // are not part of the system (eg. Dirichlet data)
    // Each matrix corresponds to an unknown, each column is for each component
    // Each row is for each Dirichlet data point
    std::vector<gsMatrix<T> > m_fixedDofs;

    // The Dof mapper is used to map patch-local DoFs to the global DoFs
    // One for each unknown, one for each patch
    // m_dofMapper[u]: DoF Mapper for unknown (and component) u
    std::vector< gsDofMapper *>  m_dofMapper; // TO DERIVED ?

    // Global matrix block structure --> offsets vector
    // gsVector<index_t> m_matrixOffset;
    // gsBlockStructure<T> m_blockStructure; // TO DERIVED ?

    // *** Information *** 

    // number of degrees of freedom
    int m_dofs;
    int m_idofs;

    // Saves the dimension of the unknown(s)
    // i.e how many components each unknown has.
    // it's size is the number of unknowns
    gsVector<int> m_unknownDim;

    // *** Outputs ***

    // Solution field(s)
    std::vector< gsField<T> *> m_solutions;

}; // class gsPdeAssembler


/// Print (as string) operator to be used by all derived classes
template<class T>
std::ostream &operator<<(std::ostream &os, const gsPdeAssembler<T>& b)
{return b.print(os); }


} // namespace gismo


#include "gsPdeAssembler.hpp"
