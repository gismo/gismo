/** @file gsStokesAssembler.h

    @brief Assembler and solver for the Stokes problem

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): J. Sogn
*/
#pragma once

#include <gsAssembler/gsPdeAssembler.h>
#include <gsAssembler/gsAssemblerOptions.h>
#include <gsPde/gsBoundaryConditions.h>
#include <gsCore/gsDivConSolution.h>

//Testing gsBasisEvaluator
#include <gsCore/gsBasisEvaluator.h>
#include <gsCore/gsBasisEvaluator.hpp>
#include <gsAssembler/gsGaussRule.h>


namespace gismo
{

/** @brief
    Implementation of an Stokes solver.

    The Stokes problem: \f{eqnarray*}{
                        -\nu\Delta\mathbf{u} - \nabla p &=&\mathbf{f} \quad \text{in} \quad \Omega,\\
                        \nabla \cdot \mathbf{u} &=& 0 \quad \text{in} \quad \Omega,\\
                        \mathbf{u} &=& \mathbf{g} \quad \text{on} \quad \Gamma_D,\\
                        \nu\nabla\mathbf{u}\cdot\mathbf{n} + p \mathbf{n} &=& \mathbf{h} \quad \text{on} \quad \Gamma_N.
                        \f}
    The Neumann condition \f$ \mathbf{h}\f$ is set by applying Neumann for the velocity.
    \note The sign pressure of the pressure is fliped to obtain a symmetric linear system!
    TODO add:
            the spaces for the weak formulation we using this solver
            inf-sup warning


    See gsStokesAssembler for general description of input.
    \param[in] patches is a gsMultiPatch object describing the geometry.
    TODO UPDATE\param[in] bconditions is a gsBoundaryConditions object.
    TODO UPDATE\param[in] bases is a vector of gsBasis for each components
    of the unknown of \f$\mathbf{u}\f$.
    \param[in] rhs is the right-hand side (source term) of the Stokes problem,
    \f$\mathbf{f}\f$.

    It sets up an assembler and assembles the system patch wise and combines
    the patch-local stiffness matrices into a global system by various methods
    (see gismo::iFace::strategy). It can also enforce Dirichlet boundary
    conditions in various ways (see gismo::gsDirichletStrategy).

    By calling gsStokesAssembler::solve(), the linear system is assembled and solved, and
    a std::vector of gsField containing the computed solution (velocity and pressure)
    is returned.

    \todo Reserve space is sparce matrix
*/


template<class T>
class gsStokesAssembler : public gsPdeAssembler<T>
{
public:
    /// Default empty constructor
    gsStokesAssembler () : gsPdeAssembler<T>() { }


    ///Constructor without basis and only BC one unknown (The velocity).
    gsStokesAssembler( gsMultiPatch<T> const & patches,
                    gsBoundaryConditions<T> const & bconditions,
                    const gsFunction<T> & rhs, T nu = 1.0)
         : gsPdeAssembler<T>(patches, bconditions), m_nu(nu)
     {
         // Add a NULL pointer for the pressure BC
         gsBoundaryConditions<T> * ptr_null = NULL;
         m_bconditions.push_back(ptr_null);

         m_unknownDim.resize(2);
         m_unknownDim[0] = rhs.targetDim(); //Target dimension for velocity
         m_unknownDim[1] = 1; //Target dimension for pressure
         m_dirStrategy = dirichlet::elimination;
         m_interfaceStrategy = iFace::glue;
         m_geoTrans = INVERSE_COMPOSITION;
         m_rhs_function = rhs.clone();
         // Collect & clone the basis of the geometry
         // \todo adjust  bases domain to match geometry
         m_bases.resize(2);
         m_bases[0].resize( m_patches.nPatches() );
         m_bases[1].resize( m_patches.nPatches() );
         for (size_t i=0; i< m_patches.nPatches(); ++i )
             {
             m_bases[0][i] = m_patches.patch(i).basis().clone();
             // Increase the basis degree for the velocity to satisfy inf-sup condition
             m_bases[0][i]->degreeElevate();
             m_bases[1][i] = m_patches.patch(i).basis().clone();
             }
     }


    ///Constructor without basis.
    gsStokesAssembler( gsMultiPatch<T> const & patches,
                    std::vector<gsBoundaryConditions<T>*> const & bconditions,
                    const gsFunction<T> & rhs, T nu = 1.0 )
        : gsPdeAssembler<T>(patches, bconditions), m_nu(nu)
    {
        m_unknownDim.resize(2);
        m_unknownDim[0] = rhs.targetDim(); //Target dimension for velocity
        m_unknownDim[1] = 1; //Target dimension for pressure
        m_dirStrategy = dirichlet::elimination;
        m_interfaceStrategy = iFace::glue;
        m_geoTrans = INVERSE_COMPOSITION;
        m_rhs_function = rhs.clone();
        // Collect & clone the basis of the geometry
        // \todo adjust  bases domain to match geometry
        m_bases.resize(2);
        m_bases[0].resize( m_patches.nPatches() );
        m_bases[1].resize( m_patches.nPatches() );
        for (size_t i=0; i< m_patches.nPatches(); ++i )
            {
            m_bases[0][i] = m_patches.patch(i).basis().clone();
            // Increase the basis degree for the velocity to satisfy inf-sup condition
            m_bases[0][i]->degreeElevate();
            m_bases[1][i] = m_patches.patch(i).basis().clone();
            }
    }


    ///Constructor with only BC one unknown (The velocity).
    gsStokesAssembler( gsMultiPatch<T> const & patches,
                    gsBoundaryConditions<T> const & bconditions,
                    std::vector< std::vector< gsBasis<T>* > > const & bases,
                    const gsFunction<T> & rhs, T nu = 1.0 )
        : gsPdeAssembler<T>(patches, bconditions, bases), m_nu(nu)
    {
        // Add a NULL pointer for the pressure BC
        gsBoundaryConditions<T> * ptr_null = NULL;
        m_bconditions.push_back(ptr_null);

        m_unknownDim.resize(2);
        m_unknownDim[0] = rhs.targetDim(); //Target dimension for velocity
        m_unknownDim[1] = 1; //Target dimension for pressure

        m_dirStrategy = dirichlet::elimination;
        m_interfaceStrategy = iFace::glue;
        m_geoTrans = INVERSE_COMPOSITION;
        m_rhs_function = rhs.clone();
    }



    ///General Constructor
    gsStokesAssembler( gsMultiPatch<T> const & patches,
                    std::vector<gsBoundaryConditions<T>*> const & bconditions,
                    std::vector< std::vector< gsBasis<T>* > > const & bases,
                    const gsFunction<T> & rhs, T nu = 1.0 )
        : gsPdeAssembler<T>(patches, bconditions, bases), m_nu(nu)
    {
        m_unknownDim.resize(2);
        m_unknownDim[0] = rhs.targetDim(); //Target dimension for velocity
        m_unknownDim[1] = 1; //Target dimension for pressure

        m_dirStrategy = dirichlet::elimination;
        m_interfaceStrategy = iFace::glue;
        m_geoTrans = INVERSE_COMPOSITION;
        m_rhs_function = rhs.clone();
    }

    ///Constructor with gsMultiBasis
    gsStokesAssembler( gsMultiPatch<T> const & patches,
                    std::vector<gsBoundaryConditions<T>*> const & bconditions,
                    std::vector< gsMultiBasis<T>* > const & bases,
                    const gsFunction<T> & rhs, T nu = 1.0 )
        : gsPdeAssembler<T>(patches, bconditions, bases), m_nu(nu)
    {
        m_unknownDim.resize(2);
        m_unknownDim[0] = rhs.targetDim(); //Target dimension for velocity
        m_unknownDim[1] = 1; //Target dimension for pressure

        m_dirStrategy = dirichlet::elimination;
        m_interfaceStrategy = iFace::glue;
        m_geoTrans = INVERSE_COMPOSITION;
        m_rhs_function = rhs.clone();
    }

    ///Constructor with only BC one unknown (The velocity).
    gsStokesAssembler( gsMultiPatch<T> const & patches,
                    gsBoundaryConditions<T> const & bconditions,
                    std::vector< gsMultiBasis<T>* > const & bases,
                    const gsFunction<T> & rhs, T nu = 1.0 )
        : gsPdeAssembler<T>(patches, bconditions, bases), m_nu(nu)
    {
        // Add a NULL pointer for the pressure BC
        gsBoundaryConditions<T> * ptr_null = NULL;
        m_bconditions.push_back(ptr_null);

        m_unknownDim.resize(2);
        m_unknownDim[0] = rhs.targetDim(); //Target dimension for velocity
        m_unknownDim[1] = 1; //Target dimension for pressure

        m_dirStrategy = dirichlet::elimination;
        m_interfaceStrategy = iFace::glue;
        m_geoTrans = INVERSE_COMPOSITION;
        m_rhs_function = rhs.clone();
    }


   ~gsStokesAssembler()
    {
        delete m_rhs_function;
    }



    /// Set the strategy for dealing with Dirichlet DoFs.\n
    /// \param strategy options are:\n
    /// Eliminate the Dirichlet DoFs from the linear system.\n
    /// \em elimination = 1 (default)\n\n
    /// Keep the Dirichlet DoFs and enforce the boundary
    /// condition weakly by a penalty term.\n
    /// \em nitsche = 2
    void setDirichletStrategy(dirichlet::strategy strategy)
    { m_dirStrategy = strategy; }

    /// Set the strategy for dealing with patch interface.\n
    /// \param strategy options are:\n
    /// Glue patches together by merging DoFs across an interface into one.
    /// This only works for conforming interfaces.\n
    /// \em glue = 1 (default)\n\n
    /// Use discontinuous Galerkin-like coupling between adjacent patches.\n
    /// \em dg = 2
    void setInterfaceStrategy(iFace::strategy strategy)
    { m_interfaceStrategy = strategy; }

    /// Set the which geometrical transformation should be used.\n
    /// \param strategy options are:\n
    /// Standard Inverse composition (Normal for IGA)
    /// \em INVERSE_COMPOSITION (default)\n\n
    /// Divergence preserving transformation (also called "Contravariant Piola transformation")
    /// \em DIV_CONFORMING \n\n
    /// No transformation
    /// \em NO_TRANSFORMATION\n\n
    void setGeometryTransformation(ValueTransformationType strategy)
    { m_geoTrans = strategy; }


    /// Set the kinematic viscosity (default is \f$ \nu = 1.0\f$)
    void setnu(T nu)
    { m_nu = nu; }


    //Virtual function
    public:

    /// @brief Initializes internal members
    /// Initializes: m_matrix, m_matrixBlocks, m_rhs, m_rhsBlocks, m_dofs, m_idofs and m_dofMapper
    void initialize();


    /// @brief Assembles the final system with all boundary conditions contained
    void assemble();



    // Solves the linear system and fills up \a m_sysSolution
    /**
     * @brief solveSystem solves the linear system of equations Ax = b
     * and fills up \a m_sysSolution.
     */
    void solveSystem();



//Accessors
public:

    /// @brief Return the right-hand side function \f$\mathbf{f}\f$
    const gsFunction<T> & rhs_function() const  { return *m_rhs_function; }

public:
    /// @brief Creates the solution fields from the given solution
    /// vector for unknown \a unk on patch \a p.
    gsFunction<T> * reconstructPatchSolution(index_t unk, int p, bool hasMatrixRhs = false) const;



private:

    /// Assembles on a single patch \a patchIndex
    void assemblePatch(int patchIndex=0);


    /// Prints the object as a string.
    std::ostream &print(std::ostream &os) const
    {   os <<"Stokes problem: -\u03BD\u0394u + \u2207p = f, \n";
        os <<"                        \u2207 \u00B7u = 0, with:\n";
        os << "* f: "<< *m_rhs_function <<"\n";
        os << "* Domain: "<< m_patches ;
        os << "* \u03BD: "<< m_nu <<"\n";
        if (m_bconditions[0]) //if bconditions[0] is not an Null pointer
            os << "Boundary condition for velocity "<< *m_bconditions[0] <<"\n";
        if (m_bconditions[1]) //if bconditions[0] is not an Null pointer
            os << "Boundary condition for pressure "<< *m_bconditions[1] <<"\n";
        return os;
    }

    //JS2: TODO Document!
    void boundaryNeumann( const int patchIndex,
                          const gsFunction<T> & f,
                          const boundary::side s);

    void boundaryNitsche( const int patchIndex,
                          const gsFunction<T> & f,
                          const boundary::side s);


    /// @brief Computes the Dirichlet Boundary function by interpolation.
    void computeDirichletDofs(dirichlet::strategy dirStrategy);



private:
    // Kinematic viscosity
    T m_nu;

    // Right hand side function
    gsFunction<T> * m_rhs_function;
    // Strategy for dealing with Dirichlet DoFs
    dirichlet::strategy m_dirStrategy;
    // Strategy for dealing with patch interface
    iFace::strategy m_interfaceStrategy;
    // Geometrical transformation type
    ValueTransformationType m_geoTrans;



    // Block structure of m_matrix
    typename gsSparseMatrix<T>::BlockView m_matrixBlocks;

    // Block structure of m_matrix
    typename gsMatrix<T>::BlockView m_rhsBlocks;

    // Members from gsPdeAssembler
    using gsPdeAssembler<T>::m_patches;
    using gsPdeAssembler<T>::m_bconditions;
    using gsPdeAssembler<T>::m_bases;

    using gsPdeAssembler<T>::m_matrix;
    using gsPdeAssembler<T>::m_rhs;
    using gsPdeAssembler<T>::m_sysSolution;
    using gsPdeAssembler<T>::m_fixedDofs;
    using gsPdeAssembler<T>::m_dofMapper;

    using gsPdeAssembler<T>::m_dofs;
    using gsPdeAssembler<T>::m_idofs;

    using gsPdeAssembler<T>::m_unknownDim;
    using gsPdeAssembler<T>::m_solutions;


}; // class gsStokesAssembler


} // namespace gismo

#include "gsStokesAssembler.hpp"

