/** @file gsPoissonAssembler.h

    @brief Provides assembler and solver for the Poisson equation, incl. adaptive refinement.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): S. Kleiss A. Mantzaflaris, J. Sogn
*/

#pragma once

#include <gsAssembler/gsPdeAssembler.h>
#include <gsPde/gsBoundaryConditions.h>

namespace gismo
{

/** @brief
    Implementation of an (Vector valued) Poisson solver.

    The Poisson equation: \f$-\Delta\mathbf{u}=\mathbf{f} \f$

    See gsPdeAssembler for general description of input.
    \param[in] patches is a gsMultiPatch object describing the geometry.
    \param[in] bconditions is a gsBoundaryConditions object.
    \param[in] bases is a vector of gsBasis for each components of the unknown
    of \f$\mathbf{u}\f$.
    \param[in] rhs is the right-hand side of the Poisson equation, \f$\mathbf{f}\f$.

    It sets up an assembler and assembles the system patch wise and combines
    the patch-local stiffness matrices into a global system by various methods
    (see gismo::gsInterfaceStrategy). It can also enforce Dirichlet boundary
    conditions in various ways (see gismo::gsDirichletStrategy).

    By calling gsPoissonAssembler::solve(), the linear system is assembled and solved, and a
    gsField containing the computed solution is returned.


*/

template<class T>
class gsPoissonAssembler : public gsPdeAssembler<T>
{
public:
    /// Default empty constructor
    gsPoissonAssembler () : gsPdeAssembler<T>() { }

    //JS2 not tested
    ///Constructor without basis.
    gsPoissonAssembler( gsMultiPatch<T> const & patches,
                      gsBoundaryConditions<T> const & bconditions,
                      const gsFunction<T> & rhs, gsFunction<T> * sol = 0 )
         : gsPdeAssembler<T>(patches, bconditions)
     {
         this->m_unknownDim.resize(1);
         this->m_unknownDim[0] = rhs.targetDim();
         m_dirStrategy = dirichlet::elimination;
         m_isSymmetric = true;
         m_interfaceStrategy = iFace::glue;
         m_rhs_function = rhs.clone();
         // Collect & clone the basis of the geometry
         // \todo adjust  bases domain to match geometry
         this->m_bases.resize(1);
         this->m_bases[0].resize( this->m_patches.nPatches() );
         for (size_t i=0; i< this->m_patches.nPatches(); ++i )
         {
             this->m_bases[0][i] = this->m_patches.patch(i).basis().clone();
         }
     }

    //JS2 not tested
    ///Constructor without basis.
    gsPoissonAssembler( gsMultiPatch<T> const & patches,
                      std::vector<gsBoundaryConditions<T>*> const & bconditions,
                      const gsFunction<T> & rhs, gsFunction<T> * sol = 0 )
         : gsPdeAssembler<T>(patches, bconditions)
     {
         std::cout << "Debugging: I am in gsPoissonAssembler Constructor  PDE:2sv "<< std::endl;
         this->m_unknownDim.resize(1);
         this->m_unknownDim[0] = rhs.targetDim();
         m_dirStrategy = dirichlet::elimination;
         m_isSymmetric = true;
         m_interfaceStrategy = iFace::glue;
         m_rhs_function = rhs.clone();
         // Collect & clone the basis of the geometry
         // \todo adjust  bases domain to match geometry
         this->m_bases.resize(1);
         this->m_bases[0].resize( this->m_patches.nPatches() );
         for (size_t i=0; i< this->m_patches.nPatches(); ++i )
         {
             this->m_bases[0][i] = this->m_patches.patch(i).basis().clone();
         }
         std::cout << "Debugging: I am in gsPoissonAssembler Constructor  End" << std::endl;
     }


    ///Constructor
    gsPoissonAssembler( gsMultiPatch<T> const & patches,
                     gsBoundaryConditions<T> const & bconditions,
                     std::vector< gsBasis<T>* > const & bases,
                     const gsFunction<T> & rhs, gsFunction<T> * sol = 0 )
         : gsPdeAssembler<T>(patches, bconditions, bases)
    {
        this->m_unknownDim.resize(1);
        this->m_unknownDim[0] = rhs.targetDim();
        m_dirStrategy = dirichlet::elimination;
        m_isSymmetric = true;
        m_interfaceStrategy = iFace::glue;
        m_rhs_function = rhs.clone();
    }

    ///Constructor.
    gsPoissonAssembler( gsMultiPatch<T> const & patches,
                     gsBoundaryConditions<T> const & bconditions,
                     std::vector< std::vector< gsBasis<T>* > > const & bases,
                     const gsFunction<T> & rhs, gsFunction<T> * sol = 0 )
        : gsPdeAssembler<T>(patches, bconditions, bases)
    {
        this->m_unknownDim.resize(1);
        this->m_unknownDim[0] = rhs.targetDim();
        m_dirStrategy = dirichlet::elimination;
        m_isSymmetric = true;
        m_interfaceStrategy = iFace::glue;
        m_rhs_function = rhs.clone();
    }

    ///Constructor.
    gsPoissonAssembler( gsMultiPatch<T> const & patches,
                     std::vector<gsBoundaryConditions<T>*> const & bconditions,
                     std::vector< gsBasis<T>* > const & bases,
                     const gsFunction<T> & rhs, gsFunction<T> * sol = 0 )
        : gsPdeAssembler<T>(patches, bconditions, bases)
    {
        this->m_unknownDim.resize(1);
        this->m_unknownDim[0] = rhs.targetDim();
        m_dirStrategy = dirichlet::elimination;
        m_isSymmetric = true;
        m_interfaceStrategy = iFace::glue;
        m_rhs_function = rhs.clone();
    }

    ///Constructor
    gsPoissonAssembler( gsMultiPatch<T> const & patches,
                     std::vector<gsBoundaryConditions<T>*> const & bconditions,
                     std::vector< std::vector< gsBasis<T>* > > const & bases,
                     const gsFunction<T> & rhs, gsFunction<T> * sol = 0 )
        : gsPdeAssembler<T>(patches, bconditions, bases)
    {
        this->m_unknownDim.resize(1);
        this->m_unknownDim[0] = rhs.targetDim();
        m_dirStrategy = dirichlet::elimination;
        m_isSymmetric = true;
        m_interfaceStrategy = iFace::glue;
        m_rhs_function = rhs.clone();
    }

    ///Constructor with gsMultiBasis
    gsPoissonAssembler( gsMultiPatch<T> const & patches,
                     gsBoundaryConditions<T> const & bconditions,
                     gsMultiBasis<T> const & bases,
                     const gsFunction<T> & rhs, gsFunction<T> * sol = 0 )
         : gsPdeAssembler<T>(patches, bconditions, bases)
    {
        this->m_unknownDim.resize(1);
        this->m_unknownDim[0] = rhs.targetDim();
        m_dirStrategy = dirichlet::elimination;
        m_isSymmetric = true;
        m_interfaceStrategy = iFace::glue;
        m_rhs_function = rhs.clone();
    }

    ///Constructor with gsMultiBasis
    gsPoissonAssembler( gsMultiPatch<T> const & patches,
                     gsBoundaryConditions<T> const & bconditions,
                     std::vector< gsMultiBasis<T>* > const & bases,
                     const gsFunction<T> & rhs, gsFunction<T> * sol = 0 )
        : gsPdeAssembler<T>(patches, bconditions, bases)
    {
        this->m_unknownDim.resize(1);
        this->m_unknownDim[0] = rhs.targetDim();
        m_dirStrategy = dirichlet::elimination;
        m_isSymmetric = true;
        m_interfaceStrategy = iFace::glue;
        m_rhs_function = rhs.clone();
    }

    ///Constructor with gsMultiBasis
    gsPoissonAssembler( gsMultiPatch<T> const & patches,
                     std::vector<gsBoundaryConditions<T>*> const & bconditions,
                     gsMultiBasis<T> const & bases,
                     const gsFunction<T> & rhs, gsFunction<T> * sol = 0 )
        : gsPdeAssembler<T>(patches, bconditions, bases)
    {
        this->m_unknownDim.resize(1);
        this->m_unknownDim[0] = rhs.targetDim();
        m_dirStrategy = dirichlet::elimination;
        m_isSymmetric = true;
        m_interfaceStrategy = iFace::glue;
        m_rhs_function = rhs.clone();
    }

    ///Constructor with gsMultiBasis
    gsPoissonAssembler( gsMultiPatch<T> const & patches,
                     std::vector<gsBoundaryConditions<T>*> const & bconditions,
                     std::vector< gsMultiBasis<T>* > const & bases,
                     const gsFunction<T> & rhs, gsFunction<T> * sol = 0 )
        : gsPdeAssembler<T>(patches, bconditions, bases)
    {
        this->m_unknownDim.resize(1);
        this->m_unknownDim[0] = rhs.targetDim();
        m_dirStrategy = dirichlet::elimination;
        m_isSymmetric = true;
        m_interfaceStrategy = iFace::glue;
        m_rhs_function = rhs.clone();
    }

    // gsPoissonAssembler( gsBVProblem<T> const & bvp, )
    // 


   ~gsPoissonAssembler()
    {
        delete m_rhs_function;
    }


    /// Set the strategy for dealing with Dirichlet dofs.\n
    /// \param strategy options are:\n
    /// Eliminate the Dirichlet dofs from the linear system.\n
    /// \em elimination = 1 (default)\n\n
    /// Keep the Dirichlet dofs and enforce the boundary
    /// condition weakly by a penalty term.\n
    /// \em Nitsche = 2
    void setDirichletStrategy(dirichlet::strategy strategy)
    { m_dirStrategy = strategy; }

    /// Set the strategy for dealing with patch interface.\n
    /// \param strategy options are:\n
    /// Glue patches together by merging dofs across an interface into one.
    /// This only works for conforming interfaces.\n
    /// \em glue = 1 (default)\n\n
    /// Use discontinuous Galerkin-like coupling between adjacent patches.\n
    /// \em dg = 2
    void setInterfaceStrategy(iFace::strategy strategy)
    { m_interfaceStrategy = strategy; }


    //Virtual function
    public:

    /// Initializes things. To be overwritten in derived classes.
    /// This will be called in solve()
    ///
    /// Typically:
    /// Computes any fixed DoFs, initializes some internal members,
    /// checks consistency
    void initialize();


    /// @brief Assembles the final system with all boundary conditions contained
    void assemble();



    // Solves the linear system and fills up \a m_sysSolution
    /**
     * @brief solveSystem solves the linear system of equations Ax = b
     * with Eigen's ConjugateGradient solver and fills up \a m_sysSolution.
     * Displays the residual error and number of iterations.
     */
    void solveSystem();

    /// @brief Solve the Poisson, e.i initialize, assemble, solve and reconstruct solution.
    void solve()
    {
        initialize ();
        assemble();
        solveSystem();
        this->reconstructSolution(true);
    }


    /// \brief Applies an error indicator to the computed discrete solution.
    ///
    /// Applies a posteriori error estimation and returns the resulting
    /// error indicator
    /// as a gsVector \em ErrInd of length \em N, where \em N is the
    /// number of patches (a.k.a. subdomains) of the considered gsMultiPatch geometry.
    ///
    /// Let <em>ErrInd_k = ErrInd[k]</em> be the entry corresponding to
    /// the patch with index \em k.
    /// <em>ErrInd_k</em> is a gsMatrix of size <em>NE</em> x <em>(1+2*d)</em>, where\n
    /// <em>NE</em> is the number of elements and\n
    /// <em>d</em> is the dimension of the parameter domain.\n
    /// \n
    /// The entries of the gsMatrix <em>ErrInd_k</em> are as follows:\n
    /// <em>ErrInd_k(i,0)</em> = estimated error on element \em i.\n
    /// <em>ErrInd_k(i,1,...,d)</em> = coordinates of lower corner of the element.\n
    /// <em>ErrInd_k(i,d+1,...,2*d)</em> = coordinates of upper corner of the element.\n
    ///
    /// This should become the main function in which the desired
    /// error indicator is selected.\n
    ///
    /// Currently (23.Apr.2014), only errIndPatch_Bubble() is implemented,
    /// so the choice is rather simple.
    ///
    /// \param[out] ErrInd gsVector of length \em N, where \em N is the
    /// number of patches of the considered gsMultiPatch geometry. Each entry of
    /// the gsVector is a gsMatrix containing the local error indicator information.
    /// See above for the format.
    ///
    /// \todo Add more error indicators.
    /// \todo Include a nice way of selecting one of these error indicators.
    gsVector< gsMatrix<T> > errorIndicator();

    /// \brief Error indicator based on bubble functions.
    ///
    /// See documentation for errorIndicator() for the details on the
    /// output data.\n
    /// \n
    /// This error indicator is based on the idea presented in\n
    /// <em>M. R. Doerfel, B. Juettler, B. Simeon. Adaptive isogeometric analysis by local
    /// h-refinement with T-splines. Computer Methods in Applied Mechanics and
    /// Engineering, 199: 264-275, 2010.</em>\n
    /// The bubble function which are used here are realized with B-splines on the
    /// unit interval with degree <em>p+1</em>, where \em p
    /// is the degree of the discrete solution.\n
    ///
    /// For details, see the mentioned reference and, as usual,
    /// the references therein.\n
    ///
    /// \param[in] patchIndex Index of the considered patch.
    /// \param[out] ErrInd gsMatrix of size <em>NE</em> x <em>(1+2*d)</em>, where\n
    /// <em>NE</em> is the number of elements and\n
    /// <em>d</em> is the dimension of the parameter domain.\n
    /// See documentation of errorIndicator() for details.
    gsMatrix<T> errIndPatch_Bubble( const index_t patchIndex = 0);


    /// \brief Apply adaptive refinement based on the errors given in \em errInd.
    ///
    /// Refinement is based on the local error estimates given in \em errInd
    /// (see documentation of errorIndicator() for details on the format).\n
    /// The parameter \em a is used for computing the threshold \em Thr. Each element
    /// with a local error estimate larger than \em Thr will be refined.\n
    /// <em>a = 0</em> corresponds to \em global refinement.\n
    /// <em>a = 1</em> corresponds to no refinement/the refinement of only one element.\n
    /// \n
    /// Three methods for computing \em Thr are implemented (as of 29.Apr.2014):\n
    /// <b>Absolute threshold</b> (<em>refCriterion=1</em>):\n
    /// <em>Thr = a * (maximum of all local error estimates)</em>.\n
    /// <b>Relative threshold</b> (<em>refCriterion=2</em>):\n
    /// \em Thr is chosen such that <em>(1-a)*100</em> percent of all cells are refined.\n
    /// <b>"Bulk-criterion"</b> (<em>refCriterion=3</em>):\n
    /// A certain number of elements with the largest error estimates are chosen for refinement, such
    /// that the sum of their local error estimates add up to <em>a * (sum of all local
    /// error estimates)</em>.
    ///
    /// \param[in] errInd See documentation for errorIndicator() for details on the format.
    /// \param[in] refCriterion Flag specifying method for computing the threshold \em Thr.
    /// \param[in] a Parameter for computing the threshold \em Thr.
    /// \param NOTE The member m_bases of gsPdeAssembler is directly modified within this function.
    ///
    ///
    void adaptiveRefine( const gsVector< gsMatrix<T> >& errInd, const int refCriterion = 1,  const T a = 0.7);


    void uniformRefine();

//Accessors
public:

    /// @brief Return the right-hand side function \f$\mathbf{f}\f$
    const gsFunction<T> & rhs_function() const  { return *m_rhs_function; }


    /// \brief Computes and returns the full stiffness matrix on patch \em int.
    gsSparseMatrix<T> * stiffnessMatrixPatch(int patchIndex = 0);

    /// \brief Computes and returns the full mass matrix on patch \em int.
    gsSparseMatrix<T> * massMatrixPatch(int patchIndex = 0);


    /**
     * \brief Fix certain DOF to certain values.
     *
     * This function fixes certain degrees of freedom (DOF) to have certain
     * values and corrects the stiffness matrix and the right-hand-side
     * of the problem.
     *
     * \param[in] Idx gsVector of length \em N, containing
     * the indices of the DOFs that will be fixed.
     * \param[in] Val gsMatrix containing the corresponding values.\n
     * The size of \em Val is \em N x \em k, where \em k is the number of
     * unknowns in the vector-valued PDE. The <em>k</em>-th component of
     * the DOF with global index \em Idx[i] will be set to the value <em>Val(i,k)</em>.
     *
     * Note that this function cannot yet be used inside the solve() function,
     * because it needs values as input!
     *
     * <b>Suggested use:</b>\n
     * 1. Do NOT specify any Dirichlet boundary conditions.
     * 2. Compute which DOF have to be fixed to which value with your own routine,
     * save them in, e.g., \em "Indices" and \em "Values", respectively.
     * 3. Instead of gsPoissonAssembler::solve(), call\n
     * gsPoissonAssembler::initialize();\n
     * gsPoissonAssembler::assemble();\n
     * gsPoissonAssembler::boundaryFixDofs( Indices, Values );\n
     * gsPoissonAssembler::solveSystem();\n
     * gsPoissonAssembler::reconstructSolution();\n
     *
     * \warning Note that, even though this function is intended for
     * fixing DOFs at the Dirichlet boundary, it is technically possible to fix
     * ANY degree of freedom to a certain value. It is NOT checked whether the
     * indices in \em Idx really correspond to DOFs at the domain boundary!
    */
    void boundaryFixDofs( const gsVector<unsigned> & Idx ,
                            const gsMatrix<T> & Val );

    /// Add contribution of boundary condition \a bc to the linear system.
    /// \param B            discretization basis for the patch
    /// \param bc           the boundary condition to apply
    /// \param mapper       contains the global Dof map for the linear system
    void applyBoundary( const gsBasis<T>   & B,
                        const boundary_condition<T> & bc,
                        const gsDofMapper& mapper);

    /// Add contribution of interface \a bi to the system \a system_matrix
    /// \param B1       discretization basis for the first patch
    /// \param B2       discretization basis for the second patch
    /// \param geo1     the first patch
    /// \param geo2     the second patch
    /// \param bi       description of the interface which joins the two patches
    /// \param mapper   contains the global Dof map for the \a system
    /// \param system_matrix   the system matrix where the contribution shall be added
    void applyDG( const gsBasis<T> & B1, const gsBasis<T> & B2,
                  const gsGeometry<T> & geo1,
                  const gsGeometry<T> & geo2,
                  const boundaryInterface & bi,
                  const gsDofMapper& mapper,
                  gsSparseMatrix<T> & system_matrix );

private:

    /// Assembles on a single patch \a patchIndex
    void assemblePatch(int patchIndex=0);
    /// Assembles on every patch by calling \a assemblePatch for each patch.
    void assembleMultipatch();


    /// Prints the object as a string.
    std::ostream &print(std::ostream &os) const
    {   os<<"Poisson's equation  -\u0394u = f ,  with:\n";
        os << "* f: "<< *m_rhs_function <<"\n";
        os << "* Domain: "<< m_patches ;
        os << *m_bconditions[0];
        return os;
    }


private:
    // Right hand side function
    gsFunction<T> * m_rhs_function;

    // Strategy for dealing with Dirichlet dofs
    dirichlet::strategy m_dirStrategy;

    // Strategy for dealing with patch interface
    iFace::strategy m_interfaceStrategy;

    //bool m_eliminateDirichlet

    // Members from gsPdeAssembler
    using gsPdeAssembler<T>::m_patches;
    using gsPdeAssembler<T>::m_bconditions;
    using gsPdeAssembler<T>::m_bases;

    using gsPdeAssembler<T>::m_matrix;
    using gsPdeAssembler<T>::m_rhs;
    using gsPdeAssembler<T>::m_sysSolution;
    using gsPdeAssembler<T>::m_fixedDofs;
    using gsPdeAssembler<T>::m_dofMapper;

    using gsPdeAssembler<T>::m_isSymmetric;

    using gsPdeAssembler<T>::m_dofs;
    using gsPdeAssembler<T>::m_idofs;

    using gsPdeAssembler<T>::m_unknownDim;
    using gsPdeAssembler<T>::m_solutions;


}; // class gsPoissonAssembler


} // namespace gismo

#include "gsPoissonAssembler.hpp"
//#ifndef GISMO_HEADERS_ONLY
//#include GISMO_HPP_HEADER(gsPoissonAssembler.hpp)
//#endif
