#pragma once

#include <gsAssembler/gsPdeAssembler.h>
#include <gsAssembler/gsSolverOptions.h>
#include <gsPde/gsBoundaryConditions.h>

namespace gismo
{

/** @brief Implementation of an scalar convection-diffusion-reaction equation solver.

    \warning Not fully tested yet (01.Jul.2014)

    The considered equation:
    \f[ - \nabla \cdot ( A \nabla u ) + b \cdot \nabla u + c\, u = f \f]
    where \f$ u: \mathbb{R}^d \to \mathbb R\f$,\n
    \f$A\f$ is a \f$ d\times d\f$-matrix,\n
    \f$b\f$ is a \f$ d\times 1\f$-vector, and\n
    \f$c\f$ is scalar.\n

    See gsPdeAssembler for general description of input.
    \param[in] patches is a gsMultiPatch object describing the geometry.
    \param[in] bconditions is a gsBoundaryConditions object.
    \param[in] bases is a vector of gsBasis for each components of the unknown
    of \f${u}\f$.
    \param[in] rhs is the right-hand side of the Poisson equation, \f${f}\f$.
    \param[in] Coeff_A is a gsFunction with 4 (d=2) or 9 (d=3) components
    (see below for details).
    \param[in] Coeff_b is a gsFunction with 2 (d=2) or 3 (d=3) components.
    \param[in] Coeff_c is a gsFunction with 1.

    <b>Coeff_A</b> has to be a gsFunction, which returns a column <em>C</em>
    for each
    evaluation point, and this column should represent \n
    <em>C</em>=\f$ ( A_{00},\ A_{10},\ A_{01},\ A_{11})^T\f$ for d=2, and\n
    <em>C</em>=\f$ ( A_{00},\ A_{10},\ A_{20},\ A_{01},\ \ldots,\ A_{22})^T\f$ for d=3,\n
    such that <em>C.resize(d,d)</em> gives back the matrix \f$A\f$ as in the PDE.

    It sets up an assembler and assembles the system patch wise and combines
    the patch-local stiffness matrices into a global system by various methods
    (see gismo::gsInterfaceStrategy). It can also enforce Dirichlet boundary
    conditions in various ways (see gismo::gsDirichletStrategy).

    SUPG-stabilization can be used for stabilization. It can be turned off by
    setting the member m_doSUPG to \m false.\n
    \warning As of now (01.Jul.2014), in the SUPG-stabilization the following
    terms are being neglected:\n
    The derivative of the convection vector \f$ b \f$, and\n
    the second derivative of the inverse geometry mapping.

    By calling gsConvDiffReactAssembler::solve(), the linear system is assembled and solved, and a
    gsField containing the computed solution is returned.

    <b>Example:</b>\n
    If you want to solve the convection-diffusion equation
    \f[ -\kappa \Delta u + b \cdot \nabla u = 0\f]
    in \f$\mathbb R^2\f$ with \f$ \kappa = 10^{-6},\ b = ( \cos(\pi/4), \sin(\pi/4))^T\f$, then
    this can be realized with the following coefficients
    \verbatim
  gsMFunctionExpr<T>  Coeff_A("10^(-6)","0","0","10^(-6)", 2);
  gsMFunctionExpr<T>  Coeff_b("1/sqrt(2)","1/sqrt(2)", 2);
  gsMFunctionExpr<T>  Coeff_c("0", 2);
  gsMFunctionExpr<T>  f("0", 2);
    \endverbatim





*/

template<class T>
class gsConvDiffReactAssembler : public gsPdeAssembler<T>
{
public:
    /// Default empty constructor
    gsConvDiffReactAssembler () : gsPdeAssembler<T>() { }

    ///Constructor
    gsConvDiffReactAssembler( gsMultiPatch<T> const & patches,
                     gsBoundaryConditions<T> const & bconditions,
                     std::vector< gsBasis<T>* > const & bases,
                     const gsFunction<T> & rhs,
                     const gsFunction<T> & coeff_A,
                     const gsFunction<T> & coeff_b,
                     const gsFunction<T> & coeff_c )
         : gsPdeAssembler<T>(patches, bconditions, bases)
    {
        this->m_unknownDim.resize(1);
        this->m_unknownDim[0] = rhs.targetDim();
        m_dirStrategy = elimination;
        m_isSymmetric = false;
        m_interfaceStrategy = glue;
        m_rhs_function = rhs.clone();
        m_coeff_A = coeff_A.clone();
        m_coeff_b = coeff_b.clone();
        m_coeff_c = coeff_c.clone();
        m_doSUPG = true;
    }


    ///Constructor.
    gsConvDiffReactAssembler( gsMultiPatch<T> const & patches,
                     std::vector<gsBoundaryConditions<T>*> const & bconditions,
                     std::vector< gsBasis<T>* > const & bases,
                     const gsFunction<T> & rhs,
                     const gsFunction<T> & coeff_A,
                     const gsFunction<T> & coeff_b,
                     const gsFunction<T> & coeff_c )
        : gsPdeAssembler<T>(patches, bconditions, bases)
    {
        this->m_unknownDim.resize(1);
        this->m_unknownDim[0] = rhs.targetDim();
        m_dirStrategy = elimination;
        m_isSymmetric = false;
        m_interfaceStrategy = glue;
        m_rhs_function = rhs.clone();
        m_coeff_A = coeff_A.clone();
        m_coeff_b = coeff_b.clone();
        m_coeff_c = coeff_c.clone();
        m_doSUPG = true;
    }

    ///Constructor with gsMultiBasis
    gsConvDiffReactAssembler( gsMultiPatch<T> const & patches,
                     gsBoundaryConditions<T> const & bconditions,
                     gsMultiBasis<T> const & bases,
                     const gsFunction<T> & rhs,
                     const gsFunction<T> & coeff_A,
                     const gsFunction<T> & coeff_b,
                     const gsFunction<T> & coeff_c )
        : gsPdeAssembler<T>(patches, bconditions, bases)
    {
        this->m_unknownDim.resize(1);
        this->m_unknownDim[0] = rhs.targetDim();
        m_dirStrategy = elimination;
        m_isSymmetric = true;
        m_interfaceStrategy = glue;
        m_rhs_function = rhs.clone();
        m_coeff_A = coeff_A.clone();
        m_coeff_b = coeff_b.clone();
        m_coeff_c = coeff_c.clone();
        m_doSUPG = true;
    }


    // gsConvDiffReactAssembler( gsBVProblem<T> const & bvp, )
    //


   ~gsConvDiffReactAssembler()
    {
        delete m_rhs_function;
        delete m_coeff_A;
        delete m_coeff_b;
        delete m_coeff_c;
    }


    /// Set the strategy for dealing with Dirichlet dofs.\n
    /// \param strategy options are:\n
    /// Eliminate the Dirichlet dofs from the linear system.\n
    /// \em elimination = 1 (default)\n\n
    /// Keep the Dirichlet dofs and enforce the boundary
    /// condition weakly by a penalty term.\n
    /// \em Nitsche = 2
    void setDirichletStrategy(gsDirichletStrategy strategy)
    { m_dirStrategy = strategy; }

    /// Set the strategy for dealing with patch interface.\n
    /// \param strategy options are:\n
    /// Glue patches together by merging dofs across an interface into one.
    /// This only works for conforming interfaces.\n
    /// \em glue = 1 (default)\n\n
    /// Use discontinuous Galerkin-like coupling between adjacent patches.\n
    /// \em dg = 2
    void setInterfaceStrategy(gsInterfaceStrategy strategy)
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


    //Error estimation and adaptive refinement
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
    /// Currently (06.Oct.2014), only errIndPatch_Bubble()
    /// and errIndPatch_Residual() are implemented.
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


    // S.Kleiss
    /** \brief Residual-based error indicator.
    *
    * Computes the simple error indicator \f$\eta_K\f$ for an element/cell \f$K\f$
    * via the residual
    \f[
     \eta_K^2 = \| f + \nabla \cdot ( A \nabla u_h ) - b \cdot \nabla u_h - c u_h \|_{L_2}^2.
    \f]
    *
    * \warning Jumps of \f$ \nabla u_h \f$ across element interfaces are \b neglected
    * (i.e., \f$C^1\f$-continuity of \f$u_h\f$ is assumed)!!!\n
    * The term corresponding to the Neumann boundary is also \b neglected !!!
    *
    * See documentation for errorIndicator() for the details on the
    * output data.\n
    *
    * \param[in] patchIndex Index of the considered patch.
    * \param[out] ErrInd gsMatrix of size <em>NE</em> x <em>(1+2*d)</em>, where\n
    * <em>NE</em> is the number of elements and\n
    * <em>d</em> is the dimension of the parameter domain.\n
    * See documentation of errorIndicator() for details.
    */
    gsMatrix<T> errIndPatch_Residual( const index_t patchIndex = 0);

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
    /// \param[in] refCriterion Flag specifying method for computing the threshold \em Thr (see above).
    /// \param[in] a Parameter for computing the threshold \em Thr.
    /// \param[in] RefExtension Extends the marked areas by a <em>RefExtension</em>-ring of cells. Note that the ring consists of cells of the finer level. <b>This refinement-extension only works with gsHTensorBasis and derived basis classes</b>. When <em>RefExtension = 0</em>, only the marked areas are refined.
    ///
    /// \remarks The member m_bases of gsPdeAssembler is directly <b>modified</b> within this function.
    ///
    ///
    void adaptiveRefine( const gsVector< gsMatrix<T> >& errInd, const int refCriterion = 1,  const T a = 0.7, const unsigned RefExtension = 0);


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
     * 3. Instead of gsConvDiffReactAssembler::solve(), call\n
     * gsConvDiffReactAssembler::initialize();\n
     * gsConvDiffReactAssembler::assemble();\n
     * gsConvDiffReactAssembler::boundaryFixDofs( Indices, Values );\n
     * gsConvDiffReactAssembler::solveSystem();\n
     * gsConvDiffReactAssembler::reconstructSolution();\n
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

    /** \brief Computes the SUPG-stabilization-parameter on a cell specified by
     * \em lo and \em up
     *
     * The parameter \f$\tau_K\f$ on a element \f$K\f$ is computed as
     * \f[ \tau_K = \frac{ h_b(K) }{2\, |b|}, \f]
     * where \f$h_b(K)\f$ is the length of the cell \f$K\f$ in direction of the
     * convection \f$b\f$, and \f$ |b| \f$ is the magnitude of the convecton
     * velocity.
     *
     * \param[in] lo gsVector of size \f$d\f$, containing the coordinates of the
     * upper corner of the element, as returned by gsDomainIterator::lowerCorner().
     * \param[in] up gsVector of size \f$d\f$, containing the coordinates of the
     * upper corner of the element, as returned by gsDomainIterator::upperCorner().
     *
     * The size \f$h_b(K)\f$ is estimated by computing some points on the boundary
     * of the element, mapping these points to the phyiscal domain, projecting these
     * points onto the vector \f$b\f$, and then taking the distance between the
     * largest and the smallest of these projected points.\n
     * Maybe not soooo efficient, but working.
     */
    T get_SUPG_parameter( const gsVector<T> & lo,
                          const gsVector<T> & up,
                          const int patchIdx = 0 );

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

public:
    bool m_doSUPG;

private:
    // Right hand side function
    gsFunction<T> * m_rhs_function;
    // Strategy for dealing with Dirichlet dofs
    gsDirichletStrategy m_dirStrategy;
    // Strategy for dealing with patch interface
    gsInterfaceStrategy m_interfaceStrategy;

    // gsFunctions for the three coefficients
    // defining the convecton-diffusion-reaction-problem
    gsFunction<T> * m_coeff_A;
    gsFunction<T> * m_coeff_b;
    gsFunction<T> * m_coeff_c;

    // Members from gsPdeAssembler
    using gsPdeAssembler<T>::m_patches;
    using gsPdeAssembler<T>::m_bconditions;
    using gsPdeAssembler<T>::m_bases;

public:
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


}; // class gsConvDiffReactAssembler


} // namespace gismo

#include "gsConvDiffReactAssembler.hpp"
//#ifndef GISMO_HEADERS_ONLY
//#include GISMO_HPP_HEADER(gsConvDiffReactAssembler.hpp)
//#endif
