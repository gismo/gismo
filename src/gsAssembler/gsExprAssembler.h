/** @file gsExprAssembler.h

    @brief Generic expressions matrix assembly

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris
*/

#pragma once

#ifdef _OPENMP
#include <omp.h>
#endif

#include <gsUtils/gsPointGrid.h>
#include <gsAssembler/gsQuadrature.h>
#include <gsExpressions/gsExprHelper.h>
#include <gsExpressions/gsFeSpaceData.h>
#include <gsDomain/gsDomain.h>

#include <gsAssembler/gsCPPInterface.h>

#include <gsMatrix/gsFiberMatrix.h>

namespace gismo
{

/**
   \brief Pattern-lock pool size for striped locking.
   Fixed power-of-two cap (1<<16) so one omp_lock_t per global DOF is never
   allocated. Threads acquire at most one lock at a time before any insertion,
   so striping cannot deadlock. Behaviour-identical to per-DOF locking.
*/
static inline index_t patternLockPoolSize(index_t numDofs)
{
    index_t cap = (index_t(1) << 16);
    index_t s = 1;
    while (s < numDofs && s < cap) s <<= 1;
    return s;
}

/**
   Assembler class for generating matrices and right-hand sides based
   on isogeometric expressions
*/
template<class T>
class gsExprAssembler
{
private:
    typedef typename gsDomain<T>::iterator elementIterator;

    typename gsExprHelper<T>::Ptr m_exprdata;
    const gsMultiPatch<T>* m_gmap;

    gsOptionList m_options;

    mutable gsSparseMatrix<T> m_matrix;
    typedef gsFiberMatrix<T,ColMajor> FiberMatrix;
    FiberMatrix m_fmatrix;
    gsMatrix<T>      m_rhs;

    std::list<gismo::expr::gsFeSpaceData<T> > m_sdata;
    std::vector<gismo::expr::gsFeSpaceData<T>*> m_vrow;
    std::vector<gismo::expr::gsFeSpaceData<T>*> m_vcol;

    int m_sparsity;//0:unknown, 1:volume, 2:boundary, 4:interface pre-allocated
    mutable bool m_modified;

    typedef typename gsExprHelper<T>::nullExpr    nullExpr;

public:

    typedef typename gsSparseMatrix<T>::BlockView matBlockView;

    typedef typename gsSparseMatrix<T>::constBlockView matConstBlockView;

    typedef typename gsBoundaryConditions<T>::bcRefList   bcRefList;
    typedef gsBoxTopology::bContainer  bContainer;
    typedef gsBoxTopology::ifContainer ifContainer;

    typedef typename gsExprHelper<T>::element     element;     ///< Current element
    typedef typename gsExprHelper<T>::geometryMap geometryMap; ///< Geometry map type
    typedef typename gsExprHelper<T>::variable    variable;    ///< Variable type
    typedef typename gsExprHelper<T>::space       space;       ///< Space type
    typedef typename expr::gsFeSolution<T>        solution;    ///< Solution type

    typedef typename gsQuadRule<T>::uPtr QuadratureRulePtr;

    /**
     * @brief Factory for an opt-in custom quadrature rule.
     *
     * The factory is called once for every rule instance required by an
     * assembly operation (in particular, once per OpenMP worker and patch for
     * volume assembly). It must return a fresh rule and, when OpenMP is enabled,
     * be safe to call concurrently. @a fixedDirection is -1 for volume
     * integration and the fixed parametric direction for boundary and interface
     * integration. Overriding gsQuadRule::mapTo() permits element-dependent
     * rules.
     */
    typedef std::function<QuadratureRulePtr(const gsBasis<T> & basis,
                                             const gsOptionList & options,
                                             index_t patch,
                                             short_t fixedDirection)>
        QuadratureFactory;

private:
    QuadratureFactory m_quadratureFactory;

public:

    void cleanUp()
    {
        m_exprdata->cleanUp();
    }

    /// Constructor
    /// \param _rBlocks Number of spaces for test functions
    /// \param _cBlocks Number of spaces for solution variables
    gsExprAssembler(index_t _rBlocks = 1, index_t _cBlocks = 1)
    : m_exprdata(gsExprHelper<T>::make()), m_gmap(nullptr), m_options(defaultOptions()),
      m_vrow(_rBlocks,nullptr), m_vcol(_cBlocks,nullptr), m_sparsity(0), m_modified(false)
    { }

    // The copy constructor replicates the same environent but does
    // not copy the expression helper

    /// @brief Returns the list of default options for assembly
    static gsOptionList defaultOptions();

    /// Returns the number of degrees of freedom (after initialization)
    index_t numDofs() const
    {
        GISMO_ASSERT( m_vcol.back()->mapper.isFinalized(),
                      "gsExprAssembler::numDofs() says: initSystem() has not been called.");
        return m_vcol.back()->mapper.firstIndex() +
	  m_vcol.back()->mapper.freeSize();
    }

    /// Returns the number of test functions (after initialization)
    index_t numTestDofs() const
    {
        GISMO_ASSERT( m_vrow.back()->mapper.isFinalized(),
                      "initSystem() has not been called.");
        return m_vrow.back()->mapper.firstIndex() +
	  m_vrow.back()->mapper.freeSize();
    }

    /// Returns the number of blocks in the matrix, corresponding to
    /// variables/components
    index_t numBlocks() const
    {
        index_t nb = 0;
        for (size_t i = 0; i!=m_vrow.size(); ++i)
            nb += m_vrow[i]->dim;
        return nb;
    }

    /// Returns a reference to the options structure
    gsOptionList & options() {return m_options;}
    const gsOptionList & options() const {return m_options;}


    /// @brief Installs a custom quadrature-rule factory.
    ///
    /// Passing an empty factory restores the standard option-driven
    /// quadrature. The factory is invoked only when a new rule is needed, never
    /// in an element or quadrature-point loop.
    void setQuadratureFactory(QuadratureFactory factory)
    { m_quadratureFactory = give(factory); }

    /// @brief Restores the standard gsQuadrature/options-based rules.
    void clearQuadratureFactory()
    { m_quadratureFactory = QuadratureFactory(); }

    /// @brief Returns whether a custom quadrature factory is installed.
    bool hasCustomQuadrature() const
    { return static_cast<bool>(m_quadratureFactory); }

    /// Returns the internally stored sparse fiber matrix
    const FiberMatrix & fiberMatrix() const
    { return m_fmatrix; }

    /// @brief Returns the left-hand global matrix
    const gsSparseMatrix<T> & matrix() const
    { return m_modified ? makeMatrix() : m_matrix; }

    /// When calling assemble, the matrix is not filled (but only stored internally).
    /// Call this function to fill the sparsematrix with all the assemblies so far
    const gsSparseMatrix<T> & makeMatrix() const
    {
        m_fmatrix.toSparseMatrix_into(m_matrix);
        m_modified = false;
        return m_matrix;
    }

    /// @brief Writes the resulting matrix in \a out. The internal matrix is moved.
    void matrix_into(gsSparseMatrix<T> & out)
    {
        matrix();
        out = give(m_matrix);
    }

    EIGEN_STRONG_INLINE gsSparseMatrix<T> giveMatrix()
    {
        matrix();
        m_modified = true;
        return give(m_matrix);
    }

    EIGEN_STRONG_INLINE FiberMatrix giveFiberMatrix()
    {
        return give(m_fmatrix);
    }

    /// @brief Returns the right-hand side vector(s)
    const gsMatrix<T> & rhs() const { return m_rhs; }

    /// @brief Writes the resulting vector in \a out. The internal data is moved.
    void rhs_into(gsMatrix<T> & out) { out = give(m_rhs); }

    /// \brief Sets the domain of integration.
    /// \warning Must be called before any computation is requested
    void setIntegrationElements(const gsMultiBasis<T> & mesh)
    { m_exprdata->setDomain(mesh.domain()); m_sparsity = 0; }

    /// \brief Sets the domain of integration.
    /// \warning Must be called before any computation is requested
    void setIntegrationDomain(typename gsDomain<T>::Ptr domain)
    { m_exprdata->setDomain(give(domain)); m_sparsity = 0; }


    /// \brief Set the geometrymap ( used for interface assembly)
    /// \warning Must be called before any computation is requested
    void setGeometryMap(const gsMultiPatch<T> & gMap)
    { m_gmap = &gMap;}

    const gsMultiPatch<T>& getGeometryMap() const
    {
        return (nullptr == m_gmap ? m_exprdata->multiPatch() : *m_gmap);
    }

#if EIGEN_HAS_RVALUE_REFERENCES
    void setIntegrationElements(const gsMultiBasis<T> &&) = delete;
    //const gsMultiBasis<T> * c++98
#endif

    /// \brief Returns the domain of integration
    const gsDomain<T> & domain() const
    { return m_exprdata->domain(); }

    const typename gsExprHelper<T>::Ptr exprData() const { return m_exprdata; }

    /// Registers \a g as an isogeometric geometry map and return a handle to it
    geometryMap getMap(const gsFunctionSet<T> & g)
    { return m_exprdata->getMap(g); }

    /// Registers \a mp as an isogeometric (both trial and test) space
    /// and returns a handle to it
    space getSpace(const gsFunctionSet<T> & mp, index_t dim = 1, index_t id = 0)
    {
        //if multiBasisSet() then check domainDom
        GISMO_ASSERT(1==mp.targetDim(), "Expecting scalar source space");
        GISMO_ASSERT(static_cast<size_t>(id)<m_vcol.size(),
                     "Given ID "<<id<<" exceeds "<<m_vcol.size()-1 );

        if (m_vcol[id]==nullptr)
        {
            m_sdata.emplace_back(mp,dim,id);
            m_vcol[id] = &m_sdata.back();
            if ((size_t)id<m_vrow.size() && nullptr==m_vrow[id]) m_vrow[id]=m_vcol[id];
        }
        else
        {
            m_vcol[id]->fs  = &mp;
            m_vcol[id]->dim = dim;
        }

        expr::gsFeSpace<T> u = m_exprdata->getSpace(mp,dim);
        u.setSpaceData(*m_vcol[id]);
        return u;
    }

    /// Registers \a mp as an isogeometric test (row) space and returns
    /// a handle to it
    space getTestSpace(const gsFunctionSet<T> & mp, index_t dim = 1, index_t id = 0)
    {
        GISMO_ASSERT(1==mp.targetDim(), "Expecting scalar source space");
        GISMO_ASSERT(static_cast<size_t>(id)<m_vrow.size(),
                     "Given ID "<<id<<" exceeds "<<m_vrow.size()-1 );

        if ( (m_vrow[id]==nullptr) ||
             ((size_t)id<m_vcol.size() && m_vrow[id]==m_vcol[id]) )
        {
            m_sdata.emplace_back(mp,dim,id);
            m_vrow[id] = &m_sdata.back();
        }
        else
        {
            m_vrow[id]->fs  = &mp;
            m_vrow[id]->dim = dim;
        }

        expr::gsFeSpace<T> s = m_exprdata->getSpace(mp,dim);
        s.setSpaceData(*m_vrow[id]);
        return s;
    }

    /// \brief Registers \a mp as an isogeometric test space
    /// corresponding to trial space \a u and return a handle to it
    ///
    /// \note Both test and trial spaces are registered at once by
    /// gsExprAssembler::getSpace.
    ///
    /// Use this function after calling gsExprAssembler::getSpace when
    /// a distinct test space is requred (eg. Petrov-Galerkin
    /// methods).
    ///
    /// \note The dimension is set to the same as \a u, unless the caller
    /// sets as a third argument a new value.
    space getTestSpace(space u, const gsFunctionSet<T> & mp, index_t dim = -1)
    { return getTestSpace( mp,(-1 == dim ? u.dim() : dim), u.id() ); }

    /// Return a variable handle (previously created by getSpace) for
    /// unknown \a id
    space trialSpace(const index_t id) const
    {
        GISMO_ASSERT(NULL!=m_vcol[id], "Not set.");
        expr::gsFeSpace<T> s = m_exprdata->
            getSpace(*m_vcol[id]->fs,m_vcol[id]->dim);
        s.setSpaceData(*m_vcol[id]);
        return s;
    }

    /// Return the trial space of a pre-existing test space \a v
    space trialSpace(space & v) const { return trialSpace(v.id()); }

    /// Return the variable (previously created by getTrialSpace) with
    /// the given \a id
    space testSpace(const index_t id)
    {
        GISMO_ASSERT(NULL!=m_vrow[id], "Not set.");
        expr::gsFeSpace<T> s = m_exprdata->
            getSpace(*m_vrow[id]->fs,m_vrow[id]->dim);
        s.setSpaceData(*m_vrow[id]);
        return s;
    }

    /// Return the test space of a pre-existing trial space \a u
    space testSpace(space u) const { return testSpace(u.id()); }

    /// Registers \a func as a variable and returns a handle to it
    ///
    variable getCoeff(const gsFunctionSet<T> & func)
    { return m_exprdata->getVar(func, 1); }

    /// Registers \a func as a variable defined on \a G and returns a
    /// handle to it
    expr::gsComposition<T> getCoeff(const gsFunctionSet<T> & func, geometryMap & G)
    { return m_exprdata->getVar(func,G); }

    /// \brief Registers a representation of a solution variable from
    /// space \a s, based on the vector \a cf.
    ///
    /// The vector \a cf should have the structure of the columns of
    /// the system matrix this->matrix(). The returned handle
    /// corresponds to a function in the space \a s
    solution getSolution(const expr::gsFeSpace<T> & s, gsMatrix<T> & cf) const
    { return solution(s, cf); }

    variable getBdrFunction() const { return m_exprdata->getMutVar(); }

    expr::gsComposition<T> getBdrFunction(geometryMap & G) const
    { return m_exprdata->getMutVar(G); }

    element getElement() const { return m_exprdata->getElement(); }

    // note: not used
    void setFixedDofVector(gsMatrix<T> & dof, short_t unk = 0);
    // note: not used
    void setFixedDofs(const gsMatrix<T> & coefMatrix, short_t unk = 0, size_t patch = 0);

    /// \brief Initializes the sparse system (sparse matrix and rhs)
    void initSystem(const index_t numRhs = 1)
    {
        // Check spaces.nPatches==mesh.patches
        initMatrix();
        m_rhs.setZero(numTestDofs(), numRhs);
    }

    /// \brief Initializes the sparse matrix only
    void initMatrix()
    {
        resetDimensions();
        clearMatrix(false);
    }

    // @hverhelst: adds explicit size, because if the RHS is moved ('given'), its sizes are lost.
    void clearRhs(const index_t numRhs = 1) { m_rhs.setZero(numTestDofs(),numRhs); }

    /**
     * @brief Re-Init Matrix (set zero by default)
     *
     * @param save_sparsety_pattern only modify values but keep sparsety
     * information by multiplying matrix by zero in-place
     */
    void clearMatrix(const bool& save_sparsety_pattern = true)
    {
        if (m_fmatrix.nonZeros() && save_sparsety_pattern)
        {
            m_fmatrix.assignZero();
        }
        else
        {
            if (m_options.askSwitch("lazyMatrix",false))
                m_fmatrix.resizeLazy(numTestDofs(), numDofs());
            else
                m_fmatrix.resize(numTestDofs(), numDofs());
            m_sparsity = 0;

            if (0 == m_fmatrix.rows() || 0 == m_fmatrix.cols())
                gsWarn << " No internal DOFs, zero sized system.\n";
            else {
                // Pick up values from options
                const T bdA = m_options.getReal("bdA");
                const index_t bdB = m_options.getInt("bdB");
                const T bdO = m_options.getReal("bdO");
                T nz = 1;
                const short_t dim = m_exprdata->domain().dim();
                for (short_t i = 0; i != dim; ++i)
                    nz *= bdA * static_cast<T>(
                                    m_exprdata->domain().degree(i)) +
                          static_cast<T>(bdB);

                m_fmatrix.reservePerColumn(numBlocks() *
                                          cast<T, index_t>(nz * (1.0 + bdO)));
            }
        }
		
        m_modified = true;
    }

    /// Initializes the pattern of the sparse matrix. The m_sparsity bit is
    /// set inside _computePattern() itself (Fix 2), only when a pattern was
    /// actually (re)computed -- i.e. only when one of \a args is a matrix
    /// expression -- not unconditionally here.
    template<class... expr> void computePattern(const expr &... args)
    {
        _computePattern(args...);
    }

    /// Initializes the pattern of the sparse matrix at boundary integrals.
    /// See computePattern() above re. where m_sparsity gets set.
    template<class... expr> void computePatternBdr(const bcRefList & BCs, const expr &... args)
    {
        _computePatternBdr(BCs, args...);
    }

    /// Initializes the pattern of the sparse matrix at boundary integrals.
    /// See computePattern() above re. where m_sparsity gets set.
    template<class... expr> void computePatternIfc(const ifContainer & iFaces, expr... args)
    {
        _computePatternIfc(iFaces, args...);
    }

    /// \brief Initializes the right-hand side vector only
    void initVector(const index_t numRhs = 1)
    {
        resetDimensions();
        m_rhs.setZero(numTestDofs(), numRhs);
    }

    /// Returns a block view of the system matrix, each block
    /// corresponding to a different space, or to different groups of
    /// dofs, in case of calar problems
    matBlockView matrixBlockView()
    {
        GISMO_ASSERT( m_vcol.back()->mapper.isFinalized(),
                      "initSystem() has not been called.");
        gsVector<index_t> rowSizes, colSizes;
        _blockDims(rowSizes, colSizes);
        matrix();
        return m_matrix.blockView(rowSizes,colSizes);
    }

    /// Returns a const block view of the system matrix, each block
    /// corresponding to a different space, or to different groups of
    /// dofs, in case of calar problems
    matConstBlockView matrixBlockView() const
    {
        GISMO_ASSERT( m_vcol.back()->mapper.isFinalized(),
                      "initSystem() has not been called.");
        gsVector<index_t> rowSizes, colSizes;
        _blockDims(rowSizes, colSizes);
        matrix();
        return m_matrix.blockView(rowSizes,colSizes);
    }

    /// Set the assembler options
    void setOptions(gsOptionList opt) { m_options = opt; } // gsOptionList opt
    // .swap(opt) todo

    /// \brief Adds the expressions \a args to the system matrix/rhs
    /// The arguments are considered as integrals over the whole domain
    /// \sa gsExprAssembler::setIntegrationElements
    template<class... expr> void assemble(const expr &... args);

    /// \brief Adds the expressions \a args to the system matrix/rhs
    /// The arguments are considered as integrals over the boundary
    /// parts in \a BCs
    template<class... expr> void assembleBdr(const bcRefList & BCs, expr&... args);

    template<class... expr> void assembleBdr(const bContainer & bnd, expr&... args);

    template<class... expr> void assembleIfc(const ifContainer & iFaces, expr... args);

    /** \brief Assembles the expressions \a args over the integration
        elements into an external \a sink instead of the internal matrix
        and right-hand side.

        A sink is any object with the member functions
        \code
        void addMatrix (const gsVector<index_t> & rows, const gsVector<index_t> & cols, const gsMatrix<T> & block);
        void addRhs    (const gsVector<index_t> & rows, const gsMatrix<T> & block);
        void addPattern(const gsVector<index_t> & rows, const gsVector<index_t> & cols); // computePattern*_into
        \endcode
        It receives one block per element (and per quadrature point, for
        expressions that are not element-wise), with the indices of the
        assembler's dof mappers. An index of -1 marks a row/column that is
        not free (eliminated or fixed); such entries must be ignored by the
        sink. With the elimination strategy, the contribution of the fixed
        dofs is passed to addRhs() by the assembler. The sink is called
        from all OpenMP threads concurrently and must synchronize itself.

        The default functions (assemble(), assembleBdr(), ...) run the same
        loops with a sink that writes into the internal fiber matrix and
        right-hand side. The *_into functions neither use nor allocate
        them, i.e. initSystem() is not required.
    */
    template<class Sink, class... expr> void assemble_into(Sink & sink, const expr &... args);

    /// \brief As assemble_into(), over the boundary parts in \a BCs
    template<class Sink, class... expr> void assembleBdr_into(Sink & sink, const bcRefList & BCs, expr&... args);

    /// \brief As assemble_into(), over the boundary sides \a bnd
    template<class Sink, class... expr> void assembleBdr_into(Sink & sink, const bContainer & bnd, expr&... args);

    /// \brief As assemble_into(), over the interfaces \a iFaces
    template<class Sink, class... expr> void assembleIfc_into(Sink & sink, const ifContainer & iFaces, expr... args);

    /// \brief Finite-difference Jacobian of \a residual with respect to \a u
    /// into \a sink (no elimination of fixed dofs)
    template<class Sink, class expr> void assembleJacobian_into(Sink & sink, const expr residual, solution & u);

    /** \brief Passes the sparsity pattern of the matrix expressions in \a
        args to \a sink (addPattern), one element at a time. Indices as in
        assemble_into(). Called from all OpenMP threads.
    */
    template<class Sink, class... expr> void computePattern_into(Sink & sink, const expr &... args);

    /// \brief As computePattern_into(), over the boundary parts in \a BCs
    template<class Sink, class... expr> void computePatternBdr_into(Sink & sink, const bcRefList & BCs, const expr &... args);

    /// \brief As computePattern_into(), over the boundary sides \a bnd
    template<class Sink, class... expr> void computePatternBdr_into(Sink & sink, const bContainer & bnd, const expr &... args);

    /// \brief As computePattern_into(), over the interfaces \a iFaces,
    /// including the couplings across the interface
    template<class Sink, class... expr> void computePatternIfc_into(Sink & sink, const ifContainer & iFaces, expr... args);
    /*
      template<class... expr> void collocate(expr... args);// eg. collocate(-ilapl(u), f)
    */

    void quPointsWeights(std::vector<gsMatrix<T> >&  cPoints, std::vector<gsVector<T> > & cWeights);

    /// \brief Assembles the Jacobian matrix of the expression \a args with
    // respect to the solution \a u
    template<class expr> void assembleJacobian(const expr residual, solution & u);

    template<class expr> void assembleJacobianIfc(const ifContainer & iFaces,
                                                  const expr residual, solution  u);

private:

    QuadratureRulePtr makeQuadratureRule(const gsBasis<T> & basis,
                                         index_t patch,
                                         short_t fixedDirection = -1) const
    {
        if (m_quadratureFactory)
        {
            QuadratureRulePtr rule =
                m_quadratureFactory(basis, m_options, patch, fixedDirection);
            GISMO_ENSURE(rule,
                         "Custom quadrature factory returned a null rule for patch "
                         << patch << ".");
            return rule;
        }

        return gsQuadrature::getPtr(basis, m_options, fixedDirection);
    }

    // Patterns of the internal fiber matrix
    template<class... expr> void _computePattern(const expr &... args);
    template<class... expr> void _computePatternBdr(const bcRefList & BCs, const expr &... args);
    template<class... expr> void _computePatternBdr(const bContainer & bnd, const expr &... args);
    template<class... expr> void _computePatternIfc(const ifContainer & iFaces, expr... args);

    // The assembly loops, shared by the internal (fiber matrix) and the
    // external (*_into) back ends
    template<class Sink, class... expr> void _volumeInto(Sink & sink, bool elim, const expr &... args);
    template<class Sink, class... expr> void _bdrInto(Sink & sink, bool elim, const bcRefList & BCs, expr&... args);
    template<class Sink, class... expr> void _bdrInto(Sink & sink, bool elim, const bContainer & bnd, expr&... args);
    template<class Sink, class... expr> void _ifcInto(Sink & sink, bool elim, const ifContainer & iFaces, expr... args);
    template<class Sink, class expr> void _jacobianInto(Sink & sink, const expr & residual, solution & u);
    template<class Sink, class expr> void _jacobianIfcInto(Sink & sink, const ifContainer & iFaces, const expr & residual, solution & u);
    template<class Sink, class... expr> void _patternVolume(Sink & sink, const expr &... args);
    template<class Sink, class... expr> void _patternBdr(Sink & sink, const bcRefList & BCs, const expr &... args);
    template<class Sink, class... expr> void _patternBdr(Sink & sink, const bContainer & bnd, const expr &... args);
    template<class Sink, class... expr> void _patternIfc(Sink & sink, const ifContainer & iFaces, expr... args);

    bool _elimination() const
    { return dirichlet::elimination==m_options.getInt("DirichletStrategy"); }

    // True if one of \a args is a matrix expression
    template<class... expr> static bool _hasMatrix(const expr &... args)
    {
        bool isMatrix = false;
        _checkMatrix CM(isMatrix);
        auto arg_tpl = std::make_tuple(args...);
        op_tuple(CM, arg_tpl);
        return isMatrix;
    }

    /// Serial warm-up for a subdomain-restricted domain's per-patch view
    /// (see gsIndexSubDomainPatchView in gsIndexSubDomain.h): touches its
    /// lazily-filled boundary-count cache once, from single-threaded code,
    /// before _computePatternBdr()/_computePatternIfc() enter their
    /// #pragma omp parallel region and construct boundary iterators on that
    /// view concurrently. Without this, first-touch of the cache would
    /// happen from multiple threads at once -- an unsynchronized write race
    /// (Fix 1). A no-op for domain types that don't restrict integration to
    /// a subdomain (subdomain() just returns the raw per-patch domain).
    void _warmupSubdomainViews(index_t patch, boxSide side)
    {
        m_exprdata->domain().subdomain(patch)->numElementsBdr(side);
    }

    void _blockDims(gsVector<index_t> & rowSizes,
                    gsVector<index_t> & colSizes)
    {
        if (1==m_vcol.size() && 1==m_vrow.size())
        {
            const gsDofMapper & dm = m_vcol.back()->mapper;
            rowSizes.resize(3);
            colSizes.resize(3);
            rowSizes[0]=colSizes[0] = dm.freeSize()-dm.coupledSize();
            rowSizes[1]=colSizes[1] = dm.coupledSize();
            rowSizes[2]=colSizes[2] = dm.boundarySize();
        }
        else
        {
            rowSizes.resize(m_vrow.size());
            for (index_t r = 0; r != rowSizes.size(); ++r) // for all row-blocks
                rowSizes[r] = m_vrow[r]->mapper.freeSize();
            colSizes.resize(m_vcol.size());
            for (index_t c = 0; c != colSizes.size(); ++c) // for all col-blocks
                colSizes[c] = m_vcol[c]->mapper.freeSize();
        }
    }

    /// \brief Reset the dimensions of all involved spaces.
    /// Called internally by the init* functions
    void resetDimensions();

    // Prints the expression to a text stream
    struct __printExpr
    {
        template <typename E> void operator() (const gismo::expr::_expr<E> & v)
        { v.print(gsInfo);gsInfo<<"\n"; }
    } _printExpr;

    // Checks validity of an expression
    struct __checkExpr
    {
        template <typename E> void operator() (const gsExprAssembler & ea,
                                               const gismo::expr::_expr<E> & ee)
        {
            auto u = ee.rowVar();
            auto v = ee.colVar();
#ifndef NDEBUG
            const bool m = E::isMatrix();
#endif
            GISMO_ASSERT(v.isValid(), "The row space is not valid");
            GISMO_ASSERT(!m || u.isValid(), "The column space is not valid");
            GISMO_ASSERT(m || (ea.numDofs()==ee.rhs().size()), "The right-hand side vector is not initialized");
        }
    };

    // Checks if an expression is a matrix
    struct _checkMatrix
    {
        bool & m_result;
        _checkMatrix( bool & _result) : m_result(_result) { }
        template <typename E> void operator() (const gismo::expr::_expr<E> &)
        { m_result |= E::isMatrix(); }
    };


    // Global indices of the active functions of space \a v (all
    // components), -1 for the ones that are not free
    static void _globalIndices(const expr::gsFeSpace<T> & v,
                               const gsMatrix<index_t> & act, index_t col,
                               index_t patch, gsVector<index_t> & idx)
    {
        const gsDofMapper & map = v.mapper();
        const index_t n = act.rows();
        idx.resize(n * v.dim());
        for (index_t r = 0; r != v.dim(); ++r)
            for (index_t i = 0; i != n; ++i)
            {
                const index_t ii = map.index(act(i,col), patch, r);
                idx[r*n+i] = map.is_free_index(ii) ? ii : -1;
            }
    }

    // The default sink: the internal fiber matrix and right-hand side
    struct _fiberSink
    {
        FiberMatrix & m_fmatrix;
        gsMatrix<T> & m_rhs;
#ifdef _OPENMP
        std::vector<omp_lock_t> * m_lock; // striped locks, for addPattern
#endif

        _fiberSink(FiberMatrix & _fmatrix, gsMatrix<T> & _rhs)
        : m_fmatrix(_fmatrix), m_rhs(_rhs)
#ifdef _OPENMP
        , m_lock(nullptr)
#endif
        { }

        // `omp atomic` only accepts scalar types: autodiff values are
        // accumulated in a critical section instead.
        template<typename U, typename std::enable_if<!gismo::is_autodiff_type<U>::value, int>::type = 0>
        static void accumulate(U & dst, const U & val)
        {
#           pragma omp atomic update
            dst += val;
        }
        template<typename U, typename std::enable_if<gismo::is_autodiff_type<U>::value, int>::type = 0>
        static void accumulate(U & dst, const U & val)
        {
#           pragma omp critical (gsExprAssembler_fiberSink)
            dst += val;
        }

        void addMatrix(const gsVector<index_t> & rows, const gsVector<index_t> & cols,
                       const gsMatrix<T> & block)
        {
            for (index_t j = 0; j != cols.size(); ++j)
            {
                const index_t jj = cols[j];
                if (jj < 0) continue;
                for (index_t i = 0; i != rows.size(); ++i)
                {
                    const index_t ii = rows[i];
                    if (ii < 0 || 0 == block(i,j)) continue;
                    accumulate(m_fmatrix.coeffRef(ii, jj), block(i,j));
                }
            }
        }

        void addRhs(const gsVector<index_t> & rows, const gsMatrix<T> & block)
        {
            for (index_t i = 0; i != rows.size(); ++i)
            {
                const index_t ii = rows[i];
                if (ii < 0) continue;
                for (index_t a = 0; a != block.cols(); ++a)
                {
                    accumulate(m_rhs(ii, a), block(i, a));
                }
            }
        }

        void addPattern(const gsVector<index_t> & rows, const gsVector<index_t> & cols)
        {
            for (index_t j = 0; j != cols.size(); ++j)
            {
                const index_t jj = cols[j];
                if (jj < 0) continue;
#ifdef _OPENMP
                GISMO_ASSERT(m_lock, "No lock pool for the pattern.");
                omp_lock_t & l = (*m_lock)[jj & (m_lock->size()-1)];
                omp_set_lock(&l);
#endif
                for (index_t i = 0; i != rows.size(); ++i)
                    if (rows[i] >= 0)
                        m_fmatrix.insertExplicitZero(rows[i], jj);
#ifdef _OPENMP
                omp_unset_lock(&l);
#endif
            }
        }
    };

    // Evaluates expressions and passes element blocks to a sink
    template<class Sink>
    struct _evalInto
    {
        Sink & m_sink;
        const gsVector<T> & m_quWeights;
        bool m_elim;
        gsMatrix<T> localMat, aux, rhsCorr;
        gsVector<index_t> rowIdx, colIdx;

        _evalInto(Sink & _sink, const gsVector<T> & _quWeights, bool _elim)
        : m_sink(_sink), m_quWeights(_quWeights), m_elim(_elim) { }

        template <typename E> void operator() (const gismo::expr::_expr<E> & ee)
        {
            GISMO_ASSERT(E::isMatrix() || E::isVector(), "Expecting a matrix or vector expression.");
            if ((ee.rowVar().data().flags & SAME_ELEMENT) &&
                ( E::isVector() || (ee.colVar().data().flags & SAME_ELEMENT) ) )
            {
                quadrature(ee, localMat);
                push<E::isMatrix()>(ee.rowVar(), ee.colVar(), 0, 0, m_elim);
            }
            else
            {
                const bool rowSame = ee.rowVar().data().flags & SAME_ELEMENT;
                const bool colSame = E::isVector() || (ee.colVar().data().flags & SAME_ELEMENT);
                const T * w = m_quWeights.data();
                for (index_t k = 0; k != m_quWeights.rows(); ++k)
                {
                    localMat.noalias() = (*(w++)) * ee.eval(k);
                    push<E::isMatrix()>(ee.rowVar(), ee.colVar(),
                                        rowSame ? 0 : k, colSame ? 0 : k, m_elim);
                }
            }
        }

        void operator() (const expr::_expr<expr::gsNullExpr<T> > &) {}

        template <typename E>
        inline void quadrature(const gismo::expr::_expr<E> & ee, gsMatrix<T> & lm)
        {
            const T * w = m_quWeights.data();
            lm.noalias() = (*w) * ee.eval(0);
            for (index_t k = 1; k != m_quWeights.rows(); ++k)
                lm.noalias() += (*(++w)) * ee.eval(k);
        }

        // Finite-difference Jacobian of the vector expression \a ee with
        // respect to \a u (fourth order central differences)
        template <typename E> void diff(const gismo::expr::_expr<E> & ee, solution & u)
        {
            GISMO_ASSERT(E::isVector(), "Expecting a vector expression.");
            static const T delta = 0.00001;

            const index_t sz = u.space().cardinality();
            localMat.setZero(sz, sz);

            for ( index_t c=0; c!= u.dim(); c++)
            {
                const index_t rls = c * u.data().actives.rows();     //local stride
                for ( index_t j = 0; j != sz/u.dim(); j++ )     // for all basis functions (col(j))
                {
                    const index_t jj = u.mapper().index(u.data().actives.at(j),
                                                        u.space().data().patchId, c);
                    if (u.mapper().is_free_index(jj) )
                    {
                        u.perturbLocal( delta  , jj, u.space().data().patchId);
                        quadrature(ee, aux);
                        localMat.col(rls+j) += 8 * aux;
                        u.perturbLocal( delta  , jj, u.space().data().patchId);
                        quadrature(ee, aux);
                        localMat.col(rls+j) -= aux;
                        u.perturbLocal(-3*delta, jj, u.space().data().patchId);
                        quadrature(ee, aux);
                        localMat.col(rls+j) -= 8 * aux;
                        u.perturbLocal( -delta , jj, u.space().data().patchId);
                        quadrature(ee, aux);
                        localMat.col(rls+j) += aux;
                        localMat.col(rls+j) /= 12*delta;
                        //Unperturb \a u
                        u.perturbLocal(2*delta, jj, u.space().data().patchId);
                    }
                }
            }
            push<true>(ee.rowVar(), ee.rowVar(), 0, 0, false);
        }

        template<bool isMatrix>
        void push(const expr::gsFeSpace<T> & v, const expr::gsFeSpace<T> & u,
                  index_t ra, index_t ca, bool elim)
        {
            _globalIndices(v, v.data().actives, ra, v.data().patchId, rowIdx);
            if (!isMatrix)
            {
                m_sink.addRhs(rowIdx, localMat);
                return;
            }
            _globalIndices(u, u.data().actives, ca, u.data().patchId, colIdx);
            GISMO_ASSERT( rowIdx.size()==localMat.rows() && colIdx.size()==localMat.cols(),
                          "Invalid local matrix (expected "<<rowIdx.size() <<"x"<< colIdx.size()
                          <<"), got\n" << localMat );
            m_sink.addMatrix(rowIdx, colIdx, localMat);

            if (!elim) return;
            // Symmetric treatment of eliminated dofs: rhs -= A_ib * g_b
            const gsDofMapper & colMap = u.mapper();
            const gsMatrix<T> & fixedDofs = u.fixedPart();
            const gsMatrix<index_t> & act = u.data().actives;
            const index_t nc = act.rows();
            bool any = false;
            for (index_t j = 0; j != colIdx.size(); ++j)
            {
                if (-1 != colIdx[j] || 0 == localMat.col(j).squaredNorm()) continue;
                const index_t jj = colMap.index(act(j % nc, ca), u.data().patchId, j / nc);
                if (colMap.is_boundary_index(jj) && !colMap.is_remote_index(jj))
                {
                    GISMO_ASSERT( colMap.boundarySize()==fixedDofs.size(),
                                  "Invalid values for fixed part " << colMap.boundarySize() <<" != "<< fixedDofs.size() );
                    if (!any) rhsCorr.setZero(localMat.rows(), 1);
                    rhsCorr.noalias() -= localMat.col(j) * fixedDofs.at(colMap.global_to_bindex(jj));
                    any = true;
                }
            }
            if (any) m_sink.addRhs(rowIdx, rhsCorr);
        }
    };

    // Passes the element-wise sparsity pattern to a sink. Either from the
    // active functions at a point of a patch (volume and boundary), or
    // from the evaluated data of the spaces (interfaces, where the row and
    // column spaces may live on different patches)
    template<class Sink>
    struct _patternInto
    {
        Sink & m_sink;
        const gsMatrix<T> & m_point;
        const index_t & m_patch;
        const bool m_fromData;
        gsMatrix<index_t> rowAct, colAct;
        gsVector<index_t> rowIdx, colIdx;

        _patternInto(Sink & _sink, const gsMatrix<T> & _point, const index_t & _patch,
                     bool fromData = false)
        : m_sink(_sink), m_point(_point), m_patch(_patch), m_fromData(fromData) { }

        template <typename E> void operator() (const gismo::expr::_expr<E> & ee)
        {
            if (!E::isMatrix()) return;
            const expr::gsFeSpace<T> & v = ee.rowVar();
            const expr::gsFeSpace<T> & u = ee.colVar();
            if (m_fromData)
            {
                _globalIndices(v, v.data().actives, 0, v.data().patchId, rowIdx);
                _globalIndices(u, u.data().actives, 0, u.data().patchId, colIdx);
            }
            else
            {
                v.source().piece(m_patch).active_into(m_point, rowAct);
                u.source().piece(m_patch).active_into(m_point, colAct);
                _globalIndices(v, rowAct, 0, m_patch, rowIdx);
                _globalIndices(u, colAct, 0, m_patch, colIdx);
            }
            m_sink.addPattern(rowIdx, colIdx);
        }

        void operator() (const expr::_expr<expr::gsNullExpr<T> > &) {}
    };

}; // gsExprAssembler

template<class T>
gsOptionList gsExprAssembler<T>::defaultOptions()
{
    gsOptionList opt;
    opt.addInt("DirichletValues"  , "Method for computation of Dirichlet DoF values [100..103]", 101);
    opt.addInt("DirichletStrategy", "Method for enforcement of Dirichlet BCs [11..14]", 11);
    opt.addReal("quA", "Number of quadrature points: quA*deg + quB; For patchRule: Regularity of the target space", 1.0  );
    opt.addInt ("quB", "Number of quadrature points: quA*deg + quB; For patchRule: Degree of the target space", 1    );
    opt.addReal("bdA", "Estimated nonzeros per column of the matrix: bdA*deg + bdB", 2.0  );
    opt.addInt ("bdB", "Estimated nonzeros per column of the matrix: bdA*deg + bdB", 1    );
    opt.addReal("bdO", "Overhead of sparse mem. allocation: (1+bdO)(bdA*deg + bdB) [0..1]", 0.333);
    opt.addInt ("quRule", "Quadrature rule used (1) Gauss-Legendre; (2) Gauss-Lobatto; (3) Patch-Rule",1);
    opt.addSwitch("overInt", "Apply over-integration on boundary elements or not?", false);
    opt.addSwitch("flipSide", "Flip side of interface where integration is performed.", false);
    opt.addSwitch("movingInterface", "Used in interface assembly when interface is not stationary.", false);
    opt.addSwitch("SameElement","Activates optimization if all quadrature points are located in the same element", true);
    opt.addSwitch("lazyMatrix","Allocate matrix columns on first use (memory-lean for subdomain assembly)", false);
    return opt;

    /// dirichlet treatment? elimination ????

    //storage of quadrature points, TP, ... non-linear assembly.

    //gsExpressions.h -> split ?

    //parallel interface assembly..

    // mpi assemly. ???
}

template<class T>
void gsExprAssembler<T>::setFixedDofVector(gsMatrix<T> & vals, short_t unk)
{
    gsMatrix<T> & fixedDofs = m_vcol[unk]->fixedDofs;
    fixedDofs.swap(vals);
    vals.resize(0, 0);
    // Assuming that the DoFs are already set by the user
    GISMO_ENSURE( fixedDofs.size() == m_vcol[unk]->mapper.boundarySize(),
                     "The Dirichlet DoFs were not provided correctly.");
}

template<class T>
void gsExprAssembler<T>::setFixedDofs(const gsMatrix<T> & coefMatrix, short_t unk, size_t patch)
{
    GISMO_ASSERT( m_options.getInt("DirichletValues") == dirichlet::user, "Incorrect options");

    expr::gsFeSpace<T> & u = *m_vcol[unk];
    //const index_t dirStr = m_options.getInt("DirichletStrategy");
    const gsMultiBasis<T> & mbasis = *dynamic_cast<const gsMultiBasis<T>* >(m_vcol[unk]->fs);

    const gsDofMapper & mapper = m_vcol[unk]->mapper;
//    const gsDofMapper & mapper =
//        dirichlet::elimination == dirStr ? u.mapper
//        : mbasis.getMapper(dirichlet::elimination,
//                           static_cast<iFace::strategy>(m_options.getInt("InterfaceStrategy")),
//                           bbc, u.id()) ;

    gsMatrix<T> & fixedDofs = m_vcol[unk]->fixedDofs;
    GISMO_ASSERT(fixedDofs.size() == m_vcol[unk]->mapper.boundarySize(),
                 "Fixed DoFs were not initialized.");

    // for every side with a Dirichlet BC
    typedef typename gsBoundaryConditions<T>::bcRefList bcRefList;
    for ( typename bcRefList::const_iterator it =  u.bc().dirichletBegin();
          it != u.bc().dirichletEnd()  ; ++it )
    {
        const index_t com = it->unkComponent();
        const index_t k = it->patch();
        if ( k == patch )
        {
            // Get indices in the patch on this boundary
            const gsMatrix<index_t> boundary =
                    mbasis[k].boundary(it->side());

            //gsInfo <<"Setting the value for: "<< boundary.transpose() <<"\n";

            for (index_t i=0; i!= boundary.size(); ++i)
            {
                // Note: boundary.at(i) is the patch-local index of a
                // control point on the patch
                const index_t ii  = mapper.bindex( boundary.at(i) , k, com );

                fixedDofs.at(ii) = coefMatrix(boundary.at(i), com);
            }
        }
    }
} // setFixedDofs


template<class T> void gsExprAssembler<T>::resetDimensions()
{
    if (!m_vcol.front()->valid()) m_vcol.front()->init();
    if (!m_vrow.front()->valid()) m_vrow.front()->init();
    for (size_t i = 1; i!=m_vcol.size(); ++i)
    {
        if (!m_vcol[i]->valid()) m_vcol[i]->init();
        m_vcol[i]->mapper.setShift(m_vcol[i-1]->mapper.firstIndex() +
                                   m_vcol[i-1]->mapper.freeSize() );

        if ( i<m_vrow.size() && m_vcol[i] != m_vrow[i] )
        {
            if (!m_vrow[i]->valid()) m_vrow[i]->init();
            m_vrow[i]->mapper.setShift(m_vrow[i-1]->mapper.firstIndex() +
                                       m_vrow[i-1]->mapper.freeSize() );
        }
    }
}

template<size_t I, class op, typename... Ts>
void op_tuple_impl (op & _op, const std::tuple<Ts...> &tuple)
{
    _op(std::get<I>(tuple));
    if (I + 1 < sizeof... (Ts))
        op_tuple_impl<(I+1 < sizeof... (Ts) ? I+1 : I)> (_op, tuple);
}

template<class op, typename... Ts>
void op_tuple (op & _op, const std::tuple<Ts...> &tuple)
{ op_tuple_impl<0>(_op,tuple); }

// ===================================================================
// Patterns of the internal fiber matrix
// ===================================================================

template<class T>
template<class... expr>
void gsExprAssembler<T>::_computePattern(const expr &... args)
{
    GISMO_ASSERT(m_fmatrix.cols()==numDofs(), "System not initialized, matrix.cols() = "<<m_fmatrix.cols()<<"!="<<numDofs()<<" = numDofs()");

    // No matrix expression, no pattern
    if (!_hasMatrix(args...)) return;

    _fiberSink fs(m_fmatrix, m_rhs);
#ifdef _OPENMP
    // Striped lock pool (power-of-two cap 1<<16) instead of one lock per
    // dof. Threads hold at most one lock at a time, so no deadlock.
    std::vector<omp_lock_t> lock(patternLockPoolSize(numDofs()));
    for (auto & l : lock) omp_init_lock(&l);
    fs.m_lock = &lock;
#endif

    _patternVolume(fs, args...);

#ifdef _OPENMP
    for (auto & l : lock) omp_destroy_lock(&l);
#endif
    // set only when a pattern was actually computed
    m_sparsity |= 1;
}

template<class T>
template<class... expr>
void gsExprAssembler<T>::_computePatternBdr(const bcRefList & BCs, const expr &... args)
{
    GISMO_ASSERT(m_fmatrix.cols()==numDofs(), "System not initialized, matrix.cols() = "<<m_fmatrix.cols()<<"!="<<numDofs()<<" = numDofs()");
    if ( BCs.empty() || 0==numDofs() || !_hasMatrix(args...) ) return;

    _fiberSink fs(m_fmatrix, m_rhs);
#ifdef _OPENMP
    std::vector<omp_lock_t> lock(patternLockPoolSize(numDofs()));
    for (auto & l : lock) omp_init_lock(&l);
    fs.m_lock = &lock;
#endif

    _patternBdr(fs, BCs, args...);

#ifdef _OPENMP
    for (auto & l : lock) omp_destroy_lock(&l);
#endif
    // Shared with the bContainer overload: both mean "boundary pattern
    // computed" (both walk boundary sides through the same domain)
    m_sparsity |= 2;
}

template<class T>
template<class... expr>
void gsExprAssembler<T>::_computePatternBdr(const bContainer & bnd, const expr &... args)
{
    GISMO_ASSERT(m_fmatrix.cols()==numDofs(), "System not initialized, matrix.cols() = "<<m_fmatrix.cols()<<"!="<<numDofs()<<" = numDofs()");
    if ( bnd.size()==0 || 0==numDofs() || !_hasMatrix(args...) ) return;

    _fiberSink fs(m_fmatrix, m_rhs);
#ifdef _OPENMP
    std::vector<omp_lock_t> lock(patternLockPoolSize(numDofs()));
    for (auto & l : lock) omp_init_lock(&l);
    fs.m_lock = &lock;
#endif

    _patternBdr(fs, bnd, args...);

#ifdef _OPENMP
    for (auto & l : lock) omp_destroy_lock(&l);
#endif
    m_sparsity |= 2;
}

template<class T>
template<class... expr>
void gsExprAssembler<T>::_computePatternIfc(const ifContainer & iFaces, expr... args)
{
    GISMO_ASSERT(m_fmatrix.cols()==numDofs(), "System not initialized");
    if ( iFaces.empty() || 0==numDofs() || !_hasMatrix(args...) ) return;

    _fiberSink fs(m_fmatrix, m_rhs);
#ifdef _OPENMP
    std::vector<omp_lock_t> lock(patternLockPoolSize(numDofs()));
    for (auto & l : lock) omp_init_lock(&l);
    fs.m_lock = &lock;
#endif

    _patternIfc(fs, iFaces, args...);

#ifdef _OPENMP
    for (auto & l : lock) omp_destroy_lock(&l);
#endif
    m_sparsity |= 4;
}

// ===================================================================
// Assembly into the internal fiber matrix and right-hand side
// ===================================================================

template<class T>
template<class... expr>
void gsExprAssembler<T>::assemble(const expr &... args)
{
    GISMO_ASSERT(m_fmatrix.cols()==numDofs(), "System not initialized, matrix.cols() = "<<m_fmatrix.cols()<<"!="<<numDofs()<<" = numDofs()");

    if ((m_sparsity & 1) == 0)
        this->_computePattern(args...);

    m_modified |= _hasMatrix(args...);
    _fiberSink fs(m_fmatrix, m_rhs);
    _volumeInto(fs, _elimination(), args...);
}

template<class T>
template<class... expr>
void gsExprAssembler<T>::assembleBdr(const bcRefList & BCs, expr&... args)
{
    GISMO_ASSERT(m_fmatrix.cols()==numDofs(), "System not initialized");
    if ( BCs.empty() || 0==numDofs() ) return;

    if ((m_sparsity & 2) == 0)
        this->_computePatternBdr(BCs, args...);

    m_modified |= _hasMatrix(args...);
    _fiberSink fs(m_fmatrix, m_rhs);
    _bdrInto(fs, true, BCs, args...);
}

template<class T>
template<class... expr>
void gsExprAssembler<T>::assembleBdr(const bContainer & bnd, expr&... args)
{
    GISMO_ASSERT(m_fmatrix.cols()==numDofs(), "System not initialized");
    if ( bnd.size()==0 || 0==numDofs() ) return;

    if ((m_sparsity & 2) == 0)
        this->_computePatternBdr(bnd, args...);

    m_modified |= _hasMatrix(args...);
    _fiberSink fs(m_fmatrix, m_rhs);
    _bdrInto(fs, true, bnd, args...);
}

template<class T> template<class... expr>
void gsExprAssembler<T>::assembleIfc(const ifContainer & iFaces, expr... args)
{
    GISMO_ASSERT(m_fmatrix.cols()==numDofs(), "System not initialized");

    if ((m_sparsity & 4) == 0)
        this->_computePatternIfc(iFaces, args...);

    m_modified |= _hasMatrix(args...);
    _fiberSink fs(m_fmatrix, m_rhs);
    _ifcInto(fs, true, iFaces, args...);
}

template<class T> template<class expr>
void gsExprAssembler<T>::assembleJacobian(const expr residual, solution & u)
{
    GISMO_ASSERT(m_fmatrix.cols()==numDofs(), "System not initialized");
    GISMO_ASSERT(expr::isVector(), "Expecting a vector expression.");

    clearMatrix();
    clearRhs();
    m_modified = true;
    _fiberSink fs(m_fmatrix, m_rhs);
    _jacobianInto(fs, residual, u);
}

template<class T> template<class expr>
void gsExprAssembler<T>::assembleJacobianIfc(const ifContainer & iFaces,
                                             const expr residual, solution  u)
{
    GISMO_ASSERT(m_fmatrix.cols()==numDofs(), "System not initialized");
    GISMO_ASSERT(expr::isVector(), "Expecting a vector expression.");

    m_modified = true;
    _fiberSink fs(m_fmatrix, m_rhs);
    _jacobianIfcInto(fs, iFaces, residual, u);
}

// ===================================================================
// Assembly into an external sink
// ===================================================================

template<class T>
template<class Sink, class... expr>
void gsExprAssembler<T>::assemble_into(Sink & sink, const expr &... args)
{
    resetDimensions();
    _volumeInto(sink, _elimination(), args...);
}

template<class T>
template<class Sink, class... expr>
void gsExprAssembler<T>::assembleBdr_into(Sink & sink, const bcRefList & BCs, expr&... args)
{
    if ( BCs.empty() ) return;
    resetDimensions();
    _bdrInto(sink, true, BCs, args...);
}

template<class T>
template<class Sink, class... expr>
void gsExprAssembler<T>::assembleBdr_into(Sink & sink, const bContainer & bnd, expr&... args)
{
    if ( bnd.size()==0 ) return;
    resetDimensions();
    _bdrInto(sink, true, bnd, args...);
}

template<class T>
template<class Sink, class... expr>
void gsExprAssembler<T>::assembleIfc_into(Sink & sink, const ifContainer & iFaces, expr... args)
{
    resetDimensions();
    _ifcInto(sink, true, iFaces, args...);
}

template<class T>
template<class Sink, class expr>
void gsExprAssembler<T>::assembleJacobian_into(Sink & sink, const expr residual, solution & u)
{
    resetDimensions();
    _jacobianInto(sink, residual, u);
}

template<class T>
template<class Sink, class... expr>
void gsExprAssembler<T>::computePattern_into(Sink & sink, const expr &... args)
{
    resetDimensions();
    _patternVolume(sink, args...);
}

template<class T>
template<class Sink, class... expr>
void gsExprAssembler<T>::computePatternBdr_into(Sink & sink, const bcRefList & BCs, const expr &... args)
{
    resetDimensions();
    _patternBdr(sink, BCs, args...);
}

template<class T>
template<class Sink, class... expr>
void gsExprAssembler<T>::computePatternBdr_into(Sink & sink, const bContainer & bnd, const expr &... args)
{
    resetDimensions();
    _patternBdr(sink, bnd, args...);
}

template<class T>
template<class Sink, class... expr>
void gsExprAssembler<T>::computePatternIfc_into(Sink & sink, const ifContainer & iFaces, expr... args)
{
    resetDimensions();
    _patternIfc(sink, iFaces, args...);
}

// ===================================================================
// The loops
// ===================================================================

template<class T>
template<class Sink, class... expr>
void gsExprAssembler<T>::_volumeInto(Sink & sink, bool elim, const expr &... args)
{
#pragma omp parallel
{
    auto arg_tpl = std::make_tuple(args...);
    m_exprdata->parse(arg_tpl);
    if (m_options.askSwitch("SameElement",true)) m_exprdata->activateFlags(SAME_ELEMENT);

    _evalInto<Sink> ee(sink, m_exprdata->weights(), elim);

    typename gsQuadRule<T>::uPtr QuRule;
    index_t QuPatch = -1;

    for ( auto & elem : m_exprdata->domain().allElements() )
    {
        if (QuPatch!=elem.patchIndex())
        {
            QuPatch = elem.patchIndex();
            // get Degree of the domain
            QuRule = makeQuadratureRule(this->trialSpace(0).source().basis(QuPatch), QuPatch);
        }

        // Map the Quadrature rule to the element
        QuRule->mapTo( elem.lowerCorner(), elem.upperCorner(),
                       m_exprdata->points(), m_exprdata->weights());

        if (m_exprdata->points().cols()==0)
            continue;

        m_exprdata->precompute( QuPatch );

        // Assemble contributions of the element
        op_tuple(ee, arg_tpl);
    }
}//omp parallel
}

template<class T>
template<class Sink, class... expr>
void gsExprAssembler<T>::_bdrInto(Sink & sink, bool elim, const bcRefList & BCs, expr&... args)
{
    m_exprdata->setMutSource(*BCs.front().get().function()); //initialize once
    auto arg_tpl = std::make_tuple(args...);
    m_exprdata->parse(arg_tpl);
    if (m_options.askSwitch("SameElement",true)) m_exprdata->activateFlags(SAME_ELEMENT);

    typename gsQuadRule<T>::uPtr QuRule;
    _evalInto<Sink> ee(sink, m_exprdata->weights(), elim);

    for (typename bcRefList::const_iterator iit = BCs.begin(); iit!= BCs.end(); ++iit)
    {
        const boundary_condition<T> * it = &iit->get();

        QuRule = makeQuadratureRule(this->trialSpace(0).source().basis(it->patch()),
                                    it->patch(), it->side().direction());

        // Update boundary function source
        m_exprdata->setMutSource(*it->function());

        typename gsBasis<T>::domainIter domIt =
            m_exprdata->domain().subdomain(it->patch())->beginBdr(it->side());
        typename gsBasis<T>::domainIter domItEnd =
            m_exprdata->domain().subdomain(it->patch())->endBdr(it->side());

        for (; domIt < domItEnd; ++domIt )
        {
            QuRule->mapTo( domIt.lowerCorner(), domIt.upperCorner(),
                           m_exprdata->points(), m_exprdata->weights());
            if (m_exprdata->points().cols()==0)
                continue;
            m_exprdata->precompute(it->patch(), it->side());
            op_tuple(ee, arg_tpl);
        }
    }
}

template<class T>
template<class Sink, class... expr>
void gsExprAssembler<T>::_bdrInto(Sink & sink, bool elim, const bContainer & bnd, expr&... args)
{
    auto arg_tpl = std::make_tuple(args...);
    m_exprdata->parse(arg_tpl);

    typename gsQuadRule<T>::uPtr QuRule;
    _evalInto<Sink> ee(sink, m_exprdata->weights(), elim);

    for (gsBoxTopology::const_biterator it = bnd.begin(); it != bnd.end(); ++it )
    {
        QuRule = makeQuadratureRule(this->trialSpace(0).source().basis(it->patch),
                                    it->patch, it->side().direction());

        typename gsBasis<T>::domainIter domIt =
            m_exprdata->domain().subdomain(it->patch)->beginBdr(it->side());
        typename gsBasis<T>::domainIter domItEnd =
            m_exprdata->domain().subdomain(it->patch)->endBdr(it->side());

        for (; domIt<domItEnd; ++domIt )
        {
            QuRule->mapTo( domIt.lowerCorner(), domIt.upperCorner(),
                           m_exprdata->points(), m_exprdata->weights());
            if (m_exprdata->points().cols()==0)
                continue;
            m_exprdata->precompute(it->patch, it->side());
            op_tuple(ee, arg_tpl);
        }
    }
}

template<class T>
template<class Sink, class... expr>
void gsExprAssembler<T>::_ifcInto(Sink & sink, bool elim, const ifContainer & iFaces, expr... args)
{
    typedef typename gsFunction<T>::uPtr ifacemap;

    auto arg_tpl = std::make_tuple(args...);
    m_exprdata->parse(arg_tpl);
    if (m_options.askSwitch("SameElement",true)) m_exprdata->activateFlags(SAME_ELEMENT); //note: SAME_ELEMENT is 0 at the opposite/mirrored patch

    typename gsQuadRule<T>::uPtr QuRule;
    _evalInto<Sink> ee(sink, m_exprdata->weights(), elim);

    const bool flipSide = m_options.askSwitch("flipSide", false);

    ifacemap interfaceMap;
    for (gsBoxTopology::const_iiterator it = iFaces.begin(); it != iFaces.end(); ++it )
    {
        // If flipSide switch is enabled, then the integration will be
        // performed on the opposite side of the interface
        const boundaryInterface & iFace =  flipSide ? it->getInverse() : *it;
        const index_t patch1 = iFace.first() .patch;
        const index_t patch2 = iFace.second().patch;

        if (iFace.type() == interaction::conforming)
            interfaceMap = gsAffineFunction<T>::make( iFace.dirMap(), iFace.dirOrientation(),
                                                      this->trialSpace(0).source().basis(patch1).support(),
                                                      this->trialSpace(0).source().basis(patch2).support() );
        else
            interfaceMap = gsCPPInterface<T>::make(getGeometryMap(), iFace);

        QuRule = makeQuadratureRule(this->trialSpace(0).source().basis(patch1),
                                    patch1, iFace.first().side().direction());

        typename gsBasis<T>::domainIter domIt =
            m_exprdata->domain().subdomain(patch1)->beginBdr(iFace.first().side());
        typename gsBasis<T>::domainIter domItEnd =
            m_exprdata->domain().subdomain(patch1)->endBdr(iFace.first().side());

        for (; domIt<domItEnd; ++domIt)
        {
            QuRule->mapTo( domIt.lowerCorner(), domIt.upperCorner(),
                           m_exprdata->points(), m_exprdata->weights());
            interfaceMap->eval_into(m_exprdata->points(), m_exprdata->pointsIfc());

            if (m_exprdata->points().cols()==0)
                continue;

            m_exprdata->precompute(iFace);
            op_tuple(ee, arg_tpl);
        }
    }
}

template<class T>
template<class Sink, class expr>
void gsExprAssembler<T>::_jacobianInto(Sink & sink, const expr & residual, solution & u)
{
#pragma omp parallel
{
    m_exprdata->parse(residual, u);
    if (m_options.askSwitch("SameElement",true)) m_exprdata->activateFlags(SAME_ELEMENT);

    _evalInto<Sink> ee(sink, m_exprdata->weights(), false);

    typename gsQuadRule<T>::uPtr QuRule;
    index_t QuPatch = -1;

    for ( auto & elem : m_exprdata->domain().allElements() )
    {
        if (QuPatch!=elem.patchIndex())
        {
            QuPatch = elem.patchIndex();
            // get Degree of the domain
            QuRule = makeQuadratureRule(this->trialSpace(0).source().basis(QuPatch), QuPatch);
        }

        QuRule->mapTo( elem.lowerCorner(), elem.upperCorner(),
                       m_exprdata->points(), m_exprdata->weights());

        if (m_exprdata->points().cols()==0)
            continue;

        m_exprdata->precompute(QuPatch);

        // the perturbation of u is shared by the threads
#       pragma omp critical (assemble_fdiffs)
        ee.diff(residual, u);
    }
}//omp parallel
}

template<class T>
template<class Sink, class expr>
void gsExprAssembler<T>::_jacobianIfcInto(Sink & sink, const ifContainer & iFaces,
                                          const expr & residual, solution & u)
{
    m_exprdata->parse(residual, u);
    if (m_options.askSwitch("SameElement",true)) m_exprdata->activateFlags(SAME_ELEMENT);

    typename gsQuadRule<T>::uPtr QuRule;
    _evalInto<Sink> ee(sink, m_exprdata->weights(), false);
    const bool flipSide = m_options.askSwitch("flipSide", false);
    const bool movingInterface = m_options.askSwitch("movingInterface", false);

    for (gsBoxTopology::const_iiterator it = iFaces.begin(); it != iFaces.end(); ++it )
    {
        const boundaryInterface & iFace =  flipSide ? it->getInverse() : *it;
        const index_t patch1 = iFace.first() .patch;

        gsCPPInterface<T> interfaceMap(getGeometryMap(), iFace);

        QuRule = makeQuadratureRule(this->trialSpace(0).source().basis(patch1),
                                    patch1, iFace.first().side().direction());

        typename gsBasis<T>::domainIter domIt =
            m_exprdata->domain().subdomain(patch1)->beginBdr(iFace.first().side());
        typename gsBasis<T>::domainIter domItEnd =
            m_exprdata->domain().subdomain(patch1)->endBdr(iFace.first().side());

        for (; domIt<domItEnd; ++domIt )
        {
            QuRule->mapTo( domIt.lowerCorner(), domIt.upperCorner(),
                           m_exprdata->points(), m_exprdata->weights());
            interfaceMap.eval_into(m_exprdata->points(), m_exprdata->pointsIfc());

            if (m_exprdata->points().cols()==0)
                continue;

            m_exprdata->precompute(iFace);

            if (!movingInterface)
            {
                ee.diff(residual, u);
            }
            else // For Moving Geometry Maps: the interface map depends on u
            {
                static const T delta = 0.00001;
                const index_t sz = residual.eval(0).rows(); // depends on .left()/.right()
                ee.localMat.setZero(sz, sz);
                auto & rowVar = residual.rowVar();

                // Perturb u by \a d, update the interface and integrate
                const auto perturbed = [&](T d, index_t jj, gsMatrix<T> & out)
                {
                    u.perturbLocal(d, jj, rowVar.data().patchId);
                    interfaceMap.updateBdr();
                    interfaceMap.eval_into(m_exprdata->points(), m_exprdata->pointsIfc());
                    m_exprdata->precompute(iFace);
                    ee.quadrature(residual, out);
                };

                for ( index_t c=0; c!= u.dim(); c++)
                {
                    const index_t rls = c * rowVar.data().actives.rows();     //local stride
                    for ( index_t j = 0; j != sz/u.dim(); j++ )     // for all basis functions (col(j))
                    {
                        const index_t jj = u.mapper().index(rowVar.data().actives(j),
                                                            rowVar.data().patchId, c);
                        if (rowVar.mapper().is_free_index(jj) )
                        {
                            perturbed(   delta, jj, ee.aux); ee.localMat.col(rls+j) += 8 * ee.aux;
                            perturbed(   delta, jj, ee.aux); ee.localMat.col(rls+j) -=     ee.aux;
                            perturbed(-3*delta, jj, ee.aux); ee.localMat.col(rls+j) -= 8 * ee.aux;
                            perturbed(  -delta, jj, ee.aux); ee.localMat.col(rls+j) +=     ee.aux;
                            ee.localMat.col(rls+j) /= 12*delta;
                            //Unperturb \a u
                            u.perturbLocal(2*delta, jj, rowVar.data().patchId);
                            interfaceMap.updateBdr();
                        }
                    }
                }
                ee.template push<true>(residual.rowVar(), residual.rowVar(), 0, 0, false);
            }
        }
    }
}

template<class T>
template<class Sink, class... expr>
void gsExprAssembler<T>::_patternVolume(Sink & sink, const expr &... args)
{
#pragma omp parallel
{
    auto arg_tpl = std::make_tuple(args...);
    m_exprdata->parsePattern(arg_tpl);
    index_t patch = 0;
    _patternInto<Sink> pp(sink, m_exprdata->points(), patch);

    for ( auto & elem : m_exprdata->domain().allElements() )
    {
        m_exprdata->points() = elem.centerPoint();
        patch = elem.patchIndex();
        op_tuple(pp, arg_tpl);
    }
}//omp parallel
}

template<class T>
template<class Sink, class... expr>
void gsExprAssembler<T>::_patternBdr(Sink & sink, const bcRefList & BCs, const expr &... args)
{
    // Touch the (lazily filled) boundary caches of subdomain views
    // serially, before the parallel region
    for (typename bcRefList::const_iterator iit = BCs.begin(); iit != BCs.end(); ++iit)
        _warmupSubdomainViews(iit->get().patch(), iit->get().side());

#pragma omp parallel
{
    auto arg_tpl = std::make_tuple(args...);
    m_exprdata->parsePattern(arg_tpl);
    index_t patch = 0;
    _patternInto<Sink> pp(sink, m_exprdata->points(), patch);

    for (typename bcRefList::const_iterator iit = BCs.begin(); iit!= BCs.end(); ++iit)
    {
        const boundary_condition<T> * it = &iit->get();
        typename gsBasis<T>::domainIter domIt =
            m_exprdata->domain().subdomain(it->patch())->beginBdr(it->side());
        typename gsBasis<T>::domainIter domItEnd =
            m_exprdata->domain().subdomain(it->patch())->endBdr(it->side());

        for (; domIt<domItEnd; ++domIt )
        {
#           pragma omp single nowait
            {
                patch = it->patch();
                m_exprdata->points() = domIt.centerPoint();
                op_tuple(pp, arg_tpl);
            }
        }
    }
}//omp parallel
}

template<class T>
template<class Sink, class... expr>
void gsExprAssembler<T>::_patternBdr(Sink & sink, const bContainer & bnd, const expr &... args)
{
    for (gsBoxTopology::const_biterator it = bnd.begin(); it != bnd.end(); ++it)
        _warmupSubdomainViews(it->patch, it->side());

#pragma omp parallel
{
    auto arg_tpl = std::make_tuple(args...);
    m_exprdata->parsePattern(arg_tpl);
    index_t patch = 0;
    _patternInto<Sink> pp(sink, m_exprdata->points(), patch);

    for (gsBoxTopology::const_biterator it = bnd.begin(); it != bnd.end(); ++it)
    {
        typename gsBasis<T>::domainIter domIt =
            m_exprdata->domain().subdomain(it->patch)->beginBdr(it->side());
        typename gsBasis<T>::domainIter domItEnd =
            m_exprdata->domain().subdomain(it->patch)->endBdr(it->side());

        for (; domIt<domItEnd; ++domIt )
        {
#           pragma omp single nowait
            {
                patch = it->patch;
                m_exprdata->points() = domIt.centerPoint();
                op_tuple(pp, arg_tpl);
            }
        }
    }
}//omp parallel
}

template<class T>
template<class Sink, class... expr>
void gsExprAssembler<T>::_patternIfc(Sink & sink, const ifContainer & iFaces, expr... args)
{
    typedef typename gsFunction<T>::uPtr ifacemap;
    const bool flipSide = m_options.askSwitch("flipSide", false);

    for (gsBoxTopology::const_iiterator it = iFaces.begin(); it != iFaces.end(); ++it)
    {
        const boundaryInterface & iFace = flipSide ? it->getInverse() : *it;
        _warmupSubdomainViews(iFace.first().patch, iFace.first().side());
    }

    // Create the interface mirror serially; inside the region below every
    // thread would otherwise race to create it.
    m_exprdata->initializeIface();

#pragma omp parallel
{
    auto arg_tpl = std::make_tuple(args...);
    m_exprdata->parsePattern(arg_tpl);
    index_t patch = 0;
    // The row and column spaces may live on different sides of the
    // interface: take the active functions from the evaluated data
    _patternInto<Sink> pp(sink, m_exprdata->points(), patch, true);

    ifacemap interfaceMap;
    for (gsBoxTopology::const_iiterator it = iFaces.begin(); it != iFaces.end(); ++it )
    {
        const boundaryInterface & iFace =  flipSide ? it->getInverse() : *it;
        const index_t patch1 = iFace.first() .patch;
        const index_t patch2 = iFace.second().patch;

        const gsBasis<T> & basis1 = this->trialSpace(0).source().basis(patch1);
        const gsBasis<T> & basis2 = this->trialSpace(0).source().basis(patch2);
        if (iFace.type() == interaction::conforming)
            interfaceMap = gsAffineFunction<T>::make( iFace.dirMap(), iFace.dirOrientation(),
                                                      basis1.support(), basis2.support() );
        else
            interfaceMap = gsCPPInterface<T>::make(getGeometryMap(), iFace);

        typename gsBasis<T>::domainIter domIt =
            m_exprdata->domain().subdomain(patch1)->beginBdr(iFace.first().side());
        typename gsBasis<T>::domainIter domItEnd =
            m_exprdata->domain().subdomain(patch1)->endBdr(iFace.first().side());

        for (; domIt<domItEnd; ++domIt )
        {
#           pragma omp single nowait
            {
                m_exprdata->points() = domIt.centerPoint();
                interfaceMap->eval_into(m_exprdata->points(), m_exprdata->pointsIfc());
                m_exprdata->precompute(iFace);
                op_tuple(pp, arg_tpl);
            }
        }
    }
}//omp parallel
}

template<class T>
void gsExprAssembler<T>::quPointsWeights(std::vector<gsMatrix<T> >&  cPoints, std::vector<gsVector<T> > & cWeights)
{
    GISMO_ASSERT(m_fmatrix.cols()==numDofs(), "System not initialized, matrix.cols() = "<<m_fmatrix.cols()<<"!="<<numDofs()<<" = numDofs()");

    //bool changeQuadrature = !m_options.askSwitch("SameQuadrature",true);
//#pragma omp parallel
{
    typename gsQuadRule<T>::uPtr QuRule; // Quadrature rule
    cPoints.resize( m_exprdata->domain().nPieces() );
    cWeights.resize( m_exprdata->domain().nPieces() );

    index_t count = 0;
    for (unsigned patchInd = 0; patchInd < m_exprdata->domain().nPieces(); ++patchInd)
    {
        auto & bb = this->trialSpace(0).source().basis(patchInd);

        QuRule = makeQuadratureRule(bb, patchInd);
        const index_t numNodes = QuRule->numNodes();

        // @hverhelst: THIS ASSUMES SAME NUMBER OF QUNODES PER ELEMENT
        const index_t sz = bb.numElements() * numNodes;
        cPoints[patchInd].resize(bb.domainDim(), sz );
        cWeights[patchInd].resize( sz );

        // Start iteration over elements of patchInd
        for ( auto & elem : bb.domain()->allElements() ) //todo: parallelize
        {
            // Map the Quadrature rule to the element
            QuRule->mapTo( elem.lowerCorner(), elem.upperCorner(),
                           m_exprdata->points(), m_exprdata->weights());

            cWeights[patchInd].segment(count, numNodes) = m_exprdata->weights();
            cPoints[patchInd].middleCols(count, numNodes) = m_exprdata->points();
            count += numNodes;
        }
    }
}//omp parallel

}



} //namespace gismo
