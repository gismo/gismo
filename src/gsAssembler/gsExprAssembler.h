/** @file gsExprAssembler.h

    @brief Generic expressions matrix assembly

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris
*/

#pragma once

#include <gsUtils/gsPointGrid.h>
#include <gsAssembler/gsQuadrature.h>
#include <gsAssembler/gsExprHelper.h>

namespace gismo
{

/**
   Assembler class for generating matrices and right-hand sides based
   on isogeometric expressions
*/
template<class T>
class gsExprAssembler
{
private:
    typename gsExprHelper<T>::Ptr m_exprdata;

    gsOptionList m_options;

    expr::gsFeElement<T> m_element;

    gsSparseMatrix<T> m_matrix;
    gsMatrix<T>       m_rhs;

    std::vector<expr::gsFeSpace<T>*> m_vrow;
    std::vector<expr::gsFeSpace<T>*> m_vcol;

    typedef typename gsExprHelper<T>::nullExpr    nullExpr;

public:

    typedef typename gsSparseMatrix<T>::BlockView matBlockView;

    typedef typename gsSparseMatrix<T>::constBlockView matConstBlockView;

    typedef typename gsBoundaryConditions<T>::bcRefList   bcRefList;
    typedef typename gsBoundaryConditions<T>::bcContainer bcContainer;
    //typedef typename gsBoundaryConditions<T>::ppContainer ifContainer;
    typedef gsBoxTopology::ifContainer ifContainer;

    typedef typename gsExprHelper<T>::element     element;     ///< Current element
    typedef typename gsExprHelper<T>::geometryMap geometryMap; ///< Geometry map type

    typedef typename gsExprHelper<T>::variable    variable;    ///< Variable type
    typedef typename gsExprHelper<T>::space       space;       ///< Space type
    typedef typename expr::gsFeSolution<T>        solution;    ///< Solution type

    /*
    typedef typename gsExprHelper<T>::function    function;    ///< Variable type
    typedef typename gsExprHelper<T>::variable    variable;    ///< Space type
    typedef typename expr::gsFeSolution<T>        solution;    ///< Solution type
    */
public:

    void cleanUp()
    {
        m_exprdata->clean();
    }

    /// Constructor
    /// \param _rBlocks Number of spaces for test functions
    /// \param _cBlocks Number of spaces for solution variables
    gsExprAssembler(index_t _rBlocks = 1, index_t _cBlocks = 1)
    : m_exprdata(gsExprHelper<T>::make()), m_options(defaultOptions()),
      m_vrow(_rBlocks,nullptr), m_vcol(_cBlocks,nullptr)
    { }

    // The copy constructor replicates the same environemnt but does
    // not copy any matrix data

    /// @brief Returns the list of default options for assembly
    static gsOptionList defaultOptions();

    /// Returns the number of degrees of freedom (after initialization)
    index_t numDofs() const
    {
        GISMO_ASSERT( m_vcol.back()->mapper().isFinalized(),
                      "gsExprAssembler::numDofs() says: initSystem() has not been called.");
        return m_vcol.back()->mapper().firstIndex() +
	  m_vcol.back()->mapper().freeSize();
    }

    /// Returns the number of test functions (after initialization)
    index_t numTestDofs() const
    {
        GISMO_ASSERT( m_vrow.back()->mapper().isFinalized(),
                      "initSystem() has not been called.");
        return m_vrow.back()->mapper().firstIndex() +
	  m_vrow.back()->mapper().freeSize();
    }

    /// Returns the number of blocks in the matrix, corresponding to
    /// variables/components
    index_t numBlocks() const
    {
        index_t nb = 0;
        for (size_t i = 0; i!=m_vrow.size(); ++i)
            nb += m_vrow[i]->dim();
        return nb;
    }

    /// Returns a reference to the options structure
    gsOptionList & options() {return m_options;}

    /// @brief Returns the left-hand global matrix
    const gsSparseMatrix<T> & matrix() const { return m_matrix; }

    /// @brief Writes the resulting matrix in \a out. The internal matrix is moved.
    void matrix_into(gsSparseMatrix<T> & out) { out = give(m_matrix); }

    EIGEN_STRONG_INLINE gsSparseMatrix<T> giveMatrix()
    {
         gsSparseMatrix<T> rvo;
         rvo.swap(m_matrix);
          return rvo;
    }

    /// @brief Returns the right-hand side vector(s)
    const gsMatrix<T> & rhs() const { return m_rhs; }

    /// @brief Writes the resulting vector in \a out. The internal data is moved.
    void rhs_into(gsMatrix<T> & out) { out = give(m_rhs); }

    /// \brief Sets the domain of integration.
    /// \warning Must be called before any computation is requested
    void setIntegrationElements(const gsMultiBasis<T> & mesh)
    { m_exprdata->setMultiBasis(mesh); }

#if EIGEN_HAS_RVALUE_REFERENCES
    void setIntegrationElements(const gsMultiBasis<T> &&) = delete;
    //const gsMultiBasis<T> * c++98
#endif

    /// \brief Returns the domain of integration
    const gsMultiBasis<T> & integrationElements() const
    { return m_exprdata->multiBasis(); }

    const typename gsExprHelper<T>::Ptr exprData() const { return m_exprdata; }

    /// Registers \a mp as an isogeometric geometry map and return a handle to it
    geometryMap getMap(const gsMultiPatch<T> & mp) //conv->tmp->error
    { return m_exprdata->getMap(mp); }

    /// Registers \a g as an isogeometric geometry map and return a handle to it
    geometryMap getMap(const gsFunction<T> & g)
    { return m_exprdata->getMap(g); }

    /// Registers \a mp as an isogeometric (both trial and test) space
    /// and return a handle to it
    space getSpace(const gsFunctionSet<T> & mp, index_t dim = 1, index_t id = 0)
    {
        //if multiBasisSet() then check domainDom
        GISMO_ASSERT(1==mp.targetDim(), "Expecting scalar source space");
        GISMO_ASSERT(static_cast<size_t>(id)<m_vrow.size(),
                     "Given ID "<<id<<" exceeds "<<m_vrow.size()-1 );
        expr::gsFeSpace<T> & u = m_exprdata->getSpace(mp,dim);
        u.setId(id);
        m_vrow[id] = m_vcol[id] = &u;
        return u;
    }

    /// \brief Registers \a mp as an isogeometric test space
    /// corresponding to trial space \a u and return a handle to it
    ///
    /// \note Both test and trial spaces are registered at once by
    /// gsExprAssembler::getSpace.
    ///
    ///Use this function after calling gsExprAssembler::getSpace when
    /// a distinct test space is requred (eg. Petrov-Galerkin
    /// methods).
    ///
    /// \note The dimension is set to the same as \a u, unless the caller
    /// sets as a third argument a new value.
    space getTestSpace(space u, const gsFunctionSet<T> & mp, index_t dim = -1)
    {
        //GISMO_ASSERT(0!=u.mapper(), "Not a space"); // done on initSystem
        expr::gsFeSpace<T> & s = m_exprdata->getSpace(mp,(-1 == dim ? u.dim() : dim));
        space uu = static_cast<space>(u);
        s.setId(uu.id());
        m_vrow[s.id()] = &s;
        return s;
    }

    /// Return the variable (previously created by getSpace) with the given \a id
    space trialSpace(const index_t id) const
    {
        GISMO_ASSERT(NULL!=m_vcol[id], "Not set.");
        return *m_vcol[id];
    }

    /// Return the trial space of a pre-existing test space \a v
    space trialSpace(space v) const { return trialSpace(v.id()); }

    /// Return the variable (previously created by getTrialSpace) with the given \a id
    space testSpace(const index_t id)
    {
        GISMO_ASSERT(NULL!=m_vrow[id], "Not set.");
        return *m_vrow[id];
    }

    /// Return the test space of a pre-existing trial space \a u
    space testSpace(space u) const { return testSpace(u.id()); }

    /// Registers \a func as a variable and returns a handle to it
    ///
    variable getCoeff(const gsFunctionSet<T> & func)
    { return m_exprdata->getVar(func, 1); }

    /// Registers \a func as a variable defined on \a G and returns a handle to it
    ///
    variable getCoeff(const gsFunctionSet<T> & func, geometryMap G)
    { return m_exprdata->getVar(func,G); }

    /// \brief Registers a representation of a solution variable from
    /// space \a s, based on the vector \a cf.
    ///
    /// The vector \a cf should have the structure of the columns of
    /// the system matrix this->matrix(). The returned handle
    /// corresponds to a function in the space \a s
    solution getSolution(space s, gsMatrix<T> & cf) const
    {
        // todo: if (m_exprdata->isSpace(u));
        //space s = static_cast<space>(u);
        return solution(s, cf);
    }

    variable getBdrFunction() const { return m_exprdata->getMutVar(); }

    element getElement() const { return m_element; }

    void computeDirichletDofs2(short_t unk);
    void computeDirichletDofsIntpl2(const expr::gsFeSpace<T> & u);
    void computeDirichletDofsL2Proj(const expr::gsFeSpace<T> & u);
    void setFixedDofVector(gsMatrix<T> & dof, short_t unk = 0);
    void setFixedDofs(const gsMatrix<T> & coefMatrix, short_t unk = 0, size_t patch = 0);

    /// \brief Initializes the sparse system (sparse matrix and rhs)
    void initSystem(bool resetFirst = true)
    {
        // Check spaces.nPatches==mesh.patches
        initMatrix(resetFirst);
        m_rhs.setZero(numDofs(), 1);

        // USE setup and gsDirichletValues instead
        //for (size_t i = 0; i!= m_vcol.size(); ++i)
        // computeDirichletDofs2(i);
    }

    /// \brief Initializes the sparse matrix only
    void initMatrix(bool resetFirst = true)
    {
        if (resetFirst)
            resetSpaces();
        resetDimensions();
        m_matrix = gsSparseMatrix<T>(numTestDofs(), numDofs());

        if ( 0 == m_matrix.rows() || 0 == m_matrix.cols() )
            gsWarn << " No internal DOFs, zero sized system.\n";
        else
        {
            // Pick up values from options
            const T bdA       = m_options.getReal("bdA");
            const index_t bdB = m_options.getInt("bdB");
            const T bdO       = m_options.getReal("bdO");
            T nz = 1;
            const short_t dim = m_exprdata->multiBasis().domainDim();
            for (short_t i = 0; i != dim; ++i)
                nz *= bdA * m_exprdata->multiBasis().maxDegree(i) + bdB;

            m_matrix.reservePerColumn(numBlocks()*cast<T,index_t>(nz*(1.0+bdO)) );
        }
    }

    /// \brief Initializes the right-hand side vector only
    void initVector(const index_t numRhs = 1, bool resetFirst = true)
    {
        if (resetFirst)
            resetSpaces();
        resetDimensions();
        m_rhs.setZero(numDofs(), numRhs);
    }

    /// Returns a block view of the system matrix, each block
    /// corresponding to a different space, or to different groups of
    /// dofs, in case of calar problems
    matBlockView matrixBlockView()
    {
        GISMO_ASSERT( m_vcol.back()->mapper().isFinalized(),
                      "initSystem() has not been called.");
        gsVector<index_t> rowSizes, colSizes;
        _blockDims(rowSizes, colSizes);
        return m_matrix.blockView(rowSizes,colSizes);
    }

    /// Returns a const block view of the system matrix, each block
    /// corresponding to a different space, or to different groups of
    /// dofs, in case of calar problems
    matConstBlockView matrixBlockView() const
    {
        GISMO_ASSERT( m_vcol.back()->mapper().isFinalized(),
                      "initSystem() has not been called.");
        gsVector<index_t> rowSizes, colSizes;
        _blockDims(rowSizes, colSizes);
        return m_matrix.blockView(rowSizes,colSizes);
    }

    /// Set the assembler options
    void setOptions(gsOptionList opt) { m_options = opt; } // gsOptionList opt
    // .swap(opt) todo

#   if(__cplusplus >= 201103L || _MSC_VER >= 1600 || defined(__DOXYGEN__))
    /// Adds the expressions \a args to the system matrix/rhs
    ///
    /// The arguments are considered as integrals over the whole domain
    /// \sa gsExprAssembler::setIntegrationElements
    template<class... expr> void assemble(expr... args);

    /// Adds the expressions \a args to the system matrix/rhs
    ///
    /// The arguments are considered as integrals over the boundary parts in \a BCs
    template<class... expr> void assemble(const bcRefList & BCs, expr... args);

    /*
      template<class... expr> void assemble(const ifContainer & iFaces, expr... args);
      template<class... expr> void collocate(expr... args);// eg. collocate(-ilapl(u), f)
    */
#else
    template<class E1> void assemble(const expr::_expr<E1> & a1)
    {assemble(a1,nullExpr(),nullExpr(),nullExpr(),nullExpr());}
    template <class E1, class E2>
    void assemble(const expr::_expr<E1> & a1, const expr::_expr<E2> & a2)
    {assemble(a1,a2,nullExpr(),nullExpr(),nullExpr());}
    template <class E1, class E2, class E3>
    void assemble(const expr::_expr<E1> & a1, const expr::_expr<E2> & a2,
                  const expr::_expr<E3> & a3)
    {assemble(a1,a2,a3,nullExpr(),nullExpr());}
    template <class E1, class E2, class E3, class E4, class E5>
    void assemble(const expr::_expr<E1> & a1, const expr::_expr<E2> & a2,
                  const expr::_expr<E3> & a3, const expr::_expr<E4> & a4)
    {assemble(a1,a2,a3,a4,nullExpr());}
    template <class E1, class E2, class E3, class E4, class E5>
    void assemble(const expr::_expr<E1> & a1, const expr::_expr<E2> & a2,
                  const expr::_expr<E3> & a3, const expr::_expr<E4> & a4,
                  const expr::_expr<E5> & a5 );

    template<class E1> void assemble(const bcRefList & BCs, const expr::_expr<E1> & a1);
#   endif

    template<class E1, class E2>
    void assembleLhsRhsBc(const expr::_expr<E1> & exprLhs,
                          const expr::_expr<E2> & exprRhs,
                          const bcContainer & BCs)
    {
        space rvar = static_cast<space>(exprLhs.rowVar());
        GISMO_ASSERT(m_exprdata->exists(rvar), "Error - inexistent variable.");
        space cvar = static_cast<space>(exprLhs.colVar());
        GISMO_ASSERT(m_exprdata->exists(cvar), "Error - inexistent variable.");
        GISMO_ASSERT(&rvar==&exprRhs.rowVar(), "Inconsistent left and right hand side");
        assembleLhsRhsBc_impl<true,true>(exprLhs, exprRhs, rvar, cvar, BCs);
    }

    template<class E1>
    void assembleRhsBc(const expr::_expr<E1> & exprRhs, const bcContainer & BCs)
    {
        space var = static_cast<space>(exprRhs.rowVar());
        GISMO_ASSERT(m_exprdata->exists(var), "Error - inexistent variable.");
        assembleLhsRhsBc_impl<false,true>(nullExpr(), exprRhs, var, var, BCs);
    }

    template<class E1>
    void assembleInterface(const expr::_expr<E1> & exprInt)
    {
        space rvar = static_cast<space>(exprInt.rowVar());
        space cvar = static_cast<space>(exprInt.colVar());
        assembleInterface_impl<true,false>(exprInt, nullExpr(), rvar, rvar, m_exprdata->multiBasis().topology().interfaces() );
    }

    template<class E1>
    void assembleRhsInterface(const expr::_expr<E1> & exprInt, const ifContainer & iFaces)
    {
        space rvar = static_cast<space>(exprInt.rowVar());
        GISMO_ASSERT(m_exprdata->exists(rvar), "Error - inexistent variable.");
        assembleInterface_impl<false,true>(nullExpr(), exprInt, rvar, rvar, iFaces);
    }

private:

    void _blockDims(gsVector<index_t> & rowSizes,
                    gsVector<index_t> & colSizes)
    {
        if (1==m_vcol.size() && 1==m_vrow.size())
        {
            const gsDofMapper & dm = m_vcol.back()->mapper();
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
                rowSizes[r] = m_vrow[r]->dim() * m_vrow[r]->mapper().freeSize();
            colSizes.resize(m_vcol.size());
            for (index_t c = 0; c != colSizes.size(); ++c) // for all col-blocks
                colSizes[c] = m_vcol[c]->dim() * m_vcol[c]->mapper().freeSize();
        }
    }

    /// \brief Reset the dimensions of all involved spaces.
    /// Called internally by the init* functions
    void resetDimensions();

    void resetSpaces();

    // template<bool left, bool right, class E1, class E2>
    // void assembleLhsRhs_impl(const expr::_expr<E1> & exprLhs,
    //                          const expr::_expr<E2> & exprRhs,
    //                          space rvar, space cvar);

    template<bool left, bool right, class E1, class E2>
    void assembleLhsRhsBc_impl(const expr::_expr<E1> & exprLhs,
                               const expr::_expr<E2> & exprRhs,
                               space rvar, space cvar,
                               const bcContainer & BCs);

    template<bool left, bool right, class E1, class E2>
    void assembleInterface_impl(const expr::_expr<E1> & exprLhs,
                                const expr::_expr<E2> & exprRhs,
                                space rvar, space cvar,
                                const ifContainer & iFaces);

#if __cplusplus >= 201103L || _MSC_VER >= 1600 // c++11
    template <class op, class E1>
    void _apply(op _op, const expr::_expr<E1> & firstArg) {_op(firstArg);}
    template <class op, class E1, class... Rest>
    void _apply(op _op, const expr::_expr<E1> & firstArg, Rest... restArgs)
    { _op(firstArg); _apply<op>(_op, restArgs...); }
#endif

    struct __setFlag
    {
        template <typename E> void operator() (const gismo::expr::_expr<E> & v)
        { v.setFlag(); }

        void operator() (const expr::_expr<expr::gsNullExpr<T> > &) {}
    } _setFlag;

    struct __printExpr
    {
        template <typename E> void operator() (const gismo::expr::_expr<E> & v)
        { v.print(gsInfo);gsInfo<<"\n"; }
    } _printExpr;

    struct _eval
    {
        gsSparseMatrix<T> & m_matrix;
        gsMatrix<T>       & m_rhs;
        const gsVector<T> & m_quWeights;
        index_t       m_patchInd;
        gsMatrix<T>         localMat;

        _eval(gsSparseMatrix<T> & _matrix,
              gsMatrix<T>       & _rhs,
              const gsVector<>  & _quWeights)
        : m_matrix(_matrix), m_rhs(_rhs),
          m_quWeights(_quWeights), m_patchInd(0)
        { }

        void setPatch(const index_t p) { m_patchInd=p; }

        template <typename E> void operator() (const gismo::expr::_expr<E> & ee)
        {
            // ------- Compute  -------
            const T * w = m_quWeights.data();
            localMat.noalias() = (*w) * ee.eval(0);
            for (index_t k = 1; k != m_quWeights.rows(); ++k)
                localMat.noalias() += (*(++w)) * ee.eval(k);

            //  ------- Accumulate  -------
            if (E::isMatrix())
                push<true>(ee.rowVar(), ee.colVar(), m_patchInd);
            else if (E::isVector())
                push<false>(ee.rowVar(), ee.colVar(), m_patchInd);
            else
            {
                GISMO_ERROR("Something went wrong at this point (rowspan: "<< E::rowSpan<< ", colSpan: "<< E::colSpan <<")");
                //GISMO_ASSERTrowSpan() && (!colSpan())
            }

        }// operator()

        void operator() (const expr::_expr<expr::gsNullExpr<T> > &) {}

        template<bool isMatrix> void push(const expr::gsFeSpace<T> & v,
                                          const expr::gsFeSpace<T> & u,
                                          const index_t patchInd)
        {
            GISMO_ASSERT(v.isValid(), "The row space is not valid");
            GISMO_ASSERT(!isMatrix || u.isValid(), "The column space is not valid");
            GISMO_ASSERT(isMatrix || (0!=m_rhs.size()), "The right-hand side vector is not initialized");

            const index_t cd            = u.dim();
            const index_t rd            = v.dim();
            const gsDofMapper  & colMap = u.mapper();
            const gsDofMapper  & rowMap = v.mapper();
            gsMatrix<index_t> & colInd0 = const_cast<gsMatrix<index_t>&>(u.data().actives);
            gsMatrix<index_t> & rowInd0 = const_cast<gsMatrix<index_t>&>(v.data().actives);
            const gsMatrix<T>  & fixedDofs = u.fixedPart();

            gsMatrix<index_t> rowInd, colInd;
            rowMap.localToGlobal(rowInd0, patchInd, rowInd);

            if (isMatrix)
            {
                //if (&rowInd0!=&colInd0)
                colMap.localToGlobal(colInd0, patchInd, colInd);
                GISMO_ASSERT( colMap.boundarySize()==fixedDofs.size(),
                              "Invalid values for fixed part");
            }
            for (index_t r = 0; r != rd; ++r)
            {
                const index_t rls = r * rowInd0.rows();     //local stride
                for (index_t i = 0; i != rowInd0.rows(); ++i)
                {
                    const index_t ii = rowMap.index(rowInd0.at(i),patchInd,r); // N_i
                    if ( rowMap.is_free_index(ii) )
                    {
                        for (index_t c = 0; c != cd; ++c)
                        {
                            if (isMatrix)
                            {
                                const index_t cls = c * colInd0.rows();     //local stride

                                for (index_t j = 0; j != colInd0.rows(); ++j)
                                {
                                    if ( 0 == localMat(rls+i,cls+j) ) continue;

                                    const index_t jj = colMap.index(colInd0.at(j),patchInd,c); // N_j
                                    if ( colMap.is_free_index(jj) )
                                    {
                                        // If matrix is symmetric, we could
                                        // store only lower triangular part
                                        //if ( (!symm) || jj <= ii )
                                        m_matrix.coeffRef(ii, jj) += localMat(rls+i,cls+j);
                                    }
                                    else // colMap.is_boundary_index(jj) )
                                    {
                                        // Symmetric treatment of eliminated BCs
                                        // GISMO_ASSERT(1==m_rhs.cols(), "-");
                                        m_rhs.at(ii) -= localMat(rls+i,cls+j) *
                                            fixedDofs.at(colMap.global_to_bindex(jj));
                                    }
                                }
                            }
                            else
                            {
                                m_rhs.row(ii) += localMat.row(rls+i);
                            }
                        }
                    }
                }
            }
        }//push

    };

}; // gsExprAssembler

template<class T>
gsOptionList gsExprAssembler<T>::defaultOptions()
{
    gsOptionList opt;
    opt.addInt("DirichletValues"  , "Method for computation of Dirichlet DoF values [100..103]", 101);
    opt.addReal("quA", "Number of quadrature points: quA*deg + quB", 1.0  );
    opt.addInt ("quB", "Number of quadrature points: quA*deg + quB", 1    );
    opt.addReal("bdA", "Estimated nonzeros per column of the matrix: bdA*deg + bdB", 2.0  );
    opt.addInt ("bdB", "Estimated nonzeros per column of the matrix: bdA*deg + bdB", 1    );
    opt.addReal("bdO", "Overhead of sparse mem. allocation: (1+bdO)(bdA*deg + bdB) [0..1]", 0.333);
    return opt;
}

template<class T>
void gsExprAssembler<T>::computeDirichletDofs2(short_t unk)
{
    expr::gsFeSpace<T> & u = *m_vcol[unk];

    //if ( m_options.getInt("DirichletStrategy") == dirichlet::nitsche)
    //    return; // Nothing to compute

    //const gsMultiBasis<T> & mbasis = dynamic_cast<const gsMultiBasis<T>&>(u.source());

    // eg. not penalize
    const gsDofMapper & mapper = u.mapper();

    switch ( m_options.getInt("DirichletValues") )
    {
    case dirichlet::homogeneous:
        // If we have a homogeneous Dirichlet problem fill boundary
        // DoFs with zeros
        u.fixedPart().setZero(mapper.boundarySize(), 1 );
        break;
    case dirichlet::interpolation:
        computeDirichletDofsIntpl2(u);
        break;
    case dirichlet::l2Projection:
        computeDirichletDofsL2Proj(u); //this->computeDirichletDofsL2Proj(mapper, mbasis,unk);
        break;
    case dirichlet::user :
        // Assuming that the DoFs are already set by the user
        GISMO_ENSURE( u.fixedPart().size() == mapper.boundarySize(),
                      "The Dirichlet DoFs are not set.");
        break;
    default:
        GISMO_ERROR("Something went wrong with Dirichlet values.");
    }

    /* Corner values -- todo
    for ( typename gsBoundaryConditions<T>::const_citerator
              it = bbc.cornerBegin();
          it != bbc.cornerEnd(); ++it )
    {
        if(it->unknown == unk)
        {
            const index_t i  = mbasis[it->patch].functionAtCorner(it->corner);
            const index_t ii = mapper.bindex( i , it->patch );
            u.fixedPart().row(ii).setConstant(it->value);
        }
        else
            continue;
    }
    */
}

template<class T>
void gsExprAssembler<T>::setFixedDofVector(gsMatrix<T> & vals, short_t unk)
{
    expr::gsFeSpace<T> & u = *m_vcol[unk];
    gsMatrix<T>        & fixedDofs = const_cast<expr::gsFeSpace<T>&>(u).fixedPart();
    fixedDofs.swap(vals);
    vals.resize(0, 0);
    // Assuming that the DoFs are already set by the user
    GISMO_ENSURE( fixedDofs.size() == u.mapper().boundarySize(),
                     "The Dirichlet DoFs were not provided correctly.");
}

template<class T>
void gsExprAssembler<T>::setFixedDofs(const gsMatrix<T> & coefMatrix, short_t unk, size_t patch)
{
    GISMO_ASSERT( m_options.getInt("DirichletValues") == dirichlet::user, "Incorrect options");

    expr::gsFeSpace<T> & u = *m_vcol[unk];
    //const index_t dirStr = m_options.getInt("DirichletStrategy");
    const gsMultiBasis<T> & mbasis = *dynamic_cast<const gsMultiBasis<T>* >(&(u).source());

    //const gsBoundaryConditions<> & bbc = u.hasBc() ? u.bc() : gsBoundaryConditions<>();

    const gsDofMapper & mapper = u.mapper();
//    const gsDofMapper & mapper =
//        dirichlet::elimination == dirStr ? u.mapper()
//        : mbasis.getMapper(dirichlet::elimination,
//                           static_cast<iFace::strategy>(m_options.getInt("InterfaceStrategy")),
//                           bbc, u.id()) ;

    gsMatrix<T> & fixedDofs = const_cast<expr::gsFeSpace<T>& >(u).fixedPart();
    GISMO_ASSERT(fixedDofs.size() == mapper.boundarySize(),
                 "Fixed DoFs were not initialized.");

    // for every side with a Dirichlet BC
    // for ( typename gsBoundaryConditions<T>::const_iterator
    //       it =  bbc.dirichletBegin();
    //       it != bbc.dirichletEnd()  ; ++it )
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

template<class T> void gsExprAssembler<T>::resetSpaces()
{
    for (size_t i = 0; i!=m_vcol.size(); ++i)
    {
        GISMO_ASSERT(NULL!=m_vcol[i], "The assembler spaces where not set.");
        m_vcol[i]->reset();

        if ( m_vcol[i] != m_vrow[i] )
        {
            GISMO_ASSERT(NULL!=m_vrow[i], "The assembler spaces where not set.");
            m_vrow[i]->reset();
        }
    }
}

template<class T> void gsExprAssembler<T>::resetDimensions()
{
    for (size_t i = 1; i!=m_vcol.size(); ++i)
    {
        m_vcol[i]->mapper().setShift(m_vcol[i-1]->mapper().firstIndex() +
                                     m_vcol[i-1]->dim()*m_vcol[i-1]->mapper().freeSize() );

        if ( m_vcol[i] != m_vrow[i] )
            m_vrow[i]->mapper().setShift(m_vrow[i-1]->mapper().firstIndex() +
                                         m_vrow[i-1]->dim()*m_vrow[i-1]->mapper().freeSize() );
    }
}

template<class T>
#if(__cplusplus >= 201103L || _MSC_VER >= 1600 || defined(__DOXYGEN__)) // c++11
template<class... expr>
void gsExprAssembler<T>::assemble(expr... args)
#else
    template <class E1, class E2, class E3, class E4, class E5>
    void gsExprAssembler<T>::assemble( const expr::_expr<E1> & a1, const expr::_expr<E2> & a2,
    const expr::_expr<E3> & a3, const expr::_expr<E4> & a4, const expr::_expr<E5> & a5)
#endif
{
    GISMO_ASSERT(matrix().cols()==numDofs(), "System not initialized");

    // initialize flags
    m_exprdata->initFlags(SAME_ELEMENT|NEED_ACTIVE, SAME_ELEMENT);
#   if __cplusplus >= 201103L || _MSC_VER >= 1600
    _apply(_setFlag, args...);
    //_apply(_printExpr, args...);
#   else
    _setFlag(a1);_setFlag(a1);_setFlag(a2);_setFlag(a4);_setFlag(a5);
#   endif
    gsQuadRule<T> QuRule;  // Quadrature rule
    gsVector<T> quWeights; // quadrature weights

    _eval ee(m_matrix, m_rhs, quWeights);

    for (unsigned patchInd = 0; patchInd < m_exprdata->multiBasis().nBases(); ++patchInd)
    {
        ee.setPatch(patchInd);
        QuRule = gsQuadrature::get(m_exprdata->multiBasis().basis(patchInd), m_options);

        // Initialize domain element iterator for current patch
        typename gsBasis<T>::domainIter domIt =  // add patchInd to domainiter ?
            m_exprdata->multiBasis().basis(patchInd).makeDomainIterator();
        m_element.set(*domIt);

        // Start iteration over elements of patchInd
        for (; domIt->good(); domIt->next() )
        {
            // Map the Quadrature rule to the element
            QuRule.mapTo( domIt->lowerCorner(), domIt->upperCorner(),
                          m_exprdata->points(), quWeights);

            // Perform required pre-computations on the quadrature nodes
            m_exprdata->precompute(patchInd);
            //m_exprdata->precompute(QuRule, *domIt); // todo

            // Assemble contributions of the element
#           if __cplusplus >= 201103L || _MSC_VER >= 1600
            _apply(ee, args...);
#           else
            ee(a1);ee(a2);ee(a3);ee(a4);ee(a5);
#           endif
        }
    }

    m_matrix.makeCompressed();
}

template<class T>
#if __cplusplus >= 201103L || _MSC_VER >= 1600 // c++11
template<class... expr>
void gsExprAssembler<T>::assemble(const bcRefList & BCs, expr... args)
#else
template <class E1>
void gsExprAssembler<T>::assemble(const bcRefList & BCs, const expr::_expr<E1> & a1)
#endif
{
    // initialize flags
    m_exprdata->initFlags(SAME_ELEMENT|NEED_ACTIVE, SAME_ELEMENT);
#   if __cplusplus >= 201103L || _MSC_VER >= 1600
    _apply(_setFlag, args...);
#   else
    _setFlag(a1);
#   endif

    gsVector<T> quWeights;// quadrature weights
    gsQuadRule<T>  QuRule;

    _eval ee(m_matrix, m_rhs, quWeights);

    for (typename bcRefList::const_iterator iit = BCs.begin(); iit!= BCs.end(); ++iit)
    {
        const boundary_condition<T> * it = &iit->get();

        QuRule = gsQuadrature::get(m_exprdata->multiBasis().basis(it->patch()), m_options, it->side().direction());

        m_exprdata->mapData.side = it->side();

        // Update boundary function source
        m_exprdata->setMutSource(*it->function(), it->parametric());
        //mutVar.registerVariable(func, mutData);

        typename gsBasis<T>::domainIter domIt =
            m_exprdata->multiBasis().basis(it->patch()).makeDomainIterator(it->side());
        m_element.set(*domIt);

        // Start iteration over elements
        for (; domIt->good(); domIt->next() )
        {
            // Map the Quadrature rule to the element
            QuRule.mapTo( domIt->lowerCorner(), domIt->upperCorner(),
                          m_exprdata->points(), quWeights);

            // Perform required pre-computations on the quadrature nodes
            m_exprdata->precompute(it->patch());

            // Assemble contributions of the element
#           if __cplusplus >= 201103L || _MSC_VER >= 1600
            _apply(ee, args...);
#           else
            ee(a1);
#           endif
        }
    }

    //this->finalize();
    m_matrix.makeCompressed();
    //g_bd.clear();
    //mutVar.clear();
}


template<class T>
template<bool left, bool right, class E1, class E2>
void gsExprAssembler<T>::assembleLhsRhsBc_impl(const expr::_expr<E1> & exprLhs,
                                               const expr::_expr<E2> & exprRhs,
                                               space rvar, space cvar,
                                               const bcContainer & BCs)
{
    //GISMO_ASSERT( exprRhs.isVector(), "Expecting vector expression");

    // initialize flags
    m_exprdata->initFlags(SAME_ELEMENT|NEED_ACTIVE, SAME_ELEMENT);
    if (left ) exprLhs.setFlag();
    if (right) exprRhs.setFlag();

    gsVector<T> quWeights;// quadrature weights
    gsQuadRule<T>  QuRule;
    _eval ee(m_matrix, m_rhs, quWeights);

    for (typename bcContainer::const_iterator it = BCs.begin(); it!= BCs.end(); ++it)
    {
        QuRule = gsQuadrature::get(m_exprdata->multiBasis().basis(it->patch()), m_options, it->side().direction());

        m_exprdata->mapData.side = it->side();

        // Update boundary function source
        m_exprdata->setMutSource(*it->function(), it->parametric());
        //mutVar.registerVariable(func, mutData);

        typename gsBasis<T>::domainIter domIt =
            m_exprdata->multiBasis().basis(it->patch()).makeDomainIterator(it->side());
        m_element.set(*domIt);

        // Start iteration over elements
        for (; domIt->good(); domIt->next() )
        {
            // Map the Quadrature rule to the element
            QuRule.mapTo( domIt->lowerCorner(), domIt->upperCorner(),
                          m_exprdata->points(), quWeights);

            // Perform required pre-computations on the quadrature nodes
            m_exprdata->precompute(it->patch());

            ee.setPatch(it->patch());
    	    ee(exprLhs);
	        ee(exprRhs);
        }
    }

    //this->finalize();
    m_matrix.makeCompressed();
    //g_bd.clear();
    //mutVar.clear();
}

template<class T>
template<bool left, bool right, class E1, class E2>
void gsExprAssembler<T>::assembleInterface_impl(const expr::_expr<E1> & exprLhs,
                                                const expr::_expr<E2> & exprRhs,
                                                space rvar, space cvar,
                                                const ifContainer & iFaces)
{
    //GISMO_ASSERT( exprRhs.isVector(), "Expecting vector expression");

    // initialize flags

    m_exprdata->initFlags(SAME_ELEMENT|NEED_ACTIVE, SAME_ELEMENT);
    if (left ) exprLhs.setFlag();
    if (right) exprRhs.setFlag();
    //m_exprdata->parse(exprLhs,exprRhs);
    //m_exprdata->parse(exprRhs);

    gsVector<T> quWeights;// quadrature weights
    gsQuadRule<T>  QuRule;
    _eval ee(m_matrix, m_rhs, quWeights);

    //gsMatrix<T> tmp;

    for (gsBoxTopology::const_iiterator it = iFaces.begin();
         it != iFaces.end(); ++it )
    {
        const boundaryInterface & iFace = *it;
        const index_t patch1 = iFace.first() .patch;
        //const index_t patch2 = iFace.second().patch;
        //const gsAffineFunction<T> interfaceMap(m_pde_ptr->patches().getMapForInterface(bi));

        QuRule = gsQuadrature::get(m_exprdata->multiBasis().basis(patch1),
                                   m_options, iFace.first().side().direction());

        m_exprdata->mapData.side = iFace.first().side(); // (!)

        typename gsBasis<T>::domainIter domIt =
            m_exprdata->multiBasis().basis(patch1).makeDomainIterator(iFace.first().side());
        m_element.set(*domIt);

        // Start iteration over elements
        for (; domIt->good(); domIt->next() )
        {
            // Map the Quadrature rule to the element
            QuRule.mapTo( domIt->lowerCorner(), domIt->upperCorner(),
                          m_exprdata->points(), quWeights);

            // Perform required pre-computations on the quadrature nodes
            m_exprdata->precompute(patch1);

            // DG: need data1, data2
            // coupling: need to know patch1/patch2

//            interfaceMap.eval_into(m_exprdata->points(), tmp);
//            m_exprdata->points().swap(tmp);
//            m_exprdata->precompute(patch2);

        ee.setPatch(patch1);
	    ee(exprLhs);
	    ee(exprRhs);
        }
    }

    m_matrix.makeCompressed();
}


template<class T> //
void gsExprAssembler<T>::computeDirichletDofsIntpl2(const expr::gsFeSpace<T> & u)
{
    const gsDofMapper  & mapper    = u.mapper();
    gsMatrix<T>        & fixedDofs = const_cast<expr::gsFeSpace<T>&>(u).fixedPart();
    fixedDofs.resize(mapper.boundarySize(), u.dim() );
    const index_t parDim = u.source().domainDim();

    const gsMultiBasis<T> & mbasis =
        *dynamic_cast<const gsMultiBasis<T>*>(&u.source());

    // Iterate over all patch-sides with Boundary conditions
    typedef typename gsBoundaryConditions<T>::bcRefList bcRefList;
    for ( typename bcRefList::const_iterator iit =  u.bc().begin();
          iit != u.bc().end()  ; ++iit )
    {
        const boundary_condition<T> * it = &iit->get();
        const index_t com = it->unkComponent();

        const index_t k = it->patch();
        if( it->unknown()!=u.id() )
            continue;
        const gsBasis<T> & basis = mbasis[k];

        // Get dofs on this boundary
        const gsMatrix<index_t> boundary = basis.boundary(it->side());

        // If the condition is homogeneous then fill with zeros
        if ( it->isHomogeneous() )
        {
            for (index_t i=0; i!= boundary.size(); ++i)
            {
                const index_t ii= mapper.bindex( boundary.at(i) , k, com );
                fixedDofs.at(ii) = 0;
            }
            continue;
        }

        // Get the side information
        short_t dir = it->side().direction( );
        index_t param = (it->side().parameter() ? 1 : 0);

        // Compute grid of points on the face ("face anchors")
        std::vector< gsVector<T> > rr;
        rr.reserve( parDim );

        for ( short_t i=0; i < parDim; ++i)
        {
            if ( i==dir )
            {
                gsVector<T> b(1);
                b[0] = ( basis.component(i).support() ) (0, param);
                rr.push_back(b);
            }
            else
            {
                rr.push_back( basis.component(i).anchors().transpose() );
            }
        }

        // GISMO_ASSERT(it->function()->targetDim() == u.dim(),
        //              "Given Dirichlet boundary function does not match problem dimension."
        //              <<it->function()->targetDim()<<" != "<<u.dim()<<"\n");

        // Compute dirichlet values
        gsMatrix<T> fpts;
        if ( it->parametric() )
            fpts = it->function()->eval( gsPointGrid<T>( rr ) );
        else
            fpts = it->function()->eval(
                m_exprdata->getMap().source().piece(it->patch()).eval(  gsPointGrid<T>( rr ) ) );

        /*
        if ( fpts.rows() != u.dim() )
        {
            // assume scalar
            gsMatrix<T> tmp(u.dim(), fpts.cols());
            tmp.setZero();
            gsDebugVar(!dir);
            tmp.row(!dir) = (param ? 1 : -1) * fpts; // normal !
            fpts.swap(tmp);
        }
        */

        // Interpolate dirichlet boundary
        typename gsBasis<T>::uPtr h = basis.boundaryBasis(it->side());
        typename gsGeometry<T>::uPtr geo = h->interpolateAtAnchors(fpts);
        const gsMatrix<T> & dVals =  geo->coefs();

        // Save corresponding boundary dofs
        for (index_t l=0; l!= boundary.size(); ++l)
        {
            const index_t ii = mapper.bindex( boundary.at(l) , it->patch(), com);
            fixedDofs.at(ii) = dVals.at(l);
        }
    }
}


/*
template<class T> //
void gsExprAssembler<T>::computeDirichletDofsIntpl3(const expr::gsFeSpace<T> & u)
{
    const gsDofMapper  & mapper    = u.mapper();
    gsMatrix<T>        & fixedDofs = const_cast<expr::gsFeSpace<T>&>(u).fixedPart();
    const index_t bsz = mapper.boundarySize();
    fixedDofs.resize(bsz, u.dim() );
    gsMatrix<T> pt, val, rhs;
    rhs.resize(bsz, u.dim() );
    gsMatrix<index_t> act;
    gsSparseMatrix<T> cmat(bsz, bsz);
    // todo: reserve

    const gsMultiBasis<T> & mbasis =
        *dynamic_cast<const gsMultiBasis<T>*>(&u.source());

    // Iterate over all patch-sides with Boundary conditions
    // Iterate over all patch-sides with Dirichlet-boundary conditions
    typedef typename gsBoundaryConditions<T>::bcRefList bcRefList;
    for ( typename bcRefList::const_iterator iit =  u.bc().begin();
          iit != u.bc().end()  ; ++iit )
    {
        const boundary_condition<T> * it = &iit->get();

        GISMO_ASSERT(it->function()->targetDim() == u.dim(),
                     "Given Dirichlet boundary function does not match problem dimension."
                     <<it->function()->targetDim()<<" != "<<u.dim()<<"\n");

        const index_t k   = it->patch();
        if( it->unknown()!=u.id() )
            continue;
        const gsBasis<T> & basis = mbasis[k];

        // Get dofs on this boundary
        gsMatrix<index_t> boundary = basis.boundary(it->side());

        // If the condition is homogeneous then fill with zeros
        if ( it->isHomogeneous() )
        {
            for (index_t i=0; i!= boundary.size(); ++i)
            {
                const index_t ii= mapper.bindex( boundary.at(i) , k );
                fixedDofs.row(ii).setZero();
            }
            continue;
        }

        // Get anchor points for the respective dofs
        for ( index_t i=0; i != boundary.size(); ++i)
            // or: preimage ?
        {
            const index_t cc = mapper.bindex( boundary.at(l) , k );
            basis.anchor_into(boundary.at(i), pt);
            basis.active_into(pt, act);
            basis.eval_into  (pt, val);

            for ( index_t l = 0; l != act.size(); ++l)
            {
                const index_t ii = mapper.index( act.at(l) , k );
                if ( mapper.is_boundary_index(ii) ) // && cmat.isExpZero
                    cmat.insert(cc, mapper.global_to_bindex(ii)) = val.at(l);
            }

            // rhs
            if ( it->parametric() )
                it->function()->eval_into(pt, val);
            else
            {
                it->function()->eval_into(
                    m_exprdata->getMap().source().piece(it->patch()).eval(pt), val );
            }
            rhs.row(cc) = val.transpose();
        }

        // Interpolate dirichlet boundary
        // Solve overconstraint using QR ?
        //fixedDofs = ..
    }
}
//*/

template<class T>
void gsExprAssembler<T>::computeDirichletDofsL2Proj(const expr::gsFeSpace<T>& u)
{
    GISMO_ASSERT(&m_exprdata->getMap().source() != NULL, "Geometry not set, call setMap(...) first!");

    const gsDofMapper & mapper = u.mapper();
    gsMatrix<T> & fixedDofs = const_cast<expr::gsFeSpace<T>& >(u).fixedPart();
    fixedDofs.resize(mapper.boundarySize(), u.dim());

    const gsMultiBasis<T> & mbasis = *dynamic_cast<const gsMultiBasis<T>* >(&u.source());

    //const gsBoundaryConditions<> & bbc = u.hasBc() ? u.bc() : gsBoundaryConditions<>();

    // Set up matrix, right-hand-side and solution vector/matrix for
    // the L2-projection
    gsSparseEntries<T> projMatEntries;
    gsMatrix<T>        globProjRhs;
    globProjRhs.setZero(mapper.boundarySize(), u.dim());

    // Temporaries
    gsVector<T> quWeights;

    gsMatrix<T> rhsVals;
    gsMatrix<index_t> globIdxAct;
    gsMatrix<T> basisVals;

    gsMapData<T> md(NEED_MEASURE | SAME_ELEMENT);

    const gsMultiPatch<T> & mp = static_cast<const gsMultiPatch<T> &>(m_exprdata->getMap().source());

    // Iterate over all patch-sides with Dirichlet-boundary conditions
    typedef typename gsBoundaryConditions<T>::bcRefList bcRefList;
    for (typename bcRefList::const_iterator iit = u.bc().begin();
         iit != u.bc().end(); ++iit)
    {
        const boundary_condition<T> * iter = &iit->get();

        const short_t unk = iter->unknown();
        if(unk != u.id())
            continue;
        const index_t patchIdx   = iter->patch();
        const gsBasis<T> & basis = mbasis[patchIdx];

        const gsGeometry<T> & patch = mp.patch(patchIdx);

        // Set up quadrature to degree+1 Gauss points per direction,
        // all lying on iter->side() except from the direction which
        // is NOT along the element
        gsGaussRule<T> bdQuRule(basis, 1.0, 1, iter->side().direction());

        // Create the iterator along the given part boundary.
        typename gsBasis<T>::domainIter bdryIter = basis.makeDomainIterator(iter->side());

        for (; bdryIter->good(); bdryIter->next())
        {
            bdQuRule.mapTo(bdryIter->lowerCorner(), bdryIter->upperCorner(),
                           md.points, quWeights);

            patch.computeMap(md);

            // the values of the boundary condition are stored
            // to rhsVals. Here, "rhs" refers to the right-hand-side
            // of the L2-projection, not of the PDE.
            rhsVals = iter->function()->eval(m_exprdata->getMap().source().piece(patchIdx).eval(md.points));

            basis.eval_into(md.points, basisVals);

            // Indices involved here:
            // --- Local index:
            // Index of the basis function/DOF on the patch.
            // Does not take into account any boundary or interface conditions.
            // --- Global Index:
            // Each DOF has a unique global index that runs over all patches.
            // This global index includes a re-ordering such that all eliminated
            // DOFs come at the end.
            // The global index also takes care of glued interface, i.e., corresponding
            // DOFs on different patches will have the same global index, if they are
            // glued together.
            // --- Boundary Index (actually, it's a "Dirichlet Boundary Index"):
            // The eliminated DOFs, which come last in the global indexing,
            // have their own numbering starting from zero.

            // Get the global indices (second line) of the local
            // active basis (first line) functions/DOFs:
            basis.active_into(md.points.col(0), globIdxAct);
            mapper.localToGlobal(globIdxAct, patchIdx, globIdxAct);

            // Out of the active functions/DOFs on this element, collect all those
            // which correspond to a boundary DOF.
            // This is checked by calling mapper.is_boundary_index( global Index )

            // eltBdryFcts stores the row in basisVals/globIdxAct, i.e.,
            // something like a "element-wise index"
            std::vector<index_t> eltBdryFcts;
            eltBdryFcts.reserve(mapper.boundarySize());
            for (index_t i = 0; i < globIdxAct.rows(); i++)
            {
                if (mapper.is_boundary_index(globIdxAct(i, 0)))
                {
                    eltBdryFcts.push_back(i);
                }
            }

            // Do the actual assembly:
            for (index_t k = 0; k < md.points.cols(); k++)
            {
                const T weight_k = quWeights[k] * md.measure(k);

                // Only run through the active boundary functions on the element:
                for (size_t i0 = 0; i0 < eltBdryFcts.size(); i0++)
                {
                    // Each active boundary function/DOF in eltBdryFcts has...
                    // ...the above-mentioned "element-wise index"
                    const index_t i = eltBdryFcts[i0];
                    // ...the boundary index.
                    const index_t ii = mapper.global_to_bindex(globIdxAct(i));

                    for (size_t j0 = 0; j0 < eltBdryFcts.size(); j0++)
                    {
                        const index_t j = eltBdryFcts[j0];
                        const index_t jj = mapper.global_to_bindex(globIdxAct(j));

                        // Use the "element-wise index" to get the needed
                        // function value.
                        // Use the boundary index to put the value in the proper
                        // place in the global projection matrix.
                        projMatEntries.add(ii, jj, weight_k * basisVals(i, k) * basisVals(j, k));
                    } // for j

                    globProjRhs.row(ii) += weight_k * basisVals(i, k) * rhsVals.col(k).transpose();

                } // for i
            } // for k
        } // bdryIter
    } // boundaryConditions-Iterator

    gsSparseMatrix<T> globProjMat(mapper.boundarySize(), mapper.boundarySize());
    globProjMat.setFrom(projMatEntries);
    globProjMat.makeCompressed();

    // Solve the linear system:
    // The position in the solution vector already corresponds to the
    // numbering by the boundary index. Hence, we can simply take them
    // for the values of the eliminated Dirichlet DOFs.
    typename gsSparseSolver<T>::CGDiagonal solver;
    fixedDofs = solver.compute(globProjMat).solve(globProjRhs);

} // computeDirichletDofsL2Proj


} //namespace gismo
