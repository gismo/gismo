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

#include <gsAssembler/gsCPPInterface.h>

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
    const gsMultiPatch<T>* m_gmap;

    gsOptionList m_options;

    gsSparseMatrix<T> m_matrix;
    gsMatrix<T>       m_rhs;

    std::list<gsFeSpaceData<T> > m_sdata;
    std::vector<gsFeSpaceData<T>*> m_vrow;
    std::vector<gsFeSpaceData<T>*> m_vcol;

    typedef typename gsExprHelper<T>::nullExpr    nullExpr;

public:

    typedef typename gsSparseMatrix<T>::BlockView matBlockView;

    typedef typename gsSparseMatrix<T>::constBlockView matConstBlockView;

    typedef typename gsBoundaryConditions<T>::bcRefList   bcRefList;
    //typedef typename gsBoundaryConditions<T>::bcContainer bcContainer;
    typedef gsBoxTopology::ifContainer ifContainer;

    typedef typename gsExprHelper<T>::element     element;     ///< Current element
    typedef typename gsExprHelper<T>::geometryMap geometryMap; ///< Geometry map type
    typedef typename gsExprHelper<T>::variable    variable;    ///< Variable type
    typedef typename gsExprHelper<T>::space       space;       ///< Space type
    typedef typename expr::gsFeSolution<T>        solution;    ///< Solution type

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
      m_vrow(_rBlocks,nullptr), m_vcol(_cBlocks,nullptr)
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
    const gsMultiBasis<T> & integrationElements() const
    { return m_exprdata->multiBasis(); }

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
        s.setSpaceData(m_sdata.back());
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
            getSpace(*m_vrow[id]->fs,m_vrow[id]->dim());
        s.setSpaceData(*m_vrow[id]);
        return *m_vrow[id];
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
        clearMatrix();
    }

    void clearRhs() { m_rhs.setZero(); }

    void clearMatrix()
    {
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
                nz *= bdA * static_cast<T>(m_exprdata->multiBasis().maxDegree(i)) + static_cast<T>(bdB);

            m_matrix.reservePerColumn(numBlocks()*cast<T,index_t>(nz*(1.0+bdO)) );
        }
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

    template<class... expr> void assembleIfc(const ifContainer & iFaces, expr... args);
    /*
      template<class... expr> void collocate(expr... args);// eg. collocate(-ilapl(u), f)
    */

    /// \brief Assembles the Jacobian matrix of the expression \a args with
    // respect to the solution \a u
    template<class expr> void assembleJacobian(const expr residual, solution & u);

    template<class expr> void assembleJacobianIfc(const ifContainer & iFaces,
                                                  const expr residual, solution  u);

private:

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
                rowSizes[r] = m_vrow[r]->dim() * m_vrow[r]->mapper.freeSize();
            colSizes.resize(m_vcol.size());
            for (index_t c = 0; c != colSizes.size(); ++c) // for all col-blocks
                colSizes[c] = m_vcol[c]->dim() * m_vcol[c]->mapper.freeSize();
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
            const bool m = E::isMatrix();
            GISMO_ASSERT(v.isValid(), "The row space is not valid");
            GISMO_ASSERT(!m || u.isValid(), "The column space is not valid");
            GISMO_ASSERT(m || (ea.numDofs()==ee.rhs().size()), "The right-hand side vector is not initialized");
        }
    } _checkExpr;

    // Evaluates expression and and assembles global matrix/rhs
    struct _eval
    {
        gsSparseMatrix<T> & m_matrix;
        gsMatrix<T>       & m_rhs;
        const gsVector<T> & m_quWeights;
        bool m_elim;
        gsMatrix<T>         localMat;
        gsMatrix<T>         aux;

        _eval(gsSparseMatrix<T> & _matrix,
              gsMatrix<T>       & _rhs,
              const gsVector<>  & _quWeights)
        : m_matrix(_matrix), m_rhs(_rhs),
          m_quWeights(_quWeights), m_elim(true)
        { }

        void setElim(bool elim) {m_elim = elim;}

        template <typename E> void operator() (const gismo::expr::_expr<E> & ee)
        {
            // ------- Compute  -------
            quadrature(ee,localMat);

            //  ------- Accumulate  -------
            if (E::isMatrix())
                if (m_elim) push<true,true>(ee.rowVar(), ee.colVar());
                else push<true,false>(ee.rowVar(), ee.colVar());
            else if (E::isVector())
                if (m_elim) push<false,true>(ee.rowVar(), ee.colVar());
                else push<false,false>(ee.rowVar(), ee.colVar());
            else
            {
                GISMO_ERROR("Something went terribly wrong at this point");
                //GISMO_ASSERTrowSpan() && (!colSpan())
            }

        }// operator()

        template <typename E>
        inline void quadrature(const gismo::expr::_expr<E> & ee,
                               gsMatrix<T> & lm)
        {
            // ------- Compute  -------
            const T * w = m_quWeights.data();
            lm.noalias() = (*w) * ee.eval(0);
            for (index_t k = 1; k != m_quWeights.rows(); ++k)
                lm.noalias() += (*(++w)) * ee.eval(k);
        }

        template <typename E> void diff(const gismo::expr::_expr<E> & ee,
                                        solution & u)
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
                        //todo: take local copy of local solution u

                        //Perturb \a u
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

            //  ------- Accumulate  -------
            push<true,false>(ee.rowVar(), ee.rowVar());
        }

        void operator() (const expr::_expr<expr::gsNullExpr<T> > &) {}

        template<bool isMatrix, bool elim = true>
        void push(const expr::gsFeSpace<T> & v,
                  const expr::gsFeSpace<T> & u)
        {
            GISMO_ASSERT(v.isValid(), "The row space is not valid");
            GISMO_ASSERT(!isMatrix || u.isValid(), "The column space is not valid");
            //GISMO_ASSERT(isMatrix || (numDofs()==m_rhs.size()), "The right-hand side vector is not initialized");

            const index_t rd            = v.dim();//row
            const index_t cd            = u.dim();//col
            //const index_t rp            = v.data().patchId;
            //const index_t cp            = (isMatrix ? u.data().patchId : 0);
            const gsDofMapper  & rowMap = v.mapper();
            const gsDofMapper  & colMap = (isMatrix ? u.mapper() : rowMap);
            gsMatrix<index_t> & rowInd0 = const_cast<gsMatrix<index_t>&>(v.data().actives);
            gsMatrix<index_t> & colInd0 = (isMatrix ? const_cast<gsMatrix<index_t>&>(u.data().actives) : rowInd0);
            const gsMatrix<T> & fixedDofs = (isMatrix ? u.fixedPart() : gsMatrix<T>());

            if (isMatrix)
            {
                GISMO_ASSERT( rowInd0.rows()*rd==localMat.rows() && colInd0.rows()*cd==localMat.cols(),
                              "Invalid local matrix (expected "<<rowInd0.rows()*rd <<"x"<< colInd0.rows()*cd <<"), got\n" << localMat );

                GISMO_ASSERT( colMap.boundarySize()==fixedDofs.size(),
                              "Invalid values for fixed part");

                //GISMO_ASSERT( colMap.boundarySize()==0 || m_rhs.cols()==1,
                //              "Invalid values for fixed part");
            }

            for (index_t r = 0; r != rd; ++r)
            {
                const index_t rls = r * rowInd0.rows();     //local stride
                for (index_t i = 0; i != rowInd0.rows(); ++i)
                {
                    const index_t ii = rowMap.index(rowInd0.at(i),v.data().patchId,r); //N_i
                    if ( rowMap.is_free_index(ii) )
                    {
                        if (isMatrix)
                        {
                            for (index_t c = 0; c != cd; ++c)
                            {
                                const index_t cls = c * colInd0.rows();     //local stride

                                for (index_t j = 0; j != colInd0.rows(); ++j)
                                {
                                    if ( 0 == localMat(rls+i,cls+j) ) continue;

                                    const index_t jj = colMap.index(colInd0.at(j),u.data().patchId,c); // N_j
                                    if ( colMap.is_free_index(jj) )
                                    {
                                        // If matrix is symmetric, we could
                                        // store only lower triangular part
                                        //if ( (!symm) || jj <= ii )
#                                       pragma omp critical (acc_m_matrix)
                                        m_matrix.coeffRef(ii, jj) += localMat(rls+i,cls+j);
                                    }
                                    else if (elim) // colMap.is_boundary_index(jj) )
                                    {
                                        // Symmetric treatment of eliminated BCs
                                        // GISMO_ASSERT(1==m_rhs.cols(), "-");
#                                       pragma omp critical (acc_m_rhs)
                                        m_rhs.at(ii) -= localMat(rls+i,cls+j) *
                                            fixedDofs.at(colMap.global_to_bindex(jj));
                                    }
                                }
                            }
                        }
                        else
                        {
                            //The right-hand side can have more than one columns
#                           pragma omp critical (acc_m_rhs)
                            m_rhs.row(ii) += localMat.row(rls+i);
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
    opt.addInt("DirichletValues"  , "Method for computation of Dirichlet DoF values [100..103]", 101);// deprecated
    opt.addInt("DirichletStrategy"  , "Strategy related to enforcement of Dirichlet DoF values [0ignore, 1:eliminate, 2:..]", 1);
    opt.addReal("quA", "Number of quadrature points: quA*deg + quB", 1.0  );
    opt.addInt ("quB", "Number of quadrature points: quA*deg + quB", 1    );
    opt.addReal("bdA", "Estimated nonzeros per column of the matrix: bdA*deg + bdB", 2.0  );
    opt.addInt ("bdB", "Estimated nonzeros per column of the matrix: bdA*deg + bdB", 1    );
    opt.addReal("bdO", "Overhead of sparse mem. allocation: (1+bdO)(bdA*deg + bdB) [0..1]", 0.333);
    opt.addSwitch("flipSide", "Flip side of interface where integration is performed.", false);
    opt.addSwitch("movingInterface", "Used in interface assembly when interface is not stationary.", false);
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
        if (!m_vcol.front()->valid()) m_vcol.front()->init();
        m_vcol[i]->mapper.setShift(m_vcol[i-1]->mapper.firstIndex() +
                                   m_vcol[i-1]->dim*m_vcol[i-1]->mapper.freeSize() );

        if ( m_vcol[i] != m_vrow[i] )
        {
            if (!m_vrow.front()->valid()) m_vrow.front()->init();
            m_vrow[i]->mapper.setShift(m_vrow[i-1]->mapper.firstIndex() +
                                       m_vrow[i-1]->dim*m_vrow[i-1]->mapper.freeSize() );
        }
    }
}

template<size_t I, class op, typename... Ts>
void op_tuple_impl (op _op, const std::tuple<Ts...> &tuple)
{
    _op(std::get<I>(tuple));
    if (I + 1 < sizeof... (Ts))
        op_tuple_impl<(I+1 < sizeof... (Ts) ? I+1 : I)> (_op, tuple);
}

template<class op, typename... Ts>
void op_tuple (op _op, const std::tuple<Ts...> &tuple)
{ op_tuple_impl<0>(_op,tuple); }

template<class T>
template<class... expr>
void gsExprAssembler<T>::assemble(const expr &... args)
{
    GISMO_ASSERT(matrix().cols()==numDofs(), "System not initialized");

#pragma omp parallel
{
#   ifdef _OPENMP
    const int tid = omp_get_thread_num();
    const int nt  = omp_get_num_threads();
#   endif
    auto arg_tpl = std::make_tuple(args...);

    m_exprdata->parse(arg_tpl);
    m_exprdata->activateFlags(SAME_ELEMENT);
    //op_tuple(__printExpr(), arg_tpl);

    typename gsQuadRule<T>::uPtr QuRule; // Quadrature rule

    gsVector<T> quWeights; // quadrature weights
    _eval ee(m_matrix, m_rhs, quWeights);
    const index_t elim = m_options.getInt("DirichletStrategy");
    ee.setElim(1==elim);

    // Note: omp thread will loop over all patches and will work on Ep/nt
    // elements, where Ep is the elements on the patch.
    for (unsigned patchInd = 0; patchInd < m_exprdata->multiBasis().nBases(); ++patchInd) //todo: distribute in parallel somehow?
    {
        QuRule = gsQuadrature::getPtr(m_exprdata->multiBasis().basis(patchInd), m_options);

        // Initialize domain element iterator for current patch
        typename gsBasis<T>::domainIter domIt =  // add patchInd to domainiter ?
            m_exprdata->multiBasis().basis(patchInd).makeDomainIterator();
        m_exprdata->getElement().set(*domIt,quWeights);

        // Start iteration over elements of patchInd
#       ifdef _OPENMP
        for ( domIt->next(tid); domIt->good(); domIt->next(nt) )
#       else
        for (; domIt->good(); domIt->next() )
#       endif
        {
            // Map the Quadrature rule to the element
            QuRule->mapTo( domIt->lowerCorner(), domIt->upperCorner(),
                           m_exprdata->points(), quWeights);

            if (m_exprdata->points().cols()==0)
                continue;

            // Perform required pre-computations on the quadrature nodes
            m_exprdata->precompute(patchInd);
            //m_exprdata->precompute(patchInd, QuRule, *domIt); // todo

            // Assemble contributions of the element
            op_tuple(ee, arg_tpl);
        }
    }

}//omp parallel
    m_matrix.makeCompressed();
}

template<class T>
template<class... expr>
void gsExprAssembler<T>::assembleBdr(const bcRefList & BCs, expr&... args)
{
    GISMO_ASSERT(matrix().cols()==numDofs(), "System not initialized");

    if ( BCs.empty() || 0==numDofs() ) return;
    m_exprdata->setMutSource(*BCs.front().get().function()); //initialize once

// #pragma omp parallel
// {
// #   ifdef _OPENMP
//     const int tid = omp_get_thread_num();
//     const int nt  = omp_get_num_threads();
// #   endif
    auto arg_tpl = std::make_tuple(args...);
    m_exprdata->parse(arg_tpl);
    m_exprdata->activateFlags(SAME_ELEMENT);

    typename gsQuadRule<T>::uPtr QuRule; // Quadrature rule
    gsVector<T> quWeights;               // quadrature weights

    _eval ee(m_matrix, m_rhs, quWeights);

//#   pragma omp parallel for
    for (typename bcRefList::const_iterator iit = BCs.begin(); iit!= BCs.end(); ++iit)
    {
        const boundary_condition<T> * it = &iit->get();

        QuRule = gsQuadrature::getPtr(m_exprdata->multiBasis().basis(it->patch()), m_options, it->side().direction());

        // Update boundary function source
        m_exprdata->setMutSource(*it->function());

        typename gsBasis<T>::domainIter domIt =
            m_exprdata->multiBasis().basis(it->patch()).
            makeDomainIterator(it->side());
        m_exprdata->getElement().set(*domIt,quWeights);

        // Start iteration over elements
        for (; domIt->good(); domIt->next() )
        {
            // Map the Quadrature rule to the element
            QuRule->mapTo( domIt->lowerCorner(), domIt->upperCorner(),
                           m_exprdata->points(), quWeights);

            if (m_exprdata->points().cols()==0)
                continue;

            // Perform required pre-computations on the quadrature nodes
            m_exprdata->precompute(it->patch(), it->side());

            // Assemble contributions of the element
            op_tuple(ee, arg_tpl);
        }
    }

//}//omp parallel

    m_matrix.makeCompressed();
}

template<class T> template<class... expr>
void gsExprAssembler<T>::assembleIfc(const ifContainer & iFaces, expr... args)
{
    GISMO_ASSERT(matrix().cols()==numDofs(), "System not initialized");

#pragma omp parallel
{
    typedef typename gsFunction<T>::uPtr ifacemap;

    auto arg_tpl = std::make_tuple(args...);

    m_exprdata->parse(arg_tpl);
    m_exprdata->activateFlags(SAME_ELEMENT); //note: SAME_ELEMENT is 0 at the opposite/mirrored patch

    typename gsQuadRule<T>::uPtr QuRule;
    gsVector<T> quWeights;// quadrature weights
    _eval ee(m_matrix, m_rhs, quWeights);

    const bool flipSide = m_options.askSwitch("flipSide", false);

    ifacemap interfaceMap;
#   pragma omp for
    for (gsBoxTopology::const_iiterator it = iFaces.begin();
         it != iFaces.end(); ++it )
    {
        // If flipSide switch is enabled, then the integration will be
        // performed on the opposite side of the interface
        const boundaryInterface & iFace =  flipSide ? it->getInverse() : *it;
        const index_t patch1 = iFace.first() .patch;
        const index_t patch2 = iFace.second().patch;

        if (iFace.type() == interaction::conforming)
            interfaceMap = gsAffineFunction<T>::make( iFace.dirMap(), iFace.dirOrientation(),
                                                      m_exprdata->multiBasis().basis(patch1).support(),
                                                      m_exprdata->multiBasis().basis(patch2).support() );
        else
            interfaceMap = gsCPPInterface<T>::make(getGeometryMap(), m_exprdata->multiBasis(), iFace);

        QuRule = gsQuadrature::getPtr(m_exprdata->multiBasis().basis(patch1),
                                   m_options, iFace.first().side().direction());

        typename gsBasis<T>::domainIter domIt =
            m_exprdata->multiBasis().basis(patch1)
            .makeDomainIterator(iFace.first().side());
        m_exprdata->getElement().set(*domIt, quWeights);

        // Start iteration over elements
        for (; domIt->good(); domIt->next() )
        {
            // Map the Quadrature rule to the element
            QuRule->mapTo( domIt->lowerCorner(), domIt->upperCorner(),
                           m_exprdata->points(), quWeights);
            interfaceMap->eval_into(m_exprdata->points(), m_exprdata->pointsIfc());

            if (m_exprdata->points().cols()==0)
                continue;

            // Perform required pre-computations on the quadrature nodes
            m_exprdata->precompute(iFace);

            //eg.
            // uL*vL/2 + uR*vL/2  - uL*vR/2 - uR*vR/2
            //[ B11 B21 ]
            //[ B12 B22 ]

            op_tuple(ee, arg_tpl);
        }
    }

}//omp parallel
    m_matrix.makeCompressed();
}

template<class T> template<class expr>
void gsExprAssembler<T>::assembleJacobian(const expr residual, solution & u)
{
    GISMO_ASSERT(matrix().cols()==numDofs(), "System not initialized");
    GISMO_ASSERT(expr::isVector(), "Expecting a vector expression.");

    clearMatrix();
    clearRhs();

#pragma omp parallel
{
#   ifdef _OPENMP
    const int tid = omp_get_thread_num();
    const int nt  = omp_get_num_threads();
#   endif

    m_exprdata->parse(residual, u);
    m_exprdata->activateFlags(SAME_ELEMENT);
    //op_tuple(__printExpr(), arg_tpl);

    typename gsQuadRule<T>::uPtr QuRule; // Quadrature rule  ---->OUT

    gsVector<T> quWeights; // quadrature weights

    _eval ee(m_matrix, m_rhs, quWeights);

    // Note: omp thread will loop over all patches and will work on Ep/nt
    // elements, where Ep is the elements on the patch.
    for (unsigned patchInd = 0; patchInd < m_exprdata->multiBasis().nBases(); ++patchInd)
    {
        QuRule = gsQuadrature::getPtr(m_exprdata->multiBasis().basis(patchInd), m_options);

        // Initialize domain element iterator for current patch
        typename gsBasis<T>::domainIter domIt =  // add patchInd to domainiter ?
            m_exprdata->multiBasis().basis(patchInd).makeDomainIterator();
        m_exprdata->getElement().set(*domIt,quWeights);

        // Start iteration over elements of patchInd
#       ifdef _OPENMP
        for ( domIt->next(tid); domIt->good(); domIt->next(nt) )
#       else
        for (; domIt->good(); domIt->next() )
#       endif
        {
            // Map the Quadrature rule to the element
            QuRule->mapTo( domIt->lowerCorner(), domIt->upperCorner(),
                           m_exprdata->points(), quWeights);

            if (m_exprdata->points().cols()==0)
                continue;

            // Evaluate at quadrature points
            m_exprdata->precompute(patchInd);

#           pragma omp critical (assemble_fdiffs)
            {
                // ee(residual); //Computes residual to m_rhs
                ee.diff(residual, u); //Computes Jacobian
            }
        }
    }

}//omp parallel
    m_matrix.makeCompressed();
}


template<class T> template<class expr>
void gsExprAssembler<T>::assembleJacobianIfc(const ifContainer & iFaces,
                                             const expr residual, solution  u)
{
    GISMO_ASSERT(matrix().cols()==numDofs(), "System not initialized");
    GISMO_ASSERT(expr::isVector(), "Expecting a vector expression.");

    // clearMatrix();

    m_exprdata->parse(residual, u);
    m_exprdata->activateFlags(SAME_ELEMENT);
    //op_tuple(__printExpr(), arg_tpl);

    typename gsQuadRule<T>::uPtr QuRule; // Quadrature rule
    gsVector<T> quWeights; // quadrature weights

    _eval ee(m_matrix, m_rhs, quWeights);
    static const T delta = 0.00001;
    const bool flipSide = m_options.askSwitch("flipSide", false);
    const bool movingInterface = m_options.askSwitch("movingInterface", false);

    for (gsBoxTopology::const_iiterator it = iFaces.begin();
         it != iFaces.end(); ++it )
    {
        // If flipSide switch is enabled, then the integration will be
        // performed on the opposite side of the interface
        const boundaryInterface & iFace =  flipSide ? it->getInverse() : *it;
        const index_t patch1 = iFace.first() .patch;
        //const index_t patch2 = iFace.second().patch;

        gsCPPInterface<T> interfaceMap(getGeometryMap(), // ! make current geometry
                                       m_exprdata->multiBasis(), iFace);

        QuRule = gsQuadrature::getPtr(m_exprdata->multiBasis().basis(patch1),
                                   m_options, iFace.first().side().direction());

        typename gsBasis<T>::domainIter domIt =
            m_exprdata->multiBasis().basis(patch1)
            .makeDomainIterator(iFace.first().side());
        m_exprdata->getElement().set(*domIt, quWeights);

        // Start iteration over elements
        for (; domIt->good(); domIt->next() )
        {
            // Map the Quadrature rule to the element
            QuRule->mapTo( domIt->lowerCorner(), domIt->upperCorner(),
                           m_exprdata->points(), quWeights);
            interfaceMap.eval_into(m_exprdata->points(), m_exprdata->pointsIfc());

            if (m_exprdata->points().cols()==0)
                continue;

            // Perform required pre-computations on the quadrature nodes
            m_exprdata->precompute(iFace);

            // ee(residual);
            if (!movingInterface)
            {
                ee.diff(residual, u);
            }
            else // For Moving Geometry Maps
            {
                GISMO_ASSERT(expr::isVector(), "Expecting a vector expression.");
                static const T delta = 0.00001;

                const index_t sz = residual.eval(0).rows(); // Caution this must be the correct size , depending if .right() or .left() are used
                ee.localMat.setZero(sz, sz);

                auto & rowVar = residual.rowVar();

                for ( index_t c=0; c!= u.dim(); c++)
                {
                    const index_t rls = c * rowVar.data().actives.rows();     //local stride
                    for ( index_t j = 0; j != sz/u.dim(); j++ )     // for all basis functions (col(j))
                    {
                        const index_t jj = u.mapper().index(rowVar.data().actives(j),
                                                            rowVar.data().patchId, c);
                        if (rowVar.mapper().is_free_index(jj) )
                        {
                            //Perturb \a u
                            u.perturbLocal( delta  , jj, rowVar.data().patchId);
                            interfaceMap.updateBdr();
                            interfaceMap.eval_into(m_exprdata->points(), m_exprdata->pointsIfc());
                            m_exprdata->precompute(iFace);
                            ee.quadrature(residual, ee.aux);
                            ee.localMat.col(rls+j) += 8 * ee.aux;


                            u.perturbLocal( delta  , jj, rowVar.data().patchId);
                            interfaceMap.updateBdr();
                            interfaceMap.eval_into(m_exprdata->points(), m_exprdata->pointsIfc());
                            m_exprdata->precompute(iFace);
                            ee.quadrature(residual, ee.aux);
                            ee.localMat.col(rls+j) -= ee.aux;

                            u.perturbLocal(-3*delta, jj, rowVar.data().patchId);
                            interfaceMap.updateBdr();
                            interfaceMap.eval_into(m_exprdata->points(), m_exprdata->pointsIfc());
                            m_exprdata->precompute(iFace);
                            ee.quadrature(residual, ee.aux);
                            ee.localMat.col(rls+j) -= 8 * ee.aux;

                            u.perturbLocal( -delta , jj, rowVar.data().patchId);
                            interfaceMap.updateBdr();
                            interfaceMap.eval_into(m_exprdata->points(), m_exprdata->pointsIfc());
                            m_exprdata->precompute(iFace);
                            ee.quadrature(residual, ee.aux);
                            ee.localMat.col(rls+j) += ee.aux;

                            ee.localMat.col(rls+j) /= 12*delta;
                            //Unperturb \a u
                            u.perturbLocal(2*delta, jj, rowVar.data().patchId);
                            interfaceMap.updateBdr();
                        }
                    }
                }

                //  ------- Accumulate  -------
                ee.template push<true,false>(residual.rowVar(), residual.rowVar());
            }
        }
    }

    m_matrix.makeCompressed();
}


} //namespace gismo
