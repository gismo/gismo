/** @file gsExprAssembler.h

    @brief Generic expressions matrix assembly

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris
*/

#pragma once

#include<gsAssembler/gsExprHelper.h>
#include <gsUtils/gsPointGrid.h>


//#include<tr1/tuple>

/*
namespace
{

template<class T>
struct __eval
{
    std::vector<gismo::gsMatrix<T> > & m;
    const gismo::gsVector<T> & w;
    index_t i;
    
    __eval(std::vector<gismo::gsMatrix<T> > & _m,
           const gismo::gsVector<T> & _w) : m(_m), w(_w), i(0) { }
    
    template<typename E>
    void operator () (const gismo::expr::_expr<E> & t)
    {
        gismo::gsMatrix<T> & mat = m[i++];        
        mat = w[0] * t.eval(0);
        for (index_t k = 1; k != w.rows(); ++k)
            mat += w[k] * t.eval(k);
    }
};

struct __setFlag
{
    template<typename E> //static
    void operator () (const gismo::expr::_expr<E> & t) { t.setFlag(); }
};

template<class T>
struct __local2global
{
    std::vector<gismo::gsMatrix<T> > & m;
    gismo::gsSparseSystem<T> & sys;
    const index_t patchInd;
    index_t i;
    
    __local2global(std::vector<gismo::gsMatrix<T> > & _m,
                   gismo::gsSparseSystem<T> & _sys,
                   const index_t _patchInd)
    : sys(_sys), m(_m), patchInd(_patchInd), i(0) { }

    template<typename E> // bool Col = E::col
    void operator () (const gismo::expr::_expr<E> & t)
    {
        
        gismo::gsMatrix<T> & mat = m[i++];
        const gismo::expr::gsFeSpace<T> & v = t.rowVar();
        //const gismo::expr::gsFeSpace<T> & u = t.colVar();
        
        //gismo::gsSpAcc<Col>(sys.matrix(), sys.rhs(), mat, patchInd)
    }
};

struct __print
{
    template<typename E> //static
    void operator () (const gismo::expr::_expr<E> & t) { t.print(gsInfo); }
};

template<std::size_t I, typename Tuple, typename Func> inline
typename gismo::internal::enable_if< I != std::tr1::tuple_size<Tuple>::value,void>::type
__for_each(const Tuple& tuple, Func func)
{
    func(std::tr1::get<I>(tuple));
    __for_each<I + 1, Tuple,Func>(tuple, func);
}

template<std::size_t I, typename Tuple, typename Func> inline
typename gismo::internal::enable_if< I == std::tr1::tuple_size<Tuple>::value,void>::type
__for_each(const Tuple& tuple, Func func) { }

template<typename Tuple, typename Func> inline void
for_each_expr(const Tuple& tuple, Func func) { __for_each<0, Tuple,Func>(tuple, func); }

}
*/

namespace gismo
{

template<class T, bool left, bool right>
void gsSpAcc(
    gsSparseMatrix<T> & m_matrix,
    gsMatrix<T> & m_rhs,
    const gsMatrix<T> & localMat,
    const gsMatrix<T> & localRhs,
    const expr::gsFeSpace<T> & v, //rowvar
    const expr::gsFeSpace<T> & u, //colvar
    const index_t patchInd)
{
    const index_t cd            = u.dim();
    const index_t rd            = v.dim();
    const gsDofMapper  & colMap = u.mapper();
    const gsDofMapper  & rowMap = v.mapper();
    gsMatrix<unsigned> & colInd0 = const_cast<gsMatrix<unsigned>&>(u.data().actives);
    gsMatrix<unsigned> & rowInd0 = const_cast<gsMatrix<unsigned>&>(v.data().actives);
    const gsMatrix<T>  & fixedDofs = u.fixedPart();

    gsMatrix<unsigned> rowInd, colInd;
    //if (&colMap==&rowMap) && (&rowInd==&colInd)
    colMap.localToGlobal(colInd0, patchInd, colInd);
    rowMap.localToGlobal(rowInd0, patchInd, rowInd);

    //gsDebugVar( localMat.dim() );
    // colMap.print();
    // rowMap.print();
    
    //GISMO_STATIC_ASSERT(left || right, "Nothing to do.");
    GISMO_ASSERT( !left || colMap.boundarySize()==fixedDofs.rows(),
                  "Invalid values for fixed part");
    
    for (index_t r = 0; r != rd; ++r)
    {
        const index_t rls = r * rowInd.rows();     //local stride
        const index_t rgs = r * rowMap.freeSize(); //global stride
        
        for (index_t i = 0; i != rowInd.rows(); ++i) // **
        {
            const int ii = rgs + rowInd.at(i); // N_i
            
            if ( rowMap.is_free_index(rowInd.at(i)) )
            {
                if (right)
                {
                    m_rhs.row(ii) += localRhs.row(rls+i);
                }
                
                if (left)
                    for (index_t c = 0; c != cd; ++c)
                    {
                        const index_t cls = c * colInd.rows();     //local stride
                        const index_t cgs = c * colMap.freeSize(); //global stride
                        
                        for (index_t j = 0; j != colInd.rows(); ++j)
                        {
                            if ( 0 == localMat(rls+i,cls+j) ) continue;
                            
                            const int jj = cgs + colInd.at(j); // N_j
                            
                            if ( colMap.is_free_index(colInd.at(j)) )
                            {
                                // If matrix is symmetric, we store only lower
                                // triangular part
                                //if ( (!symm) || jj <= ii )
                                m_matrix.coeffRef(ii, jj) +=
                                    localMat(rls+i,cls+j);
                            }
                            else // colMap.is_boundary_index(jj) )
                            {
                                // Symmetric treatment of Dirichlet BC
                                // GISMO_ASSERT(1==m_rhs.cols(), "-");
                                m_rhs.at(ii) -= localMat(rls+i,cls+j) *
                                    fixedDofs(colMap.global_to_bindex(colInd.at(j)), c);
                            }
                        }
                    }
            }
        }
    }
}


/**
   Generic assembler of (matrix or vector) isogeometric expressions
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

public:

    typedef typename gsSparseMatrix<T>::BlockView matBlockView;
    
    typedef typename gsBoundaryConditions<T>::bcContainer bcContainer;
    
    typedef std::vector< boundaryInterface > intContainer;

    typedef typename gsExprHelper<T>::element     element;
    typedef typename gsExprHelper<T>::geometryMap geometryMap;
    typedef typename gsExprHelper<T>::variable    variable;
    typedef typename gsExprHelper<T>::space       space;
    typedef typename expr::gsFeSolution<T>        solution;
    typedef typename gsExprHelper<T>::nullExpr    nullExpr;

public:

    gsExprAssembler(int _rBlocks = 1, int _cBlocks = 1)
    : m_exprdata(gsExprHelper<T>::New()), m_options(defaultOptions()), 
      m_vrow(_rBlocks,NULL), m_vcol(_cBlocks,NULL)
    { }

    /// @brief Returns the list of default options for assembly
    static gsOptionList defaultOptions();

    index_t numDofs()
    {
        GISMO_ASSERT( m_vcol.back()->mapper().isFinalized(),
                      "initSystem() has not been called.");
        return m_vcol.back()->mapper().firstIndex() + 
            m_vcol.back()->dim() * m_vcol.back()->mapper().freeSize();
    }

    index_t numTestDofs()
    {
        GISMO_ASSERT( m_vrow.back()->mapper().isFinalized(),
                      "initSystem() has not been called.");
        return m_vrow.back()->mapper().firstIndex() + 
            m_vrow.back()->dim() * m_vrow.back()->mapper().freeSize();
    }

    index_t numBlocks() const
    {
        index_t nb = 0;
        for (size_t i = 0; i!=m_vrow.size(); ++i)
            nb += m_vrow[i]->dim();
        return nb;
    }

    gsOptionList & options() {return m_options;}
    
    /// @brief Returns the left-hand global matrix
    const gsSparseMatrix<T> & matrix() const { return m_matrix; }

    /// @brief Returns the right-hand side vector(s)
    const gsMatrix<T> & rhs() const { return m_rhs; }

    void setIntegrationElements(const gsMultiBasis<T> & mesh)
    { m_exprdata->setMultiBasis(mesh); }
    
    const gsMultiBasis<T> & integrationElements() const
    { return m_exprdata->multiBasis(); }

    const typename gsExprHelper<T>::Ptr exprData() const { return m_exprdata; }
    
    geometryMap setMap(const gsMultiPatch<T> & mp) //conv->tmp->error
    { return m_exprdata->setMap(mp); }

    geometryMap setMap(const gsGeometry<T> & mp)
    { return m_exprdata->setMap(mp); }

    variable setSpace(const gsFunctionSet<T> & mp, const gsBoundaryConditions<T> & bc,
                      index_t dim = 1, index_t id = 0)
    {
        GISMO_ASSERT(1==mp.targetDim(), "Expecting scalar source space");
        GISMO_ASSERT(static_cast<size_t>(id)<m_vrow.size(),
                     "Given ID "<<id<<" exceeds "<<m_vrow.size() );
        expr::gsFeSpace<T> & u = m_exprdata->setSpace(mp,dim);
        u.setId(id);
        u.setBc(bc);
        m_vrow[id] = m_vcol[id] = &u;
        return u;
    }

    variable setSpace(const gsFunctionSet<T> & mp, index_t dim = 1, index_t id = 0)
    {
        GISMO_ASSERT(1==mp.targetDim(), "Expecting scalar source space");        
        GISMO_ASSERT(static_cast<size_t>(id)<m_vrow.size(),
                     "Given ID "<<id<<" exceeds "<<m_vrow.size()-1 );
        expr::gsFeSpace<T> & u = m_exprdata->setSpace(mp,dim);
        u.setId(id);
        m_vrow[id] = m_vcol[id] = &u;
        return u;
    }

    void computeDirichletDofs2(int unk);
    void computeDirichletDofsIntpl2(const expr::gsFeSpace<T> & u);
    void computeDirichletDofsL2Proj(const expr::gsFeSpace<T> & u);
	void setFixedDofVector(gsMatrix<T> & dof, int unk = 0);
	void setFixedDofs(const gsMatrix<T> & coefMatrix, int unk = 0, int patch = 0);
    
    // for Petrov-Galerkin
    variable setTestSpace(variable u, const gsFunctionSet<T> & mp)
    {
        //GISMO_ASSERT(0!=u.mapper(), "Not a space"); // done on initSystem
        expr::gsFeSpace<T> & s = m_exprdata->setSpace(mp,u.dim());
        space uu = static_cast<space>(u);
        s.setId(uu.id());
        s.setBc(uu.bc());
        m_vrow[s.id()] = &s;
        return s;
    }

    variable setCoeff(const gsFunctionSet<T> & func)
    { return m_exprdata->setVar(func, 1); }

    variable setCoeff(const gsFunctionSet<T> & func, geometryMap G)
    { return m_exprdata->setVar(func,G); }

    solution getSolution(variable u, gsMatrix<T> & cf) const
    {
        // todo: if (m_exprdata->isSpace(u));
        space s = static_cast<space>(u);
        return solution(s, cf);
        //GISMO_ERROR("Need a space variable");
    }

    variable getBdrFunction() const { return m_exprdata->getMutVar(); }

    element getElement() const { return m_element; }

    void initSystem()
    {
        // Check spaces.nPatches==mesh.patches
        initMatrix();
        m_rhs.setZero(numDofs(), 1);

        for (size_t i = 0; i!= m_vcol.size(); ++i)
            computeDirichletDofs2(i);
    }
    
    void initMatrix()
    {
        resetDimensions();
        m_matrix = gsSparseMatrix<T>(numTestDofs(), numDofs());
        
        if ( 0 == m_matrix.size() )
            gsWarn << " No internal DOFs, zero sized system.\n";
        else
        {
            // Pick up values from options
            const T bdA       = m_options.getReal("bdA");
            const index_t bdB = m_options.getInt("bdB");
            const T bdO       = m_options.getReal("bdO");
            T nz = 1;
            const index_t dim = m_exprdata->multiBasis().domainDim();
            for (index_t i = 0; i != dim; ++i)
                nz *= bdA * m_exprdata->multiBasis().maxDegree(i) + bdB;

            m_matrix.reservePerColumn(numBlocks()*cast<T,index_t>(nz*(1.0+bdO)) );
        }
    }
    
    void initVector(const index_t numRhs = 1)
    {
        resetDimensions();
        m_rhs.setZero(numDofs(), numRhs);
    }

    void resetDimensions()
    {
        for (size_t i = 0; i!=m_vcol.size(); ++i)
        {
            GISMO_ASSERT(NULL!=m_vcol[i], "Not set.");
            expr::gsFeSpace<T> & u = *m_vcol[i];

            const gsBoundaryConditions<> & ubc = u.hasBc() ? u.bc() : gsBoundaryConditions<>();
            
            if (const gsMultiBasis<T> * mb =
                dynamic_cast<const gsMultiBasis<T>*>(&u.source()) )
            {	
				mb->getMapper(
                    (dirichlet::strategy)(m_options.getInt("DirichletStrategy")),
                    (iFace::strategy)(m_options.getInt("InterfaceStrategy")),
                    ubc, u.mapper(), u.id(), true);
                //u.mapper().print();
            }
            else if (const gsBasis<T> * b =
                dynamic_cast<const gsBasis<T>*>(&u.source()) )
            {
                gsMultiBasis<T> mbb(*b);
                mbb.getMapper(
                    (dirichlet::strategy)(m_options.getInt("DirichletStrategy")),
                    (iFace::strategy)(m_options.getInt("InterfaceStrategy")),
                    ubc, u.mapper(), u.id(), true);
            }
            else
            {
                gsWarn<<"Problem initializing.\n";
            }
            
            if ( m_vcol[i] != m_vrow[i] )
            {
                GISMO_ASSERT(NULL!=m_vrow[i], "Not set.");
                expr::gsFeSpace<T> & v = *m_vrow[i];
                
                const gsBoundaryConditions<> & vbc = v.hasBc() ? v.bc() : gsBoundaryConditions<>();
                
                if (const gsMultiBasis<T> * mb =
                    dynamic_cast<const gsMultiBasis<T>*>(&v.source()) )
                {
                    mb->getMapper(
                        (dirichlet::strategy)(m_options.getInt("DirichletStrategy")),
                        (iFace::strategy)(m_options.getInt("InterfaceStrategy")),
                        vbc, v.mapper(), v.id(), true);
                }
                else if (const gsBasis<T> * b =
                         dynamic_cast<const gsBasis<T>*>(&v.source()) )
                {
                    gsMultiBasis<T> mbb(*b);
                    mbb.getMapper(
                        (dirichlet::strategy)(m_options.getInt("DirichletStrategy")),
                        (iFace::strategy)(m_options.getInt("InterfaceStrategy")),
                        vbc, v.mapper(), v.id(), true);
                }
                else
                {
                    gsWarn<<"Problem initializing.\n";
                }

            }
//            GISMO_ASSERT(m_vrow[i]==m_vcol[i], "to do"); // OK
        }

        for (size_t i = 1; i!=m_vcol.size(); ++i)
        {
            m_vcol[i]->mapper().setShift(m_vcol[i-1]->mapper().firstIndex() + 
                                         m_vcol[i-1]->dim()*m_vcol[i-1]->mapper().freeSize() );
                
            if ( m_vcol[i] != m_vrow[i] )
                m_vrow[i]->mapper().setShift(m_vrow[i-1]->mapper().firstIndex() +
                                             m_vrow[i-1]->dim()*m_vrow[i-1]->mapper().freeSize() );
        }
    }

    matBlockView matrixBlockView()
    {
        GISMO_ASSERT( m_vcol.back()->mapper().isFinalized(),
                      "initSystem() has not been called.");

        gsVector<index_t> rowSizes(m_vrow.size()), colSizes(m_vcol.size());

        for (index_t r = 0; r != rowSizes.size(); ++r) // for all row-blocks
            rowSizes[r] = m_vrow[r]->dim() * m_vrow[r]->mapper().freeSize();

        for (index_t c = 0; c != colSizes.size(); ++c) // for all col-blocks
            colSizes[c] = m_vcol[c]->dim() * m_vcol[c]->mapper().freeSize();

        return m_matrix.blockView(rowSizes,colSizes);
    }

    /*
    template<class Tpl>
    void assemble_tuple(const Tpl & a)
    {
        // util::tuple
        // util::tie
        gsDebugVar("-------------");
        gsDebugVar( std::tr1::tuple_size<Tpl>::value ); 
        for_each_expr(a, __print() );
        // ex.setflag();
        // lmat[i] = qu * ex.eval(0);
        gsDebugVar("-------------");
        // ex.rowVar();
        // ex.colVar();
    }
    
    template<class E1, class E2>
    void assemble1(const expr::_expr<E1> & t1, const expr::_expr<E2> & t2)
    {
        assemble_tuple( std::tr1::tie(t1,t2) );
    }
    
    template<class E1, class E2, class E3>
    void assemble1(const expr::_expr<E1> & t1,
                   const expr::_expr<E2> & t2, const expr::_expr<E3> & t3)
    {
        assemble_tuple( std::tr1::tie(t1,t2,t3) );
    }
    */
    
    void setOptions(gsOptionList opt) { m_options = opt; } // gsOptionList opt
    // .swap(opt) todo
        
    //void refresh() { }

    // assembleMixed(..) 2x2, diag blocks

    //template<class ... Types> void assemble(Types ... args);
    
    template<class E1, class E2>
    void assembleLhsRhs(const expr::_expr<E1> & exprLhs, const expr::_expr<E2> & exprRhs)
    {
        //GISMO_ASSERT(dynamic_cast<space>(exprLhs.colVar()), "Cannot deduce columns");
        space rvar = static_cast<space>(exprLhs.rowVar());
        GISMO_ASSERT(m_exprdata->exists(rvar), "Error - inexistent variable.");
        space cvar = static_cast<space>(exprLhs.colVar());
        GISMO_ASSERT(m_exprdata->exists(cvar), "Error - inexistent variable.");
        GISMO_ASSERT(&rvar==&exprRhs.rowVar(), "Inconsistent left and right hand side");
        assembleLhsRhs_impl<true,true>(exprLhs, exprRhs, rvar, cvar);
    }

    template<class E1>
    void assembleLhs(const expr::_expr<E1> & exprLhs)
    {
        space rvar = static_cast<space>(exprLhs.rowVar());
        GISMO_ASSERT(m_exprdata->exists(rvar), "Error - inexistent variable.");
        space cvar = static_cast<space>(exprLhs.colVar());
        GISMO_ASSERT(m_exprdata->exists(cvar), "Error - inexistent variable.");        
        assembleLhsRhs_impl<true,false>(exprLhs, nullExpr(), rvar, cvar);
    }

    template<class E2>
    void assembleRhs(const expr::_expr<E2> & exprRhs)
    {
        space rvar = static_cast<space>(exprRhs.rowVar());
        assembleLhsRhs_impl<false,true>(nullExpr(), exprRhs, rvar, rvar);
    }
    
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
    void assembleRhsInterface(const expr::_expr<E1> & exprInt, const intContainer & iFaces)
    {
        space rvar = static_cast<space>(exprInt.rowVar());
        GISMO_ASSERT(m_exprdata->exists(rvar), "Error - inexistent variable.");
        assembleInterface_impl<false,true>(nullExpr(), exprInt, rvar, rvar, iFaces);
    }

private:
    
    template<bool left, bool right, class E1, class E2>
    void assembleLhsRhs_impl(const expr::_expr<E1> & exprLhs,
                             const expr::_expr<E2> & exprRhs,
                             space rvar, space cvar);
    
    template<bool left, bool right, class E1, class E2>
    void assembleLhsRhsBc_impl(const expr::_expr<E1> & exprLhs,
                               const expr::_expr<E2> & exprRhs,
                               space rvar, space cvar,
                               const bcContainer & BCs);

    template<bool left, bool right, class E1, class E2>
    void assembleInterface_impl(const expr::_expr<E1> & exprLhs,
                                const expr::_expr<E2> & exprRhs,
                                space rvar, space cvar,
                                const intContainer & iFaces);
    
}; // gsExprAssembler

template<class T>
gsOptionList gsExprAssembler<T>::defaultOptions()
{
    gsOptionList opt;
    opt.addInt("DirichletStrategy", "Method for enforcement of Dirichlet BCs [11..14]", 11 );
    opt.addInt("DirichletValues"  , "Method for computation of Dirichlet DoF values [100..103]", 101);
    opt.addInt("InterfaceStrategy", "Method of treatment of patch interfaces [0..3]", 1  );
    opt.addReal("quA", "Number of quadrature points: quA*deg + quB", 1.0  );
    opt.addInt ("quB", "Number of quadrature points: quA*deg + quB", 1    );
    opt.addReal("bdA", "Estimated nonzeros per column of the matrix: bdA*deg + bdB", 2.0  );
    opt.addInt ("bdB", "Estimated nonzeros per column of the matrix: bdA*deg + bdB", 1    );
    opt.addReal("bdO", "Overhead of sparse mem. allocation: (1+bdO)(bdA*deg + bdB) [0..1]", 0.333);
    return opt;
}

template<class T>
void gsExprAssembler<T>::computeDirichletDofs2(int unk)
{
    expr::gsFeSpace<T> & u = *m_vcol[unk];
    
    if ( m_options.getInt("DirichletStrategy") == dirichlet::nitsche)
        return; // Nothing to compute

    const gsBoundaryConditions<> & bbc = u.hasBc() ? u.bc() : gsBoundaryConditions<>();
    
    const gsMultiBasis<T> & mbasis = dynamic_cast<const gsMultiBasis<T>&>(u.source());

    // eg. not penalize
    const gsDofMapper & mapper = 
            dirichlet::elimination == m_options.getInt("DirichletStrategy") ?
        u.mapper():
        mbasis.getMapper(dirichlet::elimination,
                         static_cast<iFace::strategy>(m_options.getInt("InterfaceStrategy")),
                         bbc, u.id());

    switch ( m_options.getInt("DirichletValues") )
    {
    case dirichlet::homogeneous:
        // If we have a homogeneous Dirichlet problem fill boundary
        // DoFs with zeros
        u.fixedPart().setZero(mapper.boundarySize(), u.dim() );
        break;
    case dirichlet::interpolation:
        computeDirichletDofsIntpl2(u);
        break;
    case dirichlet::l2Projection:
        computeDirichletDofsL2Proj(u); //this->computeDirichletDofsL2Proj(mapper, mbasis,unk);
        break;
    case dirichlet::user :
        // Assuming that the DoFs are already set by the user
        GISMO_ENSURE( u.fixedPart().rows() == mapper.boundarySize() &&
                      u.fixedPart().cols() == u.dim(),
                      "The Dirichlet DoFs are not set.");
        break;
    default:
        GISMO_ERROR("Something went wrong with Dirichlet values.");
    }
    
    // Corner values
    for ( typename gsBoundaryConditions<T>::const_citerator
              it = bbc.cornerBegin();
          it != bbc.cornerEnd(); ++it )
    {
        if(it->unknown == unk)
        {
            const int i  = mbasis[it->patch].functionAtCorner(it->corner);
            const int ii = mapper.bindex( i , it->patch );
            u.fixedPart().row(ii).setConstant(it->value);
        }
        else
            continue;
    }
}

template<class T>
void gsExprAssembler<T>::setFixedDofVector(gsMatrix<T> & vals, int unk)
{
	expr::gsFeSpace<T> & u = *m_vcol[unk];
	gsMatrix<T>        & fixedDofs = const_cast<expr::gsFeSpace<T>&>(u).fixedPart();
	fixedDofs.swap(vals);
	vals.resize(0, 0);
	// Assuming that the DoFs are already set by the user
	GISMO_ENSURE( fixedDofs.rows() == u.mapper().boundarySize() && 
				  fixedDofs.cols() == u.dim(),
   				  "The Dirichlet DoFs were not provided correctly.");
}

template<class T>
void gsExprAssembler<T>::setFixedDofs(const gsMatrix<T> & coefMatrix, int unk, int patch)
{
    GISMO_ASSERT( m_options.getInt("DirichletValues") == dirichlet::user, "Incorrect options");

	expr::gsFeSpace<T> & u = *m_vcol[unk];
    const int dirStr = m_options.getInt("DirichletStrategy");
    const gsMultiBasis<T> & mbasis = *dynamic_cast<const gsMultiBasis<T>* >(&(u).source());

	const gsBoundaryConditions<> & bbc = u.hasBc() ? u.bc() : gsBoundaryConditions<>();
	
    const gsDofMapper & mapper =
        dirichlet::elimination == dirStr ? u.mapper()
        : mbasis.getMapper(dirichlet::elimination,
                           static_cast<iFace::strategy>(m_options.getInt("InterfaceStrategy")),
                           bbc, u.id()) ;

	gsMatrix<T> & fixedDofs = const_cast<expr::gsFeSpace<T>& >(u).fixedPart();
	GISMO_ASSERT(fixedDofs.rows() == mapper.boundarySize() &&
				 fixedDofs.cols() == u.dim(),
				 "Fixed DoFs were not initialized.");

    // for every side with a Dirichlet BC
    for ( typename gsBoundaryConditions<T>::const_iterator
          it =  bbc.dirichletBegin();
          it != bbc.dirichletEnd()  ; ++it )
    {
        const int k = it->patch();
        if ( k == patch )
        {
            // Get indices in the patch on this boundary
            const gsMatrix<unsigned> boundary =
                    mbasis[k].boundary(it->side());
            
            //gsInfo <<"Setting the value for: "<< boundary.transpose() <<"\n";
            
            for (index_t i=0; i!= boundary.size(); ++i)
            {
                // Note: boundary.at(i) is the patch-local index of a
                // control point on the patch
                const int ii  = mapper.bindex( boundary.at(i) , k );
                
                fixedDofs.row(ii) = coefMatrix.row(boundary.at(i));
            }
        }
    }
} // setFixedDofs 

template<class T>
template<bool left, bool right, class E1, class E2>
void gsExprAssembler<T>::assembleLhsRhs_impl(const expr::_expr<E1> & exprLhs,
                                             const expr::_expr<E2> & exprRhs,
                                             space rvar, space cvar)
{
    //GISMO_ASSERT( exprLhs.isMatrix(), "Expecting matrix expression"); // or null
    //GISMO_ASSERT( exprRhs.isVector(), "Expecting vector expression");

    /*
    exprLhs.print(gsInfo);
    if (right) exprRhs.print(gsInfo);
    */

    GISMO_ASSERT(matrix().cols()==numDofs(), "System not initialized");
    
    // initialize flags
    m_exprdata->initFlags(SAME_ELEMENT|NEED_ACTIVE, SAME_ELEMENT);
    if (left ) exprLhs.setFlag();
    if (right) exprRhs.setFlag();
    
    gsMatrix<T> localMat, localRhs; // Local matrix and local rhs        
    gsQuadRule<T> QuRule;  // Quadrature rule
    gsVector<T> quWeights; // quadrature weights
    
    for (unsigned patchInd = 0; patchInd < m_exprdata->multiBasis().nBases(); ++patchInd) //**1
    {
        QuRule = gsGaussRule<T>(m_exprdata->multiBasis().basis(patchInd), m_options);
		
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
                
            // Assemble on element
            if (left ) localMat = quWeights[0] * exprLhs.eval(0);
            if (right) localRhs = quWeights[0] * exprRhs.eval(0);
            for (index_t k = 1; k != quWeights.rows(); ++k)
            {
                if (left ) localMat += quWeights[k] * exprLhs.eval(k);
                if (right) localRhs += quWeights[k] * exprRhs.eval(k);
            }

            // gsDebugVar(localMat.dim());
            //gsDebugVar(localRhs.transpose() );
            
            // Add contributions to the system matrix and right-hand side
            gsSpAcc<T,left,right>(m_matrix, m_rhs,
                                  localMat, localRhs, rvar, cvar, patchInd);
        }
    }
    
    //this->finalize();
    m_matrix.makeCompressed();
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
    
    // Local matrix and local rhs
    gsMatrix<T> localMat, localRhs;
    gsVector<T> quWeights;// quadrature weights
    gsQuadRule<T>  QuRule;

    for (typename bcContainer::const_iterator it = BCs.begin(); it!= BCs.end(); ++it)
    {
        QuRule = gsGaussRule<T>(m_exprdata->multiBasis().basis(it->patch()), m_options,
                                it->side().direction());

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
            
            if (left)  localMat = quWeights[0] * exprLhs.eval(0);
            if (right) localRhs = quWeights[0] * exprRhs.eval(0);
            for (index_t k = 1; k != quWeights.rows(); ++k)
            {
                if (left ) localMat += quWeights[k] * exprLhs.eval(k);
                if (right) localRhs += quWeights[k] * exprRhs.eval(k);
            }
            
            // Add contributions to the system matrix and right-hand side
            gsSpAcc<T,left,right>(m_matrix, m_rhs,
                                  localMat, localRhs, rvar, cvar, it->patch() );
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
                                                const intContainer & iFaces)
{
    //GISMO_ASSERT( exprRhs.isVector(), "Expecting vector expression");
    
    // initialize flags
    m_exprdata->initFlags(SAME_ELEMENT|NEED_ACTIVE, SAME_ELEMENT);
    m_exprdata->parse(exprLhs,exprRhs);
    
    // Local matrix
    gsMatrix<T> localMat, localRhs;
    gsVector<T> quWeights;// quadrature weights
    gsQuadRule<T>  QuRule;
    
    gsMatrix<T> tmp;

    for (gsBoxTopology::const_iiterator it = iFaces.begin(); 
         it != iFaces.end(); ++it )
    {
        const boundaryInterface & iFace = *it;
        const int patch1 = iFace.first() .patch;
        const int patch2 = iFace.second().patch;
        //const gsAffineFunction<T> interfaceMap(m_pde_ptr->patches().getMapForInterface(bi));

        QuRule = gsGaussRule<T>(m_exprdata->multiBasis().basis(patch1), m_options,
                                iFace.first().side().direction());

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

            if (left ) localMat = quWeights[0] * exprLhs.eval(0);
            if (right) localRhs = quWeights[0] * exprRhs.eval(0);
            for (index_t k = 1; k != quWeights.rows(); ++k)
            {
                if (left ) localMat += quWeights[k] * exprLhs.eval(k);
                if (right) localRhs += quWeights[k] * exprRhs.eval(k);
            }
            
            // Add contributions to the system matrix and right-hand side
            gsSpAcc<T,left,right>(m_matrix, m_rhs,
                                  localMat, localRhs, rvar, cvar, patch1);
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
                    
    const gsBoundaryConditions<> & bbc = u.hasBc() ? u.bc() : gsBoundaryConditions<>();
    
    // Iterate over all patch-sides with Dirichlet-boundary conditions
    for ( typename gsBoundaryConditions<T>::const_iterator
          it = bbc.dirichletBegin(); it != bbc.dirichletEnd(); ++it )
    {
        const int k = it->patch();
        if( it->unknown()!=u.id() )
            continue;
        const gsBasis<T> & basis = mbasis[k];

        // Get dofs on this boundary
        const gsMatrix<unsigned> boundary = basis.boundary(it->side());

        // If the condition is homogeneous then fill with zeros
        if ( it->isHomogeneous() )
        {
            for (index_t i=0; i!= boundary.size(); ++i)
            {
                const int ii= mapper.bindex( boundary.at(i) , k );
                fixedDofs.row(ii).setZero();
            }
            continue;
        }

        // Get the side information
        int dir = it->side().direction( );
        index_t param = (it->side().parameter() ? 1 : 0);

        // Compute grid of points on the face ("face anchors")
        std::vector< gsVector<T> > rr;
        rr.reserve( parDim );

        for ( int i=0; i < parDim; ++i)
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


        if ( fpts.rows() != u.dim() )
        {
            // assume scalar 
            gsMatrix<T> tmp(u.dim(), fpts.cols());
            tmp.setZero();
            gsDebugVar(!dir);
            tmp.row(!dir) = (param ? 1 : -1) * fpts; // normal !
            fpts.swap(tmp);
        }
        
        // Interpolate dirichlet boundary
        typename gsBasis<T>::uPtr h = basis.boundaryBasis(it->side());
        typename gsGeometry<T>::uPtr geo = h->interpolateAtAnchors(fpts);
        const gsMatrix<T> & dVals =  geo->coefs();

        // Save corresponding boundary dofs
        for (index_t l=0; l!= boundary.size(); ++l)
        {
            const int ii = mapper.bindex( boundary.at(l) , it->patch() );

            fixedDofs.row(ii) = dVals.row(l);
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
    gsMatrix<unsigned> act;
    gsSparseMatrix<T> cmat(bsz, bsz);
    // todo: reserve
    
    const gsMultiBasis<T> & mbasis =
        *dynamic_cast<const gsMultiBasis<T>*>(&u.source());
                    
    const gsBoundaryConditions<> & bbc = u.hasBc() ? u.bc() : gsBoundaryConditions<>();
    
    // Iterate over all patch-sides with Dirichlet-boundary conditions
    for ( typename gsBoundaryConditions<T>::const_iterator
          it = bbc.dirichletBegin(); it != bbc.dirichletEnd(); ++it )
    {
        GISMO_ASSERT(it->function()->targetDim() == u.dim(),
                     "Given Dirichlet boundary function does not match problem dimension."
                     <<it->function()->targetDim()<<" != "<<u.dim()<<"\n");

        const int k   = it->patch();
        if( it->unknown()!=u.id() )
            continue;
        const gsBasis<T> & basis = mbasis[k];

        // Get dofs on this boundary
        gsMatrix<unsigned> boundary = basis.boundary(it->side());

        // If the condition is homogeneous then fill with zeros
        if ( it->isHomogeneous() )
        {
            for (index_t i=0; i!= boundary.size(); ++i)
            {
                const int ii= mapper.bindex( boundary.at(i) , k );
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

	const gsBoundaryConditions<> & bbc = u.hasBc() ? u.bc() : gsBoundaryConditions<>(); 

    // Set up matrix, right-hand-side and solution vector/matrix for
    // the L2-projection
    gsSparseEntries<T> projMatEntries;
    gsMatrix<T>        globProjRhs;
    globProjRhs.setZero( mapper.boundarySize(), u.dim() );

    // Temporaries
    gsMatrix<T> quNodes;
    gsVector<T> quWeights;

    gsMatrix<T> rhsVals;
    gsMatrix<unsigned> globIdxAct;
    gsMatrix<T> basisVals;

	const gsMultiPatch<T> * mp = (const gsMultiPatch<T>*)(& m_exprdata->getMap().source());	
	
    // Iterate over all patch-sides with Dirichlet-boundary conditions
    for ( typename gsBoundaryConditions<T>::const_iterator
          iter = bbc.dirichletBegin();
          iter != bbc.dirichletEnd(); ++iter )
    {
        const int unk = iter->unknown();
        if(unk != u.id())
            continue;
        const int patchIdx   = iter->patch();
        const gsBasis<T> & basis = mbasis[patchIdx];

        typename gsGeometry<T>::Evaluator geoEval( mp->patches()[patchIdx]->evaluator(NEED_MEASURE));

        // Set up quadrature to degree+1 Gauss points per direction,
        // all lying on iter->side() except from the direction which
        // is NOT along the element
        gsGaussRule<T> bdQuRule(basis, 1.0, 1, iter->side().direction());

        // Create the iterator along the given part boundary.
        typename gsBasis<T>::domainIter bdryIter = basis.makeDomainIterator(iter->side());

        for(; bdryIter->good(); bdryIter->next() )
        {
            bdQuRule.mapTo( bdryIter->lowerCorner(), bdryIter->upperCorner(),
                            quNodes, quWeights);

            geoEval->evaluateAt( quNodes );

            // the values of the boundary condition are stored
            // to rhsVals. Here, "rhs" refers to the right-hand-side
            // of the L2-projection, not of the PDE.
            rhsVals = iter->function()->eval( m_exprdata->getMap().source().piece(patchIdx).eval( quNodes ) );

            basis.eval_into( quNodes, basisVals);

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
            basis.active_into(quNodes.col(0), globIdxAct );
            mapper.localToGlobal( globIdxAct, patchIdx, globIdxAct);

            // Out of the active functions/DOFs on this element, collect all those
            // which correspond to a boundary DOF.
            // This is checked by calling mapper.is_boundary_index( global Index )

            // eltBdryFcts stores the row in basisVals/globIdxAct, i.e.,
            // something like a "element-wise index"
            std::vector<index_t> eltBdryFcts;
            eltBdryFcts.reserve(mapper.boundarySize());
            for( index_t i=0; i < globIdxAct.rows(); i++)
                if( mapper.is_boundary_index( globIdxAct(i,0)) )
                    eltBdryFcts.push_back( i );

            // Do the actual assembly:
            for( index_t k=0; k < quNodes.cols(); k++ )
            {
                const T weight_k = quWeights[k] * geoEval->measure(k);
				
                // Only run through the active boundary functions on the element:
                for( size_t i0=0; i0 < eltBdryFcts.size(); i0++ )
                {
                    // Each active boundary function/DOF in eltBdryFcts has...
                    // ...the above-mentioned "element-wise index"
                    const unsigned i = eltBdryFcts[i0];
                    // ...the boundary index.
                    const unsigned ii = mapper.global_to_bindex( globIdxAct( i ));

                    for( size_t j0=0; j0 < eltBdryFcts.size(); j0++ )
                    {
                        const unsigned j = eltBdryFcts[j0];
                        const unsigned jj = mapper.global_to_bindex( globIdxAct( j ));

                        // Use the "element-wise index" to get the needed
                        // function value.
                        // Use the boundary index to put the value in the proper
                        // place in the global projection matrix.
                        projMatEntries.add(ii, jj, weight_k * basisVals(i,k) * basisVals(j,k));
                    } // for j

                    globProjRhs.row(ii) += weight_k *  basisVals(i,k) * rhsVals.col(k).transpose();

                } // for i
            } // for k
        } // bdryIter
    } // boundaryConditions-Iterator

    gsSparseMatrix<T> globProjMat( mapper.boundarySize(), mapper.boundarySize() );
    globProjMat.setFrom( projMatEntries );
    globProjMat.makeCompressed();

    // Solve the linear system:
    // The position in the solution vector already corresponds to the
    // numbering by the boundary index. Hence, we can simply take them
    // for the values of the eliminated Dirichlet DOFs.
    typename gsSparseSolver<T>::CGDiagonal solver;
    fixedDofs = solver.compute( globProjMat ).solve ( globProjRhs );
    
} // computeDirichletDofsL2Proj

    
} //namespace gismo
