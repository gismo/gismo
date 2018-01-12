/** @file gsExpressions.h

    @brief Defines different expressions 

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris
*/

#pragma once

#include <gsCore/gsFuncData.h>
#include <gsUtils/gsSortedVector.h>

namespace gismo
{

template<class T> class gsExprHelper;

// Adaptor to compute Hessian
template<typename T>
void secDerToHessian(typename gsMatrix<T>::constRef & secDers,
                     const index_t dim,
                     Eigen::Matrix<T,-1,-1> & hessian)
{
    const index_t sz = dim*(dim+1)/2;
    const gsAsConstMatrix<T> ders(secDers.data(), sz, secDers.size() / sz );
    hessian.resize(dim*dim, ders.cols() );

    switch ( dim )
    {
    case 1:
        hessian = secDers;
        break;
    case 2:
        hessian.row(0)=ders.row(0);//0,0
        hessian.row(1)=//1,0
        hessian.row(2)=ders.row(2);//0,1
        hessian.row(3)=ders.row(1);//1,1
        break;
    case 3:
        hessian.row(0)=ders.row(0);//0,0
        hessian.row(3)=//0,1
        hessian.row(1)=ders.row(3);//1,0
        hessian.row(6)=//0,2
        hessian.row(2)=ders.row(4);//2,0
        hessian.row(4)=ders.row(1);//1,1
        hessian.row(7)=//1,2
        hessian.row(5)=ders.row(5);//2,1
        hessian.row(8)=ders.row(2);//2,2
        break;
    default:
        break;
    }
}

namespace expr
{

#if(__cplusplus >= 201402L) // c++14
#  define MatExprType  auto
#  define AutoReturn_t auto
//#elif(__cplusplus >= 201103L) // c++11
//note: in c++11 auto-return requires -> decltype(.)
#else // 199711L
#  define MatExprType typename gsMatrix<Scalar>::constRef
#  define AutoReturn_t typename util::conditional<ScalarValued,Scalar,MatExprType>::type
#endif

/**
   Traits class for expressions
 */
template <typename E> struct expr_traits
{
    typedef real_t Scalar;//todo
    typedef const E Nested_t;
    //typedef E type;
};

template<class E> class symm_expr;
template<class E> class trace_expr;
template<class E> class norm_expr;
template<class E> class sqNorm_expr;
template<class E> class value_expr;
template<class E> class asdiag_expr;
template<class T> class meas_expr;
template<class E> class inv_expr;
template<class E> class tr_expr;
template<class E> class temp_expr;
//bool = E1::ScalarValued && E2::ScalarValued
template<class E1, class E2, bool = E1::ColBlocks> class mult_expr
{using E1::GISMO_ERROR_mult_expr_has_invalid_template_arguments;};

template<class T> class gsFeVariable;
template<class T> class gsFeSolution;

/**
   Base class for all expressions
 */
template <typename E>
class _expr
{
private:
     _expr(){}
     _expr(const _expr&) { }
    friend E;
public:
    
    enum {ScalarValued = 0, ColBlocks = 0}; // ValueType=0,1,2 (scalar,vector,matrix)
    
    typedef typename expr_traits<E>::Nested_t Nested_t;
    typedef typename expr_traits<E>::Scalar   Scalar;

    void print(std::ostream &os) const
    { static_cast<E const&>(*this).print(os); os<<"\n"; }
    /*
    {
        std::string tmp(__PRETTY_FUNCTION__);
        tmp.erase(0,74);
        tmp.erase(tmp.size()-42,42);
        size_t pos = 0;
        while((pos=tmp.find(", false",0))!=std::string::npos) tmp.erase(pos,7);
        while((pos=tmp.find(", true",0))!=std::string::npos) tmp.erase(pos,6);
        while((pos=tmp.find("gismo::expr::",0))!=std::string::npos) tmp.erase(pos,13);
        while((pos=tmp.find("_expr",0))!=std::string::npos) tmp.erase(pos,5);
        while((pos=tmp.find("<double>",0))!=std::string::npos) tmp.erase(pos,8);
        // while((pos=tmp.find("<long double>",0))!=std::string::npos) tmp.erase(pos,13);
        // while((pos=tmp.find("<float>",0))!=std::string::npos) tmp.erase(pos,7);
        tmp.erase(std::remove_if(tmp.begin(),tmp.end(),::isspace),tmp.end());
        os<<tmp<<"\n";
    }
    //*/
    
    std::ostream & printDetail(std::ostream &os) const
    {
        os << (isScalar() ? "Scalar " :
               (isVector() ? "Vector " :
                (isMatrix() ? "Matrix " :
                 "Unknown ") ) )
           <<"expression of size "<< rows() // big: this is not fixed, or may not be known
           << " x "<<cols()<<"\n";
        print(os);
        return os;
    }

    MatExprType eval(const index_t k) const
    { return static_cast<E const&>(*this).eval(k); }

    tr_expr<E> tr() const
    { return tr_expr<E>(static_cast<E const&>(*this)); }

    temp_expr<E> temp() const
    { return temp_expr<E>(static_cast<E const&>(*this)); }

    inv_expr<E> const inv() const
    { return inv_expr<E>(static_cast<E const&>(*this)); }

    trace_expr<E> trace() const
    { return trace_expr<E>(static_cast<E const&>(*this)); }

    norm_expr<E> norm() const
    { return norm_expr<E>(static_cast<E const&>(*this)); }

    sqNorm_expr<E> sqNorm() const
    { return sqNorm_expr<E>(static_cast<E const&>(*this)); }

    mult_expr<E,E,0> sqr() const { return (*this)*(*this); }
    
    symm_expr<E> symm() const
    { return symm_expr<E>(static_cast<E const&>(*this)); }

    value_expr<E> val() const
    { return value_expr<E>(static_cast<E const&>(*this)); }

    asdiag_expr<E> asDiag() const
    { return asdiag_expr<E>(static_cast<E const&>(*this)); }

    index_t rows() const
    { return static_cast<E const&>(*this).rows(); }
    
    index_t cols() const
    { return static_cast<E const&>(*this).cols(); }

    // note: runtime check. For compile-time check use E::ScalarValued
    bool isScalar() const { return rows()*cols()<=1; } //!rowSpan() && !colSpan()

    bool rowSpan() const // remove
    { return static_cast<E const&>(*this).rowSpan(); }
    
    bool colSpan() const // remove
    { return static_cast<E const&>(*this).colSpan(); }

    bool isVector() const { return rowSpan() && (!colSpan()); }

    bool isMatrix() const { return rowSpan() && colSpan(); }
    
    void setFlag() const
    { static_cast<E const&>(*this).setFlag(); }

    void parse(gsSortedVector<const gsFunctionSet<Scalar>*> & evList) const
    { static_cast<E const&>(*this).parse(evList); }

    const gsFeVariable<Scalar> & rowVar() const
    {
        // assert ValueType!=0
        return static_cast<E const&>(*this).rowVar();
    }
    
    const gsFeVariable<Scalar> & colVar() const
    {
        // assert ValueType==2
        return static_cast<E const&>(*this).colVar();
    }

    //E const & derived() const { return static_cast<const E&>(*this); }

    // Overload conversions to E, e.g., for _expr<mult_expr>, converts
    // to mult_expr.
    operator E&()             { return static_cast<      E&>(*this); }
    operator E const&() const { return static_cast<const E&>(*this); }
};

template <typename E>
std::ostream &operator<<(std::ostream &os, const _expr<E> & b)
{b.print(os); return os; }

/** 
    Expression for a constant value
 */
template <>
class _expr<real_t> : public _expr<_expr<real_t> >
{
    const real_t _c;
public:
    typedef real_t Scalar;
    typedef const _expr<real_t> Nested_t;
    
    _expr(const Scalar & c) : _c(c) { }
    
public:
    enum {ScalarValued = 1};
    
    inline Scalar eval(const index_t k) const { return _c; }

    inline _expr val() const { return *this; }
    index_t rows() const { return 0; }
    index_t cols() const { return 0; }
    void setFlag() const { }
    void parse(gsSortedVector<const gsFunctionSet<Scalar>*> & evList) const { GISMO_UNUSED(evList); }
    bool rowSpan() const {return false;}
    bool colSpan() const {return false;}
    const gsFeVariable<real_t> & rowVar() const { GISMO_ERROR("scalar: rowVar"); }
    const gsFeVariable<real_t> & colVar() const { GISMO_ERROR("scalar: colVar"); }

    void print(std::ostream &os) const { os<<_c; }
};

/// Returns the unit as an expression
EIGEN_STRONG_INLINE _expr<real_t> one() { return _expr<real_t>(1); }

/**
   Geometry map expression
 */
template<class T>
class gsGeometryMap : public _expr<gsGeometryMap<T> >
{
    const gsFunctionSet<T> * m_fs;
    const gsMapData<T> *  m_fd;

    //index_t d, n;
    
public:
    typedef T Scalar;

    friend class gismo::gsExprHelper<T>;

    const gsFunctionSet<T> & source() const {return *m_fs;}

    const gsMapData<T> & data() const  { return *m_fd; }

    void print(std::ostream &os) const { os << "G"; }

    // was protected
    MatExprType eval(const index_t k) const { return m_fd->values[0].col(k); }

    // was protected
    void setFlag() const
    {
        GISMO_ASSERT(NULL!=m_fd, "GeometryMap not registered");
        m_fd->flags |= NEED_VALUE;
    }

protected:
    
    gsGeometryMap() : m_fs(NULL), m_fd(NULL) { }

    void registerData(const gsFunctionSet<T> & fs, const gsMapData<T> & val)
    {
        m_fs = &fs;
        m_fd = &val;
    }

    bool isValid() const { return NULL!=m_fs; }
    
    index_t rows() const { return m_fd->dim.second; }
    index_t cols() const { return 1; }
    
    bool rowSpan() const {return false;}
    bool colSpan() const {return false;}
    
    const gsFeVariable<T> & rowVar() const { GISMO_ERROR("GeometryMap: rowVar"); }
    const gsFeVariable<T> & colVar() const { GISMO_ERROR("GeometryMap: colVar"); }

    void parse(gsSortedVector<const gsFunctionExpr<T>*> & evList) const 
    { 
        GISMO_ASSERT(NULL!=m_fd, "GeometryMap not registered");
        evList.push_unique(m_fs);
        m_fd->flags |= NEED_VALUE;
    }
};

/**
   Null expression is a compatibility expression invalid at runtime
 */
template<class T = real_t>
class gsNullExpr : public _expr<gsNullExpr<T> >
{
public:
    typedef T Scalar;
    gsMatrix<T> eval(const index_t k) const { GISMO_ERROR("gsNullExpr"); }
    inline index_t rows() const { GISMO_ERROR("gsNullExpr"); }    
    inline index_t cols() const { GISMO_ERROR("gsNullExpr"); }
    inline void setFlag() const { }
    void parse(gsSortedVector<const gsFunctionSet<T>*> & evList) const { GISMO_UNUSED(evList); }

    const gsFeVariable<T> & rowVar() const { GISMO_ERROR("gsNullExpr"); }
    const gsFeVariable<T> & colVar() const { GISMO_ERROR("gsNullExpr"); }
    bool rowSpan() const {return false; }
    bool colSpan() const {return false;}

    void print(std::ostream &os) const { os << "NullExpr"; }
};

template<class T> class cdiam_expr;

/**
   An element object collecting relevant expressions
*/
template<class T>
class gsFeElement
{
    friend class cdiam_expr<T>;
    
    const gsDomainIterator<T> * m_di;

    cdiam_expr<T> cd;
public:    
    typedef T Scalar;
    
    gsFeElement() : m_di(NULL), cd(*this) { }
    
    void set(const gsDomainIterator<T> & di)
    { m_di = &di; }

    const cdiam_expr<T> & diam() const
    {
        return cd;
        //return cdiam_expr<T>(*this);
        // return cdiam_expr<T>(m_di); // bug ..
    }

    void print(std::ostream &os) const { os << "e"; }
};

/**
   An expression of the element diameter
 */
template<class T>
class cdiam_expr : public _expr<cdiam_expr<T> >
{
    const gsFeElement<T> & _e;
public:
    typedef T Scalar;

    enum {ScalarValued = 1};
    
    explicit cdiam_expr(const gsFeElement<T> & el) : _e(el) { }
    
    T eval(const index_t k) const { return _e.m_di->getCellSize(); }

    inline cdiam_expr<T> val() const { return *this; }
    inline index_t rows() const { return 0; }    
    inline index_t cols() const { return 0; }
    inline void setFlag() const { }
    void parse(gsSortedVector<const gsFunctionSet<Scalar>*> & evList) const { GISMO_UNUSED(evList); }

    const gsFeVariable<T> & rowVar() const { GISMO_ERROR("diam problem rowVar"); }
    const gsFeVariable<T> & colVar() const { GISMO_ERROR("diam problem colVar"); }
    bool rowSpan() const {return false; }
    bool colSpan() const {return false;}

    void print(std::ostream &os) const
    { os << "diam(e)"; }
};

/*
template<class T>
class parNv_expr : public _expr<parNv_expr<T> >
{
    const gsFeElement<T> & _e;
public:
    typedef T Scalar;

    enum {ScalarValued = 1};
    
    explicit parNv_expr(const gsFeElement<T> & el) : _e(el) { }
    
    T eval(const index_t k) const { return _e.m_di->getCellSize(); }

    inline cdiam_expr<T> val() const { return *this; }
    inline index_t rows() const { return 0; }    
    inline index_t cols() const { return 0; }
    inline void setFlag() const { }
    void parse(gsSortedVector<const gsFunctionSet<Scalar>*> & evList) const { GISMO_UNUSED(evList); }
    const gsFeVariable<T> & rowVar() const { GISMO_ERROR("diam problem rowVar"); }
    const gsFeVariable<T> & colVar() const { GISMO_ERROR("diam problem colVar"); }
    bool rowSpan() const {return false; }
    bool colSpan() const {return false;}

    void print(std::ostream &os) const
    { os << "diam(e)"; }
};
*/

template <typename T> 
struct expr_traits<gsGeometryMap<T> >
{
    typedef T Scalar;
    typedef gsGeometryMap<T> const & Nested_t;
};

/**
   Expression for a finite element variable
*/
template<class T>
class gsFeVariable  : public _expr<gsFeVariable<T> >
{
protected:
    const gsFunctionSet<T> * m_fs;
    const gsFuncData<T>    * m_fd;
    index_t m_d;

public:
    typedef T Scalar;

    friend class gismo::gsExprHelper<T>;

    const gsFuncData<T> & data() const {return *m_fd;}
    
    const gsFunctionSet<T> & source() const {return *m_fs;}
    void setSource(const gsFunctionSet<T> & fs) { m_fs = &fs;}
    
    void clear() { m_fs = NULL; }
    
protected:
    
    explicit gsFeVariable(index_t _d = 1) : m_fs(NULL), m_fd(NULL), m_d(_d)
    { }
    
    void registerData(const gsFunctionSet<T> & fs, const gsFuncData<T> & val, index_t d)
    {
        GISMO_ASSERT(NULL==m_fs, "gsFeVariable: already registered");
        m_fs = &fs ;
        m_fd = &val;
        m_d  = d;
    }

    void setData(const gsFuncData<T> & val) { m_fd = &val;}

    gsFuncData<T> & data() {return *m_fd;}

    bool isValid() const { return NULL!=m_fd && NULL!=m_fs; }

public:

    // component
    // expr comp(const index_t i) const { return comp_expr<T>(*this,i); }
    // eval(k).col(i)

    // rows: functions. cols: components
    MatExprType eval(const index_t k) const
    { return m_fd->values[0].col(k).blockDiag(m_d); }
    //{ return m_fd->values[0].block(0,k,rows(),1); }

    const gsFeVariable<T> & rowVar() const {return *this;}
    const gsFeVariable<T> & colVar() const
    {GISMO_ERROR("Invalid call to gsFeVariable::colVar()");}

    bool rowSpan() const {return true; } // note: gsFunction ??
    bool colSpan() const {return false;}
    
    index_t rows() const
    {
        GISMO_ASSERT(NULL!=m_fs, "FeVariable: Function member not registered");
        GISMO_ASSERT(NULL!=m_fd, "FeVariable: FuncData member not registered");

        //gsDebugVar(*m_fd);
        /*
        gsDebugVar(m_fd->flags & NEED_VALUE);
        gsDebugVar(m_fd->values[0].rows());

        gsDebugVar(m_fd->flags & NEED_ACTIVE);
        gsDebugVar(m_fd->actives.rows());

        gsDebugVar(m_fd->flags & NEED_DERIV);
        */
        
        // note: precomputation is needed
        if (m_fd->flags & NEED_VALUE)
        {return m_d * m_fd->values[0].rows();}
        if (m_fd->flags & NEED_ACTIVE) // coeffs ??
        {return m_d * m_fd->actives.rows();}
        if (m_fd->flags & NEED_DERIV)
        {return m_d * m_fd->values[0].rows();}
        GISMO_ERROR("Cannot deduce row size.");
    }
    
    index_t cols() const { return m_d; }
    
    void setFlag() const
    {
        GISMO_ASSERT(NULL!=m_fd, "FeVariable: FuncData member not registered");
        m_fd->flags |= NEED_VALUE;
    }

    void parse(gsSortedVector<const gsFunctionSet<Scalar>*> & evList) const
    {
        GISMO_ASSERT(NULL!=m_fd, "FeVariable: FuncData member not registered");
        evList.push_sorted_unique(m_fs);
        m_fd->flags |= NEED_VALUE;
    }
    
    index_t dim() const { return m_d;}

    index_t parDim() const
    {
        return m_fs->domainDim();
        //return m_fd->dim.first;
    }
    index_t cSize()  const { return m_fd->values[0].rows(); } // coordinate size

    void print(std::ostream &os) const { os << "u"; }
};

template<class T>
class gsFeSpace :public gsFeVariable<T>
{
protected:
    friend gsFeSolution<T>;

    typedef gsFeVariable<T> Base;

    index_t             m_id;

    // C^r coupling
    mutable index_t m_r; // iFace::handling

    typedef typename gsBoundaryConditions<T>::bcRefList bcRefList;
    mutable bcRefList m_bcs;
    
    gsDofMapper m_mapper;
    gsMatrix<T> m_fixedDofs;

public:
    typedef T Scalar;

    gsDofMapper & mapper() {return m_mapper;}
    const gsDofMapper & mapper() const {return m_mapper;}

    inline const gsMatrix<T> & fixedPart() const {return m_fixedDofs;}
    gsMatrix<T> & fixedPart() {return m_fixedDofs;}

    const bcRefList & bc() const { return m_bcs; }
    void addBc(bcRefList bc) const { m_bcs = bc; }
    void clearBc() const { m_bcs.clear(); }
    std::size_t bcSize() const { return m_bcs.size(); }

    index_t   id() const {return m_id;}
    index_t & setId(const index_t _id) {return m_id = _id;}

    index_t   interfaceCont() const {return m_r;}
    index_t & setInterfaceCont(const index_t _r) const
    {
        GISMO_ASSERT(_r>-2 && _r<1, "Invalid or not implemented (r="<<_r<<").");
        return m_r = _r;
    }

    //void getFunction(const gsMatrix<T>& solVector, gsMultiPatch<T>& result);
    //gsAsFunction<T> asFunction (const gsMatrix<T>& solVector);

    void getCoeffs(const gsMatrix<T>& solVector, gsMatrix<T> & result,
                   const index_t p = 0) const
    {
        // GISMO_ASSERT(p<
        const index_t dim = this->dim();

        const gsMultiBasis<T> & mb = static_cast<const gsMultiBasis<T>&>(this->source());
        GISMO_ASSERT( dynamic_cast<const gsMultiBasis<T>*>(&this->source()), "error");

        // Reconstruct solution coefficients on patch p
        const int sz  = mb[p].size();
        result.resize(sz, dim); // (!)
        
        for (index_t i = 0; i < sz; ++i)
        {
            const int ii = m_mapper.index(i, p);
            
            if ( m_mapper.is_free_index(ii) ) // DoF value is in the solVector
            {
                for (index_t k = 0; k != dim; ++k)
                {
                    const index_t cgs = k * m_mapper.freeSize(); //global stride
                    result(i,k) = solVector.at( cgs + ii );
                }
            }
            else // eliminated DoF: fill with Dirichlet data
            {
                result.row(i) = m_fixedDofs.row( m_mapper.global_to_bindex(ii) ).head(dim);
            }
        }
    }

    // space restrictTo(boundaries);
    // space restrictTo(bcRefList domain);

    void reset()
    {
        if (const gsMultiBasis<T> * mb =
            dynamic_cast<const gsMultiBasis<T>*>(&this->source()) )
        {
            m_mapper = gsDofMapper(*mb);
            //m_mapper.init(*mb); //bug
            if ( 0==this->interfaceCont() ) // Conforming boundaries ?
            {
                for ( gsBoxTopology::const_iiterator it = mb->topology().iBegin();
                      it != mb->topology().iEnd(); ++it )
                {
                    mb->matchInterface(*it, m_mapper);
                }
            }

            gsMatrix<unsigned> bnd;
            for (typename bcRefList::const_iterator
                     it = this->bc().begin() ; it != this->bc().end(); ++it )
            {
                GISMO_ASSERT( it->get().ps.patch < static_cast<index_t>(this->mapper().numPatches()),
                              "Problem: a boundary condition is set on a patch id which does not exist.");

                bnd = mb->basis(it->get().ps.patch).boundary( it->get().ps.side() );
                m_mapper.markBoundary(it->get().ps.patch, bnd);
            }
        }
/*
        else if (const gsBasis<T> * b =
                 dynamic_cast<const gsBasis<T>*>(&u.source()) )
        {
            gsMultiBasis<T> mbb(*b);
            mbb.getMapper(
                dirichlet::elimination,
                0==u.interfaceCont() ? iFace::conforming : iFace::none,
                ubc, u.mapper(), u.id(), true);
        }
        else
        {
            gsWarn<<"Problem initializing.\n";
        }
*/
        this->mapper().finalize();
        //this->mapper().print();
    }

protected:
    friend class gismo::gsExprHelper<T>;    
    explicit gsFeSpace(index_t _d = 1) : Base(_d), m_id(-1), m_r(-1)
    { }
};


template<class T>
class gsFeSolution :public _expr<gsFeSolution<T> >
{
protected:
    const gsFeSpace<T> & _u;

    gsMatrix<T> * _Sv;
public:
    typedef T Scalar;

    explicit gsFeSolution(const gsFeSpace<T> & u)
    : //Base(u.dim())
      //Base(u)
      _u(u), _Sv(NULL)
    {
        //registerData(u.source(), u.data());
    }
    
    gsFeSolution(const gsFeSpace<T> & u, gsMatrix<T> & Sv) : _u(u), _Sv(&Sv)
    { }

    mutable gsMatrix<T> res;
    const gsMatrix<T> & eval(index_t k) const
    {
        GISMO_ASSERT(1==_u.data().actives.cols(), "Single actives expected");
        
        res.setZero(_u.dim(), 1);
        const gsDofMapper & map = _u.mapper();
        for (index_t i = 0; i!=_u.data().actives.size(); ++i)
        {
            const index_t ii = map.index(_u.data().actives.at(i), _u.data().patchId);
            if ( map.is_free_index(ii) ) // DoF value is in the solVector
            {
                for (index_t r = 0; r != res.rows(); ++r)
                {
                    const index_t cgs = r * map.freeSize();
                    res.at(r) += _Sv->at(cgs+ii) * _u.data().values[0](i,k);
                }
            }
            else
                res.noalias() += _u.data().values[0](i,k) *
                    _u.fixedPart().row( map.global_to_bindex(ii) ).transpose(); //head(_u.dim());
        }
        return res;
    }

    bool rowSpan() const {return false; }
    bool colSpan() const {return false;}

    index_t rows() const {return _u.dim(); }
    index_t cols() const {return 1; }

    void setFlag() const
    {
        _u.data().flags |= NEED_VALUE | NEED_ACTIVE;
    }
        
    void parse(gsSortedVector<const gsFunctionSet<Scalar>*> & evList) const
    {
        //GISMO_ASSERT(NULL!=m_fd, "FeVariable: FuncData member not registered");
        evList.push_sorted_unique(&_u.source());
        _u.data().flags |= NEED_VALUE | NEED_ACTIVE;
    }

    void print(std::ostream &os) const { os << "s"; }
    
public:
    index_t dim() const { return _u.m_d;}
    
    index_t parDim() const
    {
        return _u.m_fs->domainDim();
        //return _u.m_fd->dim.first;
    }

    gsDofMapper & mapper() {return _u.m_mapper;}
    const gsDofMapper & mapper() const {return _u.m_mapper;}

    inline const gsMatrix<T> & fixedPart() const {return _u.m_fixedDofs;}
    gsMatrix<T> & fixedPart() {return _u.m_fixedDofs;}

    gsFuncData<T> & data() {return *_u.m_fd;}
    const gsFuncData<T> & data() const {return *_u.m_fd;}
    
    void setSolutionVector(const gsMatrix<T>& solVector)
    { _Sv = & solVector; }

    const gsMatrix<T> & solutionVector() const { return *_Sv; }

    void getCoeffs(gsMatrix<T> & result, const index_t p = 0) const
    {
        GISMO_ASSERT(NULL!=_Sv, "Solution vector not set, call setSolutionVector");
        _u.getCoeffs(*_Sv, result, p);
    }

    void extract(gsMultiPatch<T> & result) const
    {
		result.clear();

		const gsMultiBasis<T>* basis = dynamic_cast<const gsMultiBasis<T>* >(&_u.source());
		for (size_t i = 0; i < basis->nBases(); i++)
		{
			memory::unique_ptr<gsGeometry<T> > p(this->extractPiece(i));
			result.addPatch(*p);
		}
    }

    memory::unique_ptr<gsGeometry<T> > extractPiece(const index_t p) const
    {
        if ( const gsBasis<T> * b = dynamic_cast<const gsBasis<T>*>(&_u.source().piece(p)) )
        {
            gsMatrix<T> cf;
            getCoeffs(cf, p);
            return b->makeGeometry(give(cf));
        }
        GISMO_ERROR("gsFeSolution: Extraction error");
    }

    // insert g-coefficients to solution vector
    void insert(const gsGeometry<T> & g, const index_t p = 0) const
    {
        const index_t dim = _u.dim();
        const gsMatrix<T> & cf = g.coefs();
        gsMatrix<T> & sol = *_Sv;
        //gsMatrix<T> & fixedPart = _u.fixedPart();
        const gsDofMapper & mapper = _u.mapper();
        for (index_t i = 0; i != cf.rows(); ++i)
        {
            const index_t ii = mapper.index(i, p);
            if ( mapper.is_free_index(ii) ) // DoF value is in the solVector
            {
                for (index_t k = 0; k != dim; ++k)
                {
                    const index_t cgs = k * mapper.freeSize();
                    sol.at(cgs+ii) = cf(i, k);
                }
            }
            /*
              else
              {
              fixedPart.row(m_mapper.global_to_bindex(ii)) = cf.row(i);
              }
            */
        }
    }
};


template<class T>
class solGrad_expr :public _expr<solGrad_expr<T> >
{
protected:
    const gsFeSolution<T> & _u;

public:
    typedef T Scalar;

    explicit solGrad_expr(const gsFeSolution<T> & u) : _u(u) { }
    
    mutable gsMatrix<T> res;
    const gsMatrix<T> & eval(index_t k) const
    {
        GISMO_ASSERT(1==_u.data().actives.cols(), "Single actives expected");
            
        res.setZero(_u.dim(), _u.parDim());
        const gsDofMapper & map = _u.mapper();
        for (index_t i = 0; i!=_u.data().actives.size(); ++i)
        {
            const index_t ii = map.index(_u.data().actives.at(i), _u.data().patchId);
            if ( map.is_free_index(ii) ) // DoF value is in the solVector
            {
                for (index_t r = 0; r != res.rows(); ++r)
                {
                    const index_t cgs = r * map.freeSize();
                    res.row(r) += _u.solutionVector().at(cgs+ii) *
                        _u.data().values[1].block(i*_u.parDim(),k,_u.parDim(),1).transpose();
                        //col(k).segment(i*_u.parDim(), _u.parDim()).transpose();
                }
            }
            else
            {
                const index_t j = map.global_to_bindex(ii);
                for (index_t r = 0; r != res.rows(); ++r)
                {
                    res.row(r) += _u.fixedPart()(j,r) *
                        _u.data().values[1].block(i*_u.parDim(),k,_u.parDim(),1).transpose();
                        //.col(k).segment(i*_u.parDim(), _u.parDim()).transpose();
                }
                
                /*
                res.noalias() +=
                    _u.fixedPart().row( map.global_to_bindex(ii) ).asDiagonal() *
                    _u.data().values[1].col(k).segment(i*_u.parDim(), _u.parDim())
                    .transpose().replicate(_u.dim(),1);
                */
            }
        }
        return res;
    }

    bool rowSpan() const {return false;}
    bool colSpan() const {return false;}

    index_t rows() const {return _u.dim();}

    index_t cols() const {return _u.parDim(); }

    void setFlag() const
    {
        _u.data().flags |= NEED_ACTIVE;
        _u.data().flags |= NEED_GRAD | NEED_VALUE;
    }

    void parse(gsSortedVector<const gsFunctionSet<Scalar>*> & evList) const
    {
        //GISMO_ASSERT(NULL!=m_fd, "FeVariable: FuncData member not registered");
        evList.push_sorted_unique(&_u.source());
        _u.data().flags |= NEED_GRAD;
    }

    void print(std::ostream &os) const { os << "grad(s)"; }
};

/// The gradient of a solution variable
template<class T> solGrad_expr<T> grad(const gsFeSolution<T> & u)
{ return solGrad_expr<T>(u); }


//template<class T>
//class gsFeMutVariable  : public _expr<gsFeMutVariable<T> >
// member gsFuncData

template <typename T> 
struct expr_traits<gsFeVariable<T> >
{
    typedef T Scalar;
    typedef gsFeVariable<T> const & Nested_t;
};

/**
   Expression for the transpose of an expression
 */
template<class E>
//class tr_expr<E,false>  : public _expr<tr_expr<E,false> >
class tr_expr : public _expr<tr_expr<E> >
{
    typename E::Nested_t _u;
    
public:
    
    typedef typename E::Scalar Scalar;
    
    tr_expr(_expr<E> const& u)
    : _u(u) { }
    
public:
    enum {ColBlocks = E::ColBlocks};

    // template<bool S  = ColBlocks>
    // typename util::enable_if<S,MatExprType>::type
    MatExprType eval(const index_t k) const
    {
        if (E::ColBlocks)
            return _u.eval(k).blockTranspose(_u.cols()/_u.rows());
        else
            //return _u.eval(k).transpose(); // auto
            return _u.eval(k).blockTranspose(1);
    }

    index_t rows() const
    {
        if (ColBlocks) // note: sq. blocks assumed
            return _u.rows();
        else
            return _u.cols();
        
    }
    index_t cols() const
    {
        if (ColBlocks) // note: sq. blocks assumed
            return _u.cols();
        else
            return _u.rows();
    }
    void setFlag() const { _u.setFlag(); }

    void parse(gsSortedVector<const gsFunctionSet<Scalar>*> & evList) const
    { _u.parse(evList); }

    const gsFeVariable<Scalar> & rowVar() const { return _u.colVar(); }
    const gsFeVariable<Scalar> & colVar() const { return _u.rowVar(); }

    bool rowSpan() const {return _u.colSpan();}
    bool colSpan() const {return _u.rowSpan();}

    void print(std::ostream &os) const { os<<"("; _u.print(os); os <<")'"; }
};

template<class E>
class temp_expr : public _expr<temp_expr<E> >
{
    typename E::Nested_t _u;
    typedef typename E::Scalar Scalar;
    mutable gsMatrix<Scalar> tmp;

public:
    temp_expr(_expr<E> const& u)
    : _u(u) { }
    
public:
    enum {ColBlocks = E::ColBlocks};

    // template<bool S  = ColBlocks>
    // typename util::enable_if<S,MatExprType>::type
    const gsMatrix<Scalar> & eval(const index_t k) const
    {
        //GISMO_ERROR("NO!");
        //If not printing then it inlines, otheriwise NOT
        //gsInfo<<"* * * * * * * * * * * * * * * * Making tmp\n"; 
        tmp = _u.eval(k);
        return tmp;
    }

    index_t rows() const
    {
            return _u.rows();
    }
    index_t cols() const
    {
            return _u.cols();
    }

    void setFlag() const { _u.setFlag(); }

    void parse(gsSortedVector<const gsFunctionSet<Scalar>*> & evList) const
    { _u.parse(evList); }

    const gsFeVariable<Scalar> & rowVar() const { return _u.rowVar(); }
    const gsFeVariable<Scalar> & colVar() const { return _u.colVar(); }

    bool rowSpan() const {return _u.rowSpan();}
    bool colSpan() const {return _u.colSpan();}

    void print(std::ostream &os) const { _u.print(os); }
};

/**
   Expression for the trace of a (matrix) expression
 */
template<class E>
class trace_expr  : public _expr<trace_expr<E> >
{
public:
    typedef typename E::Scalar Scalar;
    enum {ScalarValued = 0};
private:
    typename E::Nested_t _u;
    mutable gsMatrix<Scalar> res;
    
public:
    trace_expr(_expr<E> const& u) : _u(u)
    { GISMO_ASSERT(0== _u.cols()%_u.rows(), "Expecting square-block expression, got "
                   << _u.rows() <<" x "<< _u.cols() ); }
    
    // choose if ColBlocks
    const gsMatrix<Scalar> & eval(const index_t k) const
    {
        // Assume mat ??
        MatExprType tmp = _u.eval(k);
        const index_t cb = _u.rows();
        const index_t r  = _u.cols() / cb;
        res.resize(r, 1);
        for (index_t i = 0; i!=r; ++i)
            res(i,0) = tmp.middleCols(i*cb,cb).trace();
        return res;
    }

    // choose if !ColBlocks
    //todo: Scalar eval(const index_t k) const
    
    index_t rows() const { return _u.cols() / _u.rows(); }
    index_t cols() const { return 1; }
    void setFlag() const { _u.setFlag(); }

    void parse(gsSortedVector<const gsFunctionSet<Scalar>*> & evList) const
    { _u.parse(evList); }

    const gsFeVariable<Scalar> & rowVar() const { return _u.rowVar(); }
    const gsFeVariable<Scalar> & colVar() const { return _u.colVar(); }

    bool rowSpan() const {return _u.rowSpan();}
    bool colSpan() const {return _u.colSpan();}

    void print(std::ostream &os) const { os << "trace("; _u.print(os); os<<")"; }
};

/**
   Expression for the diagonal(s) of a (matrix) expression

template<class E>
class diag_expr  : public _expr<diag_expr<E> >
{
public:
    typedef typename E::Scalar Scalar;
    enum {ScalarValued = 0};
private:
    typename E::Nested_t _u;
    mutable gsMatrix<Scalar> res;
    
public:
    diag_expr(_expr<E> const& u) : _u(u)
    { GISMO_ASSERT(0== _u.cols()%_u.rows(), "Expecting square-block expression, got "
                   << _u.rows() <<" x "<< _u.cols() ); }
    
    // choose if ColBlocks
    const gsMatrix<Scalar> & eval(const index_t k) const
    {
        // Assume mat ??
        MatExprType tmp = _u.eval(k);
        const index_t cb = _u.rows();
        const index_t r  = _u.cols() / cb;
        res.resize(r, cb);
        for (index_t i = 0; i!=r; ++i)
            res.row(i) = tmp.middleCols(i*cb,cb).diagonal();
        return res;
    }

    // choose if !ColBlocks
    //todo: Scalar eval(const index_t k) const
    
    index_t rows() const { return _u.cols() / _u.rows(); }
    index_t cols() const { return 1; }
    void setFlag() const { _u.setFlag(); }

    void parse(gsSortedVector<const gsFunctionSet<Scalar>*> & evList) const
    { _u.parse(evList); }

    const gsFeVariable<Scalar> & rowVar() const { return _u.rowVar(); }
    const gsFeVariable<Scalar> & colVar() const { return _u.colVar(); }

    bool rowSpan() const {return _u.rowSpan();}
    bool colSpan() const {return _u.colSpan();}

    void print(std::ostream &os) const { os << "trace("; _u.print(os); os<<")"; }
};
 */

#define GISMO_EXPR_VECTOR_EXPRESSION(name, mname, isSv)                 \
template<class E> class name##_##expr  : public _expr<name##_##expr<E> > { \
    typename E::Nested_t _u;                                            \
public:                                                                 \
    typedef typename E::Scalar Scalar;                                  \
    enum {ScalarValued = isSv};                                         \
    name##_##expr(_expr<E> const& u) : _u(u) { }                        \
    AutoReturn_t eval(const index_t k) const { return _u.eval(k).mname();} \
    index_t rows() const { return isSv ? 0 : _u.rows(); }               \
    index_t cols() const { return isSv ? 0 : _u.cols(); }               \
    void setFlag() const { _u.setFlag(); }                              \
    void parse(gsSortedVector<const gsFunctionSet<Scalar>*> & evList) const { _u.parse(evList); } \
    const gsFeVariable<Scalar> & rowVar() const {GISMO_ERROR("Invalid call to rowVar()");} \
    const gsFeVariable<Scalar> & colVar() const {GISMO_ERROR("Invalid call to colVar()");} \
    void print(std::ostream &os) const                                  \
    { os << #name <<"("; _u.print(os); os <<")"; }                      \
    bool rowSpan() const {return _u.rowSpan();}                         \
    bool colSpan() const {return _u.colSpan();} };
    // const gsFeVariable<Scalar> & rowVar() const { return _u.rowVar(); }
    // const gsFeVariable<Scalar> & colVar() const { return _u.colVar(); }

GISMO_EXPR_VECTOR_EXPRESSION(norm,norm,1)
GISMO_EXPR_VECTOR_EXPRESSION(sqNorm,squaredNorm,1)
GISMO_EXPR_VECTOR_EXPRESSION(normalized,normalized,0) // (!) mem.
GISMO_EXPR_VECTOR_EXPRESSION(inv,inverse,0)
// GISMO_EXPR_VECTOR_EXPRESSION(cwSqr,array().square,0)
// GISMO_EXPR_VECTOR_EXPRESSION(sum,array().sum,1)
// GISMO_EXPR_VECTOR_EXPRESSION(sqrt,array().sqrt,0)
// GISMO_EXPR_VECTOR_EXPRESSION(abs,array().abs,0)

#undef GISMO_EXPR_VECTOR_EXPRESSION

/**
   Expression for turning a vector into a diagonal matrix
 */
template<class E>
class asdiag_expr : public _expr<asdiag_expr<E> >
{
public:
    typedef typename E::Scalar Scalar;
private:
    typename E::Nested_t _u;
    mutable gsMatrix<Scalar> res;

public:
    asdiag_expr(_expr<E> const& u) : _u(u) { }
    
public:

    MatExprType eval(const index_t k) const
    {
        MatExprType m = _u.eval(k);
        const index_t r = m.rows();
        const index_t c = m.cols();
        res.resize(r,r*c);
        for (index_t i = 0; i!=c; ++i)
            res.middleCols(i*r,r) = m.col(i).asDiagonal();
        return res;
    }

    const gsFeVariable<Scalar> & rowVar() const { return _u.rowVar(); }
    const gsFeVariable<Scalar> & colVar() const { return _u.colVar(); }

    bool rowSpan() const {return _u.rowSpan();}
    bool colSpan() const {return _u.colSpan();}
    
    index_t rows() const { return _u.rows(); }
    index_t cols() const { return _u.rows() * _u.cols(); }
    void setFlag() const { _u.setFlag(); }
    void parse(gsSortedVector<const gsFunctionSet<Scalar>*> & evList) const
    { _u.parse(evList); }

    void print(std::ostream &os) const { os << "diag("; _u.print(os); os <<")";}
};

/**
   Expression for turning a vector into a diagonal matrix
 */
class idMat_expr : public _expr<idMat_expr >
{
public:
    typedef real_t Scalar;
private:
    index_t _dim;

public:
    idMat_expr(const index_t dim) : _dim(dim) { }
    
public:

    MatExprType eval(const index_t k) const
    {
        return gsMatrix<Scalar>::Identity(_dim,_dim);
    }

    index_t rows() const { return _dim; }
    index_t cols() const { return  _dim; }
    void setFlag() const { }
    void parse(gsSortedVector<const gsFunctionSet<Scalar>*> & evList) const {  }

    bool rowSpan() const {return false;}
    bool colSpan() const {return false;}

    const gsFeVariable<Scalar> & rowVar() const
    {GISMO_ERROR("Invalid call to value::rowVar()");}
    const gsFeVariable<Scalar> & colVar() const
    {GISMO_ERROR("Invalid call to value::colVar()");}

    void print(std::ostream &os) const { os << "id("<<_dim <<")";}
};

/// The identity matrix of dimension \a dim
EIGEN_STRONG_INLINE idMat_expr id(const index_t dim) { return idMat_expr(dim); }

/**
   Adaptor for scalar-valued expression
 */
template<class E>
class value_expr  : public _expr<value_expr<E> >
{
    typename E::Nested_t _u;

public:
    typedef typename E::Scalar Scalar;    
    value_expr(_expr<E> const& u) : _u(u)
    {
        // rows/cols not known at construction time
        //GISMO_ASSERT(u.rows()*u.cols()<=1, "Expression\n"<<u<<"is not a scalar.");
    }
    
public:
    enum {ScalarValued = 1};

    Scalar eval(const index_t k) const { return eval_impl(_u,k); }
    
    // enables expr.val().val()
    inline value_expr<E> val() const { return *this; }
    index_t rows() const { return 0; }
    index_t cols() const { return 0; }
    void setFlag() const { _u.setFlag(); }
    void parse(gsSortedVector<const gsFunctionSet<Scalar>*> & evList) const
    { _u.parse(evList); }

    static bool isScalar() { return true; }
    
    bool rowSpan() const {return false;}
    bool colSpan() const {return false;}

    const gsFeVariable<Scalar> & rowVar() const
    {GISMO_ERROR("Invalid call to value::rowVar()");}
    const gsFeVariable<Scalar> & colVar() const
    {GISMO_ERROR("Invalid call to value::colVar()");}

    void print(std::ostream &os) const { _u.print(os); }

    // Math functions. eg
    // sqrt_mexpr<T> sqrt() { return sqrt_mexpr<T>(*this); }
private:
    template<class U> static inline
    typename util::enable_if<U::ScalarValued,Scalar>::type
    eval_impl(const U & u, const index_t k) { return u.eval(k); }
    
    template<class U> static inline
    typename util::enable_if<!U::ScalarValued,Scalar>::type
    eval_impl(const U & u, const index_t k) { return u.eval(k).value(); }
};

/**
   Expression for the gradient of a finite element variable

   Transposed gradient vectors are returned as a matrix
 */
template<class T>
class grad_expr : public _expr<grad_expr<T> >
{
    typename gsFeVariable<T>::Nested_t _u;
    
public:    
    typedef T Scalar;
    
    grad_expr(const gsFeVariable<T> & u) : _u(u)
    { GISMO_ASSERT(1==u.dim(),"grad(.) requires 1D variable, use jac(.) instead.");}

    MatExprType eval(const index_t k) const
    {
        // numActive x dim
        return _u.data().values[1].reshapeCol(k, cols(), rows()).transpose();
    }

    index_t rows() const { return _u.data().values[0].rows(); }
    //index_t rows() const { return _u.data().actives.size(); }
    //index_t rows() const { return _u.rows(); }
    //index_t rows() const { return _u.data().dim.second; }

    index_t cols() const { return _u.data().dim.first; }
    void setFlag() const { _u.data().flags |= NEED_GRAD|NEED_ACTIVE; }

    void parse(gsSortedVector<const gsFunctionSet<Scalar>*> & evList) const
    {
        //GISMO_ASSERT(NULL!=m_fd, "FeVariable: FuncData member not registered");
        evList.push_sorted_unique(&_u.source());
        _u.data().flags |= NEED_GRAD|NEED_ACTIVE;
    }

    const gsFeVariable<T> & rowVar() const { return _u.rowVar(); }
    const gsFeVariable<T> & colVar() const
    {GISMO_ERROR("Invalid call to gsFeVariable::colVar()");}

    bool rowSpan() const {return true; }
    bool colSpan() const {return false;}

    void print(std::ostream &os) const { os << "grad("; _u.print(os); os <<")"; }
};

/// The gradient of a finite element variable
template<class T>
grad_expr<T> grad(const gsFeVariable<T> & u)
{ return grad_expr<T>(u); }

/**
   Expression for the nabla of a finite element variable,
 */
template<class T>
class nabla_expr : public _expr<nabla_expr<T> >
{
    typename gsFeVariable<T>::Nested_t u;

public:
    
    typedef T Scalar;

    /* // todo
    nabla_expr(const gsGeometryMap<T> & G)
    : m_data(G.data()) { }
    */
    
    nabla_expr(const gsFeVariable<T> & _u) : u(_u)
    {
        //GISMO_ASSERT(u.parDim()==u.dim(),"nabla(.) requires tarDim==parDim:"
        //             << u.parDim() <<"!="<< u.dim() <<"\n" );
    }

    mutable gsMatrix<Scalar> res;
    
    const gsMatrix<Scalar> & eval(const index_t k) const
    {
        const index_t d = u.cols();
        const index_t n = rows() / d;
        res.setZero(rows(), d);
        
        for (index_t i = 0; i!=d; ++i)
            res.col(i).segment(i*n,n) = u.data().values[1].reshapeCol(k, d, n).row(i);
        return res;
    }

    index_t rows() const { return u.rows(); }
    index_t cols() const { return u.cols(); }
    void setFlag() const { u.data().flags |= NEED_GRAD; }

    void parse(gsSortedVector<const gsFunctionSet<Scalar>*> & evList) const
    {
        //GISMO_ASSERT(NULL!=m_fd, "FeVariable: FuncData member not registered");
        evList.push_sorted_unique(&u.source());
        u.data().flags |= NEED_GRAD;
    }

    const gsFeVariable<T> & rowVar() const { return u.rowVar(); }
    const gsFeVariable<T> & colVar() const
    {GISMO_ERROR("Invalid call to nabla::colVar()");}

    bool rowSpan() const {return true; }
    bool colSpan() const {return false;}

    void print(std::ostream &os) const { os << "nabla("; u.print(os); os <<")"; }
};

/// The nabla (nabla) of a finite element variable
template<class T>
nabla_expr<T> nabla(const gsFeVariable<T> & u)
{ return nabla_expr<T>(u); }

/**
   Expression for the nabla^2 (or Del^2) of a finite element variable,
   see also https://en.wikipedia.org/wiki/Del

   Transposed pure second derivatives are returned as a matrix
 */
template<class T>
class nabla2_expr : public _expr<nabla2_expr<T> >
{
    typename gsFeVariable<T>::Nested_t u;

public:
    
    typedef T Scalar;

    /* // todo
    nabla2_expr(const gsGeometryMap<T> & G)
    : m_data(G.data()) { }
    */
    
    nabla2_expr(const gsFeVariable<T> & _u)
    : u(_u)
    { }

    MatExprType eval(const index_t k) const
    {
        // numActive x parDim
        return u.data().values[2]
            .reShapeCol(k, u.data().values[2].rows()/u.cSize(), u.cSize() )
            .topRows(u.parDim()).transpose();
    }

    index_t rows() const { return u.rows();   }
    index_t cols() const { return u.parDim(); }
    void setFlag() const { u.data().flags |= NEED_DERIV2; }
    void parse(gsSortedVector<const gsFunctionSet<Scalar>*> & evList) const
    {
        //GISMO_ASSERT(NULL!=m_fd, "FeVariable: FuncData member not registered");
        evList.push_sorted_unique(&u.source());
        u.data().flags |= NEED_DERIV2;
    }

    const gsFeVariable<T> & rowVar() const { return u.rowVar(); }
    const gsFeVariable<T> & colVar() const
    {GISMO_ERROR("Invalid call to nabla2::colVar()");}

    bool rowSpan() const {return true; }
    bool colSpan() const {return false;}
};

/// The nabla2 (nabla^2) of a finite element variable
template<class T>
nabla2_expr<T> nabla2(const gsFeVariable<T> & u)
{ return nabla2_expr<T>(u); }
// #define lapl(x) nabla2(x).sum() // assume tarDim==1

/**
   Expression for the outer pointing normal of of a geometry map. This
   expression is valid only at the boundaries of a geometric patch
 */
template<class T>
class onormal_expr : public _expr<onormal_expr<T> >
{
    typename gsGeometryMap<T>::Nested_t _G;

public:
    typedef T Scalar;
    
    onormal_expr(const gsGeometryMap<T> & G) : _G(G) { }

    MatExprType eval(const index_t k) const
    {
        return _G.data().outNormals.col(k);
    }
    
    index_t rows() const { return _G.data().dim.second; }
    index_t cols() const { return 1; }

    bool rowSpan() const {GISMO_ERROR("onormal");}
    bool colSpan() const {GISMO_ERROR("onormal");}
    
    void setFlag() const { _G.data().flags |= NEED_OUTER_NORMAL; }

    void parse(gsSortedVector<const gsFunctionSet<Scalar>*> & evList) const
    {
        //GISMO_ASSERT(NULL!=m_fd, "FeVariable: FuncData member not registered");
        evList.push_sorted_unique(&_G.source());
        _G.data().flags |= NEED_OUTER_NORMAL;
    }

    // Normalized to unit length
    normalized_expr<onormal_expr<T> > normalized()
    { return normalized_expr<onormal_expr<T> >(*this); }

    void print(std::ostream &os) const { os << "nv("; _G.print(os); os <<")"; }
};

/// The (outer pointing) boundary normal of a geometry map
// nv(.), nor(.)
template<class T> EIGEN_STRONG_INLINE
onormal_expr<T> nv(const gsGeometryMap<T> & u) { return onormal_expr<T>(u); }

/**
   Expression for the tangent vector of a geometry map. This
   expression is valid only at the boundaries of a geometric patch
 */
template<class T>
class tangent_expr : public _expr<tangent_expr<T> >
{
    typename gsGeometryMap<T>::Nested_t _G;

public:
    typedef T Scalar;
    
    tangent_expr(const gsGeometryMap<T> & G) : _G(G) { }

    mutable gsMatrix<Scalar> res;

    MatExprType eval(const index_t k) const
    {
        res = _G.data().outNormals.col(k);//2x1
        std::swap( res(0,0), res(1,0) );
        res(0,0) *= -1;
        return res;
    }
    
    index_t rows() const { return _G.data().dim.second; }
    index_t cols() const { return 1; }

    bool rowSpan() const {GISMO_ERROR("tangent");}
    bool colSpan() const {GISMO_ERROR("tangent");}
    
    void setFlag() const { _G.data().flags |= NEED_OUTER_NORMAL; }

    void parse(gsSortedVector<const gsFunctionSet<Scalar>*> & evList) const
    {
        //GISMO_ASSERT(NULL!=m_fd, "FeVariable: FuncData member not registered");
        evList.push_sorted_unique(&_G.source());
        _G.data().flags |= NEED_OUTER_NORMAL;
    }

    // Normalized to unit length
    normalized_expr<tangent_expr<T> > normalized()
    { return normalized_expr<tangent_expr<T> >(*this); }

    void print(std::ostream &os) const { os << "tv("; _G.print(os); os <<")"; }
};

/// The tangent boundary normal of a geometry map
template<class T> EIGEN_STRONG_INLINE
tangent_expr<T> tv(const gsGeometryMap<T> & u) { return tangent_expr<T>(u); }

/**
   Expression for the Laplacian of a finite element variable
 */
template<class T>
class lapl_expr : public _expr<lapl_expr<T> >
{
    typename gsFeVariable<T>::Nested_t _u;

public:
    typedef T Scalar;
    
    lapl_expr(const gsFeVariable<T> & u) : _u(u) { }

    MatExprType eval(const index_t k) const
    {
        // numActive x 1
        return _u.data().laplacians.col(k);
        //todo: replace by
        // NEED_DERIV2
        // ..nabla2.sum()
    }

    index_t rows() const { return _u.data().values[0].rows(); }
    index_t cols() const { return 1; }

    bool rowSpan() const {return true; }
    bool colSpan() const {return false;}
    
    void setFlag() const { _u.data().flags |= NEED_LAPLACIAN; }

    void parse(gsSortedVector<const gsFunctionSet<Scalar>*> & evList) const
    {
        //GISMO_ASSERT(NULL!=m_fd, "FeVariable: FuncData member not registered");
        evList.push_sorted_unique(&_u.source());
        _u.data().flags |= NEED_LAPLACIAN;
    }

    void print(std::ostream &os) const { os << "lap("; _u.print(os); os <<")"; }
};

/// The laplacian of a finite element variable
template<class T> EIGEN_STRONG_INLINE
lapl_expr<T> lapl(const gsFeVariable<T> & u)
{ return lapl_expr<T>(u); }

/*
   Expression for the (precomputed) first fundamental form of a geometry map
*/
template<class T>
class fform_expr  : public _expr<fform_expr<T> >
{
    typename gsGeometryMap<T>::Nested_t _G;
public:
    typedef T Scalar;

    fform_expr(const gsGeometryMap<T> & G) : _G(G) { }
    
    MatExprType eval(const index_t k) const // todo: fix funcdata
    { return _G.data().fundForm(k).transpose() * _G.data().fundForm(k) ; }
    
    index_t rows() const { return _G.data().dim.second; }
    index_t cols() const { return _G.data().dim.first ; }

    bool rowSpan() const {return false; }
    bool colSpan() const {return false;}

    void setFlag() const { _G.data().flags |= NEED_GRAD_TRANSFORM; }
    
    void parse(gsSortedVector<const gsFunctionSet<Scalar>*> & evList) const
    {
        //GISMO_ASSERT(NULL!=m_fd, "FeVariable: FuncData member not registered");
        evList.push_sorted_unique(&_G.source());
        _G.data().flags |= NEED_GRAD_TRANSFORM;
    }

    const gsFeVariable<Scalar> & rowVar() const {GISMO_ERROR("Invalid call to rowVar()");}
    const gsFeVariable<Scalar> & colVar() const {GISMO_ERROR("Invalid call to colVar()");}

    void print(std::ostream &os) const { os << "fform("; _G.print(os); os <<")"; }
};

/// The first fundamental form of \a G
template<class T> EIGEN_STRONG_INLINE // or G.fform() ?
fform_expr<T> fform(const gsGeometryMap<T> & G)
{ return fform_expr<T>(G); }

/*
   Expression for the (precomputed) inverse of the Jacobian matrix of
   a geometry map
*/
template<class T>
class jacGinv_expr  : public _expr<jacGinv_expr<T> >
{
    typename gsGeometryMap<T>::Nested_t _G;
public:
    typedef T Scalar;

    jacGinv_expr(const gsGeometryMap<T> & G) : _G(G)
    {
        // Note: for non-square Jacobian matrices, generalized inverse, i.e.: (J^t J)^{-t} J^t
        //GISMO_ASSERT(rows() == cols(), "The Jacobian matrix is not square");
    }
    
    MatExprType eval(const index_t k) const { return _G.data().fundForm(k).transpose(); }
    
    index_t rows() const { return _G.data().dim.first;  }
    index_t cols() const { return _G.data().dim.second; }

    bool rowSpan() const {return false; }
    bool colSpan() const {return false;}

    void setFlag() const { _G.data().flags |= NEED_GRAD_TRANSFORM; }
    
    void parse(gsSortedVector<const gsFunctionSet<Scalar>*> & evList) const
    {
        //GISMO_ASSERT(NULL!=m_fd, "FeVariable: FuncData member not registered");
        evList.push_sorted_unique(&_G.source());
        _G.data().flags |= NEED_GRAD_TRANSFORM;
    }

    const gsFeVariable<Scalar> & rowVar() const {GISMO_ERROR("Invalid call to rowVar()");}
    const gsFeVariable<Scalar> & colVar() const {GISMO_ERROR("Invalid call to colVar()");}

    // todo mat_expr ?
    // tr() const --> _G.data().fundForm(k)
    
    void print(std::ostream &os) const { os << "jacInv("; _G.print(os); os <<")"; }
};

// #define GISMO_FUNCDATA_EXPRESSION( ...

/**
   Expression for the Jacobian matrix of a geometry map
 */
template<class T>
class jacG_expr : public _expr<jacG_expr<T> >
{
    typename gsGeometryMap<T>::Nested_t _G;
    
public:
    typedef T Scalar;
    
    jacG_expr(const gsGeometryMap<T> & G) : _G(G) { }

    MatExprType eval(const index_t k) const
    {
        // TarDim x ParDim
        return _G.data().values[1]
            .reshapeCol(k, _G.data().dim.first, _G.data().dim.second).transpose();
    }

    index_t rows() const { return _G.data().dim.second; }
    index_t cols() const { return _G.data().dim.first; }

    bool rowSpan() const {return false; }
    bool colSpan() const {return false;}
    
    void setFlag() const { _G.data().flags |= NEED_DERIV; }

    void parse(gsSortedVector<const gsFunctionSet<Scalar>*> & evList) const
    {
        //GISMO_ASSERT(NULL!=m_fd, "FeVariable: FuncData member not registered");
        evList.push_sorted_unique(&_G.source());
        _G.data().flags |= NEED_DERIV;
    }

    meas_expr<T> det() const
    {
        GISMO_ASSERT(rows() == cols(), "The Jacobian matrix is not square");
        return meas_expr<T>(_G);
    }

    jacGinv_expr<T> inv() const
    {
        GISMO_ASSERT(rows() == cols(), "The Jacobian matrix is not square");
        return jacGinv_expr<T>(_G);
    }

    /// The generalized Jacobian matrix inverse, i.e.: (J^t J)^{-t} J^t
    jacGinv_expr<T> ginv() const { return jacGinv_expr<T>(_G); }

    void print(std::ostream &os) const { os << "jac_("; _G.print(os); os <<")"; }
};

/// The Jacobian matrix of a geometry map
template<class T> EIGEN_STRONG_INLINE
jacG_expr<T> jac(const gsGeometryMap<T> & G)
{ return jacG_expr<T>(G); }

/**
   Expression for the Jacobian matrix of a FE variable
 */
template<class T>
class jac_expr : public _expr<jac_expr<T> >
{
    const gsFeVariable<T> & m_fev;
public:
    enum {ColBlocks = 1};
    
    typedef T Scalar;
    
    jac_expr(const gsFeVariable<T> & _u)
    : m_fev(_u) { }

    MatExprType eval(const index_t k) const
    {
        // (Dim x ParDim) * (Dim*numActive)
        return m_fev.data().values[1].col(k).transpose().blockDiag(m_fev.dim());
    }
    
    const gsFeVariable<T> & rowVar() const { return m_fev; }
    const gsFeVariable<T> & colVar() const { GISMO_ERROR("jac is a row expression?"); }

    index_t rows() const { return m_fev.dim(); }
    index_t cols() const
    {   //bug
        return m_fev.dim() * m_fev.data().actives.rows() * m_fev.data().dim.first;
    }

    bool rowSpan() const {return true; }
    bool colSpan() const {return false;}
    
    void setFlag() const
    {
        m_fev.data().flags |= NEED_DERIV;
        m_fev.data().flags |= NEED_ACTIVE;// rows() depend on this
    }
    
    void parse(gsSortedVector<const gsFunctionSet<Scalar>*> & evList) const
    {
        //GISMO_ASSERT(NULL!=m_fd, "FeVariable: FuncData member not registered");
        evList.push_sorted_unique(& m_fev.source());
         m_fev.data().flags |= NEED_DERIV;
         m_fev.data().flags |= NEED_ACTIVE;// rows() depend on this
    }

    //jacInv_expr<T> inv()
    //meas_expr<T> det() const

    void print(std::ostream &os) const { os << "jac("; m_fev.print(os);os <<")"; }
};

/// The Jacobian matrix of a FE variable
template<class T> EIGEN_STRONG_INLINE
jac_expr<T> jac(const gsFeVariable<T> & u)
{ return jac_expr<T>(u); }


template<class T>
class fjac_expr : public _expr<fjac_expr<T> >
{
    const gsFeVariable<T> & m_fev;
public:
    enum {ColBlocks = 0};// main difference from jac
    
    typedef T Scalar;
    
    fjac_expr(const gsFeVariable<T> & _u)
    : m_fev(_u)
    { }

    MatExprType eval(const index_t k) const
    {
        return m_fev.data().values[1].reshapeCol(k, cols(), rows()).transpose();
    }
    
    const gsFeVariable<T> & rowVar() const { return m_fev; }
    const gsFeVariable<T> & colVar() const { GISMO_ERROR("jac is a row expression?"); }

    index_t rows() const { return m_fev.data().dim.first; }
    index_t cols() const {return m_fev.data().dim.second; }

    bool rowSpan() const {return true; }
    bool colSpan() const {return false;}
    
    void setFlag() const
    {
        m_fev.data().flags |= NEED_DERIV;
    }
    
    void parse(gsSortedVector<const gsFunctionSet<Scalar>*> & evList) const
    {
        //GISMO_ASSERT(NULL!=m_fd, "FeVariable: FuncData member not registered");
        evList.push_sorted_unique(& m_fev.source());
        m_fev.data().flags |= NEED_DERIV;
    }

    //jacInv_expr<T> inv()
    //meas_expr<T> det() const

    void print(std::ostream &os) const { os << "fjac("; m_fev.print(os);os <<")"; }
};

/// The Jacobian matrix of a FE variable
template<class T> EIGEN_STRONG_INLINE
fjac_expr<T> fjac(const gsFeVariable<T> & u)
{ return fjac_expr<T>(u); }

/**
   Expression for the Hessian matrix of (a coordinate of) the geometry
   map or a finite element variable
 */
template<class T>
class hess_expr : public _expr<hess_expr<T> >
{
public:
    typedef T Scalar;
    enum {ScalarValued = 0, ColBlocks = 1};

private:
    const gsFuncData<T> & m_data;
    mutable gsMatrix<Scalar> res;

public:
    hess_expr(const gsGeometryMap<T> & G)
    : m_data(G.data()) { } //ColBlocks=0 ?

    hess_expr(const gsFeVariable<T> & _u)
    : m_data(_u.data())
    { GISMO_ASSERT(1==_u.dim(),"hess(.) requires 1D variable");}

    const gsMatrix<Scalar> & eval(const index_t k) const
    {
        const index_t sz = m_data.values[0].rows();
        res.resize(m_data.dim.first, sz*m_data.dim.first);
        secDerToHessian(m_data.values[2].col(k), m_data.dim.first, res);
        res.resize(m_data.dim.first, res.cols()*m_data.dim.first);
        // Note: auto returns by value here
        return res;
    }

    index_t rows() const { return m_data.dim.first; }
    index_t cols() const { return m_data.values[0].rows() * m_data.dim.first; }
    void setFlag() const { m_data.flags |= NEED_2ND_DER; }

    void parse(gsSortedVector<const gsFunctionSet<Scalar>*> & evList) const
    {
        GISMO_ERROR("error 1712");
    }

    bool rowSpan() const {return true; }
    bool colSpan() const {return false;}

    const gsFeVariable<T> & rowVar() const { GISMO_ERROR("hess problem rowVar"); }
    const gsFeVariable<T> & colVar() const { GISMO_ERROR("hess problem colVar"); }

    void print(std::ostream &os) const
    //    { os << "hess("; _u.print(os);os <<")"; }
    { os << "hess(U)"; }
};

/// The Hessian matrix of (all coordianates of) the geometry map \a G
template<class T> EIGEN_STRONG_INLINE
hess_expr<T> hess(const gsGeometryMap<T> & G)
{ return hess_expr<T>(G); }

/// The Hessian matrix of a finite element variable
template<class T> EIGEN_STRONG_INLINE
hess_expr<T> hess(const gsFeVariable<T> & u)
{ return hess_expr<T>(u); }


/**
   Expression for the partial derivative (matrices) of the Jacobian
   matrix of the geometry map
 */
template<class T>
class dJacG_expr : public _expr<dJacG_expr<T> >
{
    typename gsGeometryMap<T>::Nested_t _G;

    mutable gsMatrix<T> res;
public:
    typedef T Scalar;
    
    dJacG_expr(const gsGeometryMap<T> & G) : _G(G) { }

    MatExprType eval(const index_t k) const
    {
        const index_t sz = _G.data().values[0].rows();
        const index_t s = _G.data().derivSize(); //dim.first*(_G.data().dim.first+1)/2;
        (void)s;
        res.resize(_G.data().dim.second, sz*_G.data().dim.first);
        res.setOnes();//_G.data().values[2].segment(i*k,k); // todo
        return res;
    }

    index_t rows() const { return _G.data().dim.second; }
    index_t cols() const { return _G.data().dim.first; }
    void setFlag() const { _G.data().flags |= NEED_2ND_DER; }

    void parse(gsSortedVector<const gsFunctionSet<Scalar>*> & evList) const
    {
        //GISMO_ASSERT(NULL!=m_fd, "FeVariable: FuncData member not registered");
        evList.push_sorted_unique(&_G.source());
        _G.data().flags |= NEED_2ND_DER;
    }

    bool rowSpan() const {return false;}
    bool colSpan() const {return false;}
};


/// The partial derivatives of the Jacobian matrix of a geometry map
template<class T> EIGEN_STRONG_INLINE
dJacG_expr<T> dJac(const gsGeometryMap<T> & G)
{ return dJacG_expr<T>(G); }

/**
   Expression for the measure of a geometry map
 */
template<class T>
class meas_expr : public _expr<meas_expr<T> >
{
    typename gsGeometryMap<T>::Nested_t _G;

public:
    enum {ScalarValued = 1};
    
    typedef T Scalar;
    
    meas_expr(const gsGeometryMap<T> & G) : _G(G) { }

    T eval(const index_t k) const { return _G.data().measures.at(k); }

    index_t rows() const { return 0; }    
    index_t cols() const { return 0; }
    void setFlag() const { _G.data().flags |= NEED_MEASURE; }
    void parse(gsSortedVector<const gsFunctionSet<Scalar>*> & evList) const
    {
        //GISMO_ASSERT(NULL!=m_fd, "FeVariable: FuncData member not registered");
        evList.push_sorted_unique(&_G.source());
        _G.data().flags |= NEED_MEASURE;
    }

    const gsFeVariable<T> & rowVar() const { GISMO_ERROR("meas problem rowVar"); }
    const gsFeVariable<T> & colVar() const { GISMO_ERROR("meas problem colVar"); }

    bool rowSpan() const {return false; }
    bool colSpan() const {return false;}

    void print(std::ostream &os) const { os << "meas("; _G.print(os); os <<")"; }
};

/// The measure of a geometry map
template<class T> EIGEN_STRONG_INLINE
meas_expr<T> meas(const gsGeometryMap<T> & G)
{ return meas_expr<T>(G); }

/**
   Expression for multiplication operation

   mult_v1

   bool ScalarValued --> if not then TMP
 */
template <typename E1, typename E2>
class mult_expr<E1,E2,false> : public _expr<mult_expr<E1, E2, false> >
{
    typename E1::Nested_t _u;
    typename E2::Nested_t _v;

public:
    enum {ScalarValued = E1::ScalarValued && E2::ScalarValued,
          ColBlocks = E2::ColBlocks};
    
    typedef typename E1::Scalar Scalar;

    typedef typename
    util::conditional<ScalarValued,Scalar,typename gsMatrix<Scalar>::Base >::type Temporary_t;
    
    mutable Temporary_t tmp;    

    mult_expr(_expr<E1> const& u,
              _expr<E2> const& v)
    : _u(u), _v(v) { }

    const Temporary_t & //MatExprType
    eval(const index_t k) const
    {
        // _u.printDetail(gsInfo);
        // _v.printDetail(gsInfo);
        GISMO_ASSERT(0==_u.cols()*_v.rows() || _u.cols() == _v.rows(),
                     "Wrong dimensions "<<_u.cols()<<"!="<<_v.rows()<<" in * operation:\n"
                     << _u <<" times \n" << _v );
        // Note: a * b * c --> (a*b).eval()*c
        tmp = _u.eval(k) * _v.eval(k); return tmp; // assume result not scalarvalued
        //return ( _u.eval(k) * _v.eval(k) ); 
    }
    
    index_t rows() const { return E1::ScalarValued ? _v.rows()  : _u.rows(); }
    index_t cols() const { return E2::ScalarValued ? _u.cols()  : _v.cols(); }
    void setFlag() const { _u.setFlag(); _v.setFlag(); }
    void parse(gsSortedVector<const gsFunctionSet<Scalar>*> & evList) const
    { _u.parse(evList); _v.parse(evList); }

    bool rowSpan() const { return E1::ScalarValued ? _v.rowSpan() : _u.rowSpan(); }
    bool colSpan() const { return E2::ScalarValued ? _u.colSpan() : _v.colSpan(); }
    
    const gsFeVariable<Scalar> & rowVar() const
    { return E1::ScalarValued ? _v.rowVar() : _u.rowVar(); }
    const gsFeVariable<Scalar> & colVar() const
    {
        return E2::ScalarValued ? _u.colVar() : _v.colVar();
        /*
        if ( _v.isScalar() )
            return _u.colVar();
        else
            return _v.colVar();
        */
    }

    void print(std::ostream &os) const { _u.print(os); os<<"*"; _v.print(os); }
};

/**
   Partial specialization for (right) blockwise multiplication
   [A1 A2 A3] * B = [A1*B  A2*B  A3*B]

   mult_v2
*/
template <typename E1, typename E2>
class mult_expr<E1, E2, true> : public _expr<mult_expr<E1, E2, true> >
{
public:
    typedef typename E2::Scalar Scalar;
private:
    typename E1::Nested_t _u;
    typename E2::Nested_t _v;

    mutable gsMatrix<Scalar> res;
public:
    enum {ScalarValued = 0, ColBlocks = 1};
    
    mult_expr(_expr<E1> const& u,
              _expr<E2> const& v)
    : _u(u), _v(v) 
    { 

    }
    
    const gsMatrix<Scalar> & eval(const index_t k) const
    {
        const index_t uc = _u.cols();
        const index_t ur = _u.rows();
        const index_t nb = uc / ur;
        const MatExprType tmpA = _u.eval(k);
        const MatExprType tmpB = _v.eval(k);

        if ( _v.cols() == ur )
        {
            res.resize(ur, uc);
            for (index_t i = 0; i!=nb; ++i)
                res.middleCols(i*ur,ur).noalias()
                    = tmpA.middleCols(i*ur,ur) * tmpB;
        }
        else
        {
            GISMO_ASSERT( 0 == _v.cols() % nb, "Invalid dimensions");
            const index_t vc = _v.cols() / nb;
            res.resize(ur, _v.cols());
            for (index_t i = 0; i!=nb; ++i)
                res.middleCols(i*vc,vc).noalias()
                    = tmpA.middleCols(i*ur,ur) * tmpB.middleCols(i*vc,vc);
        }

        return res;
    }

    index_t rows() const { return _u.rows(); }
    index_t cols() const { return _v.cols() * (_u.cols()/_u.rows()); }
    void setFlag() const { _u.setFlag(); _v.setFlag(); }
    void parse(gsSortedVector<const gsFunctionSet<Scalar>*> & evList) const
    { _u.parse(evList); _v.parse(evList); }

    bool rowSpan() const { return _u.rowSpan(); }
    bool colSpan() const { return false; }
    
    const gsFeVariable<Scalar> & rowVar() const { return _u.rowVar(); }
    const gsFeVariable<Scalar> & colVar() const
    {
        /*
        if ( _v.isScalar() )
            return _u.colVar();
        else
        */
            return _v.colVar();
    }

    void print(std::ostream &os) const
    { os << "("; _u.print(os);os <<"*";_v.print(os);os << ")"; }
};

/// Multiplication operator for expressions
template <typename E1, typename E2> EIGEN_STRONG_INLINE
mult_expr<E1,E2> const operator*(_expr<E1> const& u, _expr<E2> const& v) 
{ return mult_expr<E1, E2>(u, v); }

/**
   Scalar multiplication

   mult_v3
*/
template <typename E2> 
class mult_expr<typename E2::Scalar, E2, false>
    : public _expr<mult_expr<typename E2::Scalar, E2, false> >
// template <typename E> class scmult_expr : public _expr<scmult_expr<E> >
{
public:
    typedef typename E2::Scalar Scalar;
private:
    Scalar const _c;
    typename E2::Nested_t _v;

    //mult_expr(const mult_expr&);
public:
    enum {ScalarValued = E2::ScalarValued, ColBlocks = E2::ColBlocks};
    
    mult_expr(Scalar const & c, _expr<E2> const& v)
    : _c(c), _v(v) { }

    AutoReturn_t eval(const index_t k) const
    {
        return ( _c * _v.eval(k) );
    }

    index_t rows() const { return _v.rows(); }    
    index_t cols() const { return _v.cols(); }
    void setFlag() const { _v.setFlag(); }
    void parse(gsSortedVector<const gsFunctionSet<Scalar>*> & evList) const
    { _v.parse(evList); }

    bool rowSpan() const { return _v.rowSpan(); }
    bool colSpan() const { return _v.colSpan(); }
    
    const gsFeVariable<Scalar> & rowVar() const { return _v.rowVar(); }
    const gsFeVariable<Scalar> & colVar() const { return _v.colVar(); }

    void print(std::ostream &os) const { os << _c <<"*";_v.print(os); }
};

template <typename E2> EIGEN_STRONG_INLINE
mult_expr<typename E2::Scalar,E2,false> const
operator*(typename E2::Scalar const& u, _expr<E2> const& v)
{ return mult_expr<typename E2::Scalar, E2, false>(u, v); }

template <typename E1> EIGEN_STRONG_INLINE
mult_expr<typename E1::Scalar,E1,false> const 
operator*(_expr<E1> const& v, typename E1::Scalar const& u)
{ return mult_expr<typename E1::Scalar,E1, false>(u, v); }
//{ return mult_expr<expr_<typename E1::Scalar>,E1, false>(expr_<typename E1::Scalar>(v), u); }

template <typename E1> EIGEN_STRONG_INLINE
mult_expr<typename E1::Scalar,E1,false> const 
operator-(_expr<E1> const& u)
{ return mult_expr<typename E1::Scalar,E1, false>(-1, u); }
//{ return mult_expr<expr_<typename E1::Scalar>,E1, false>(expr_<typename E1::Scalar>(-1), u); }

/*
template <typename E1> mult_expr<gsMatrix<typename E1::Scalar>,E1,false> const
operator*(gsMatrix<typename E1::Scalar> const& u, _expr<E1> const& v)
{ return mult_expr<gsMatrix<typename E1::Scalar>,E1, false>(u, v); }
*/

/**
   Expression for the Frobenius matrix product, 
   Also block-wise

   [A1 A2 A3] . [B1 B2 B3] = [A1.B1  A2.B2  A3.B3]
*/
template <typename E1, typename E2, bool = E2::ColBlocks>
class frprod_expr : public _expr<frprod_expr<E1, E2> >
{
public:
    typedef typename E1::Scalar Scalar;
    enum {ScalarValued = 0};

private:
    typename E1::Nested_t _u;
    typename E2::Nested_t _v;
    
    mutable gsMatrix<Scalar> res;
    
public:

    frprod_expr(_expr<E1> const& u, _expr<E2> const& v)
    : _u(u), _v(v)
    {
        GISMO_ASSERT(_u.rows() == _v.rows(),
                     "Wrong dimensions "<<_u.rows()<<"!="<<_v.rows()<<" in % operation");
        GISMO_ASSERT(_u.cols() == _v.cols(),
                     "Wrong dimensions "<<_u.cols()<<"!="<<_v.cols()<<" in % operation");
    }
    
    const gsMatrix<Scalar> & eval(const index_t k) const //todo: specialize for nb==1
    {
        // assert _u.size()==_v.size()
        const index_t rb = _u.rows(); //==cb
        const index_t nb = _u.cols() / rb;
        MatExprType A = _u.eval(k);
        MatExprType B = _v.eval(k);
        res.resize(nb, nb);
        for (index_t i = 0; i!=nb; ++i) // all with all
            for (index_t j = 0; j!=nb; ++j)
                res(i,j) =
                (A.middleCols(i*rb,rb).array() * B.middleCols(j*rb,rb).array()).sum();
        return res;
    }

    index_t rows() const { return _u.cols() / _u.rows(); }
    index_t cols() const { return _u.cols() / _u.rows(); }
    void setFlag() const { _u.setFlag(); _v.setFlag(); }
    void parse(gsSortedVector<const gsFunctionSet<Scalar>*> & evList) const
    { _u.parse(evList); _v.parse(evList); }

    bool rowSpan() const { return _u.rowSpan(); }
    bool colSpan() const { return _v.rowSpan(); }
    
    const gsFeVariable<Scalar> & rowVar() const { return _u.rowVar(); }
    const gsFeVariable<Scalar> & colVar() const { return _v.rowVar(); }

    void print(std::ostream &os) const
    { os << "("; _u.print(os); os<<" % "; _v.print(os); os<<")";}
};

/**
   Expression for the Frobenius matrix product, 
   When left hand only side is block-wise

   [A1 A2 A3] . B = [A1.B  A2.B  A3.B]
*/
template <typename E1, typename E2>
class frprod_expr<E1,E2,false> : public _expr<frprod_expr<E1, E2,false> >
{
public:
    typedef typename E1::Scalar Scalar;
    enum {ScalarValued = 0};

private:
    typename E1::Nested_t _u;
    typename E2::Nested_t _v;
    
    mutable gsMatrix<Scalar> res;
    
public:

    frprod_expr(_expr<E1> const& u, _expr<E2> const& v)
    : _u(u), _v(v)
    {
        // GISMO_ASSERT(_u.rows() == _v.rows(),
        //              "Wrong dimensions "<<_u.rows()<<"!="<<_v.rows()<<" in % operation");
        // GISMO_ASSERT(_u.cols() == _v.cols(),
        //              "Wrong dimensions "<<_u.cols()<<"!="<<_v.cols()<<" in % operation");
    }
    
    const gsMatrix<Scalar> & eval(const index_t k) const //todo: specialize for nb==1
    {
        // assert _u.size()==_v.size()
        const index_t rb = _u.rows(); //==cb
        const index_t nb = _u.cols() / rb;
        MatExprType A = _u.eval(k);
        MatExprType B = _v.eval(k);

        res.resize(nb, 1);
        for (index_t i = 0; i!=nb; ++i) // all with all
                res(i,0) =
                (A.middleCols(i*rb,rb).array() * B.array()).sum();
        return res;
    }

    index_t rows() const { return _u.cols() / _u.rows(); }
    index_t cols() const { return _u.cols() / _u.rows(); }
    void setFlag() const { _u.setFlag(); _v.setFlag(); }
    void parse(gsSortedVector<const gsFunctionSet<Scalar>*> & evList) const
    { _u.parse(evList); _v.parse(evList); }

    bool rowSpan() const { return _u.rowSpan(); }
    bool colSpan() const { return _v.rowSpan(); }
    
    const gsFeVariable<Scalar> & rowVar() const { return _u.rowVar(); }
    const gsFeVariable<Scalar> & colVar() const { return _v.rowVar(); }

    void print(std::ostream &os) const
    { os << "("; _u.print(os); os<<" % "; _v.print(os); os<<")";}
};


/// Frobenious product operator for expressions
template <typename E1, typename E2> EIGEN_STRONG_INLINE
frprod_expr<E1,E2> const  operator%(_expr<E1> const& u, _expr<E2> const& v)
{ return frprod_expr<E1, E2>(u, v); }

/**
   Expression for scalar division operation
 */
template <typename E1, typename E2>
class divide_expr : public _expr<divide_expr<E1,E2> >
{
    typename E1::Nested_t _u;
    typename E2::Nested_t _v;

public:
    typedef typename E1::Scalar Scalar;

public:
    enum {ScalarValued = E1::ScalarValued};
    
    divide_expr(_expr<E1> const& u, _expr<E2> const& v)
    : _u(u), _v(v)
    {
        GISMO_STATIC_ASSERT(E2::ScalarValued, "The denominator needs to be scalar valued.");
    }

    AutoReturn_t eval(const index_t k) const
    { return ( _u.eval(k) / _v.eval(k) ); }
    
    index_t rows() const { return _u.rows(); }    
    index_t cols() const { return _u.cols(); }
    void setFlag() const { _u.setFlag(); _v.setFlag();}
    void parse(gsSortedVector<const gsFunctionSet<Scalar>*> & evList) const
    { _u.parse(evList); _v.parse(evList); }

    bool rowSpan() const { return _u.rowSpan(); }
    bool colSpan() const { return _u.colSpan(); }
    
    const gsFeVariable<Scalar> & rowVar() const { return _u.rowVar(); }
    const gsFeVariable<Scalar> & colVar() const { return _u.colVar(); }

    void print(std::ostream &os) const
    { os << "("; _u.print(os);os <<" / ";_v.print(os);os << ")"; }
};

/**
   Division specialization for constant value denominator
 */
template <typename E1>
class divide_expr<E1,typename E1::Scalar>
    : public _expr<divide_expr<E1,typename E1::Scalar> >
{
public:
    typedef typename E1::Scalar Scalar;

private:
    typename E1::Nested_t _u;
    Scalar  const   _c;

public:
    enum {ScalarValued = E1::ScalarValued};
    
    divide_expr(_expr<E1> const& u, Scalar const  c)
    : _u(u), _c(c) { }

    AutoReturn_t eval(const index_t k) const
    { return ( _u.eval(k) / _c ); }

    index_t rows() const { return _u.rows(); }    
    index_t cols() const { return _u.cols(); }
    void setFlag() const { _u.setFlag(); }
    void parse(gsSortedVector<const gsFunctionSet<Scalar>*> & evList) const
    { _u.parse(evList); }

    bool rowSpan() const { return _u.rowSpan(); }
    bool colSpan() const { return _u.colSpan(); }

    void print(std::ostream &os) const
    { os << "("; _u.print(os);os <<"/"<< _c << ")"; }
};

/**
   Division specialization for constant value numerator
 */
template <typename E2>
class divide_expr<typename E2::Scalar,E2>
    : public _expr<divide_expr<typename E2::Scalar,E2> >
{
public:
    typedef typename E2::Scalar Scalar;

private:
    Scalar  const   _c;
    typename E2::Nested_t _u;
public:
    enum {ScalarValued = 1};
    
    divide_expr(Scalar const c, _expr<E2> const& u)
    : _c(c), _u(u)
    { GISMO_STATIC_ASSERT(E2::ScalarValued, "The denominator needs to be scalar valued."); }

    Scalar eval(const index_t k) const
    { return ( _c / _u.val().eval(k) ); }

    index_t rows() const { return 0; }    
    index_t cols() const { return 0; }
    void setFlag() const { _u.setFlag(); }
    void parse(gsSortedVector<const gsFunctionSet<Scalar>*> & evList) const
    { _u.parse(evList); }

    bool rowSpan() const { return false; }
    bool colSpan() const { return false; }

    void print(std::ostream &os) const
    { os << "("<< _c <<"/";_u.print(os);os << ")";}
};

/// Scalar division operator for expressions
template <typename E1, typename E2> EIGEN_STRONG_INLINE
divide_expr<E1,E2> const operator/(_expr<E1> const& u, _expr<E2> const& v)
{ return divide_expr<E1,E2>(u, v); }

template <typename E> EIGEN_STRONG_INLINE
divide_expr<E,typename E::Scalar> const
operator/(_expr<E> const& u, const typename E::Scalar v)
{ return divide_expr<E,typename E::Scalar>(u, v); }

template <typename E> EIGEN_STRONG_INLINE
divide_expr<typename E::Scalar,E> const
operator/(const typename E::Scalar u, _expr<E> const& v)
{ return divide_expr<typename E::Scalar,E>(u, v); }

/**
   Expression for addition operation
 */
template <typename E1, typename E2>
class add_expr : public _expr<add_expr<E1, E2> >
{
    typename E1::Nested_t _u;
    typename E2::Nested_t _v;

public:
    enum {ScalarValued = E1::ScalarValued && E2::ScalarValued,
          ColBlocks = E1::ColBlocks && E2::ColBlocks };
    
    typedef typename E1::Scalar Scalar;
    
    add_expr(_expr<E1> const& u, _expr<E2> const& v)
    : _u(u), _v(v)
    {
        // false on c++98 ?
        //GISMO_STATIC_ASSERT((int)(E1::ColBlocks)==(int)(E2::ColBlocks), "Cannot add if the number of colums do not agree.");
        GISMO_ASSERT((int)E1::ColBlocks == (int)E2::ColBlocks, "Error: "<< E1::ColBlocks
                     <<"!="<< E2::ColBlocks);
    }

    AutoReturn_t eval(const index_t k) const
    {
        // (?) Wrong dimensions 0!=1 in + operation ?
        GISMO_ASSERT(_u.rows() == _v.rows(),
                     "Wrong dimensions "<<_u.rows()<<"!="<<_v.rows()<<" in + operation:\n"
                     << _u <<" plus \n" << _v );
        GISMO_ASSERT(_u.cols() == _v.cols(),
                     "Wrong dimensions "<<_u.cols()<<"!="<<_v.cols()<<" in + operation:\n"
                     << _u <<" plus \n" << _v );
        return _u.eval(k) + _v.eval(k);
    }
    
    index_t rows() const { return _u.rows(); }    
    index_t cols() const { return _u.cols(); }
    void setFlag() const { _u.setFlag(); _v.setFlag(); }
    void parse(gsSortedVector<const gsFunctionSet<Scalar>*> & evList) const
    { _u.parse(evList); _v.parse(evList); }

    bool rowSpan() const { return _u.rowSpan(); }
    bool colSpan() const { return _u.colSpan(); }
    
    const gsFeVariable<Scalar> & rowVar() const { return _u.rowVar(); }
    const gsFeVariable<Scalar> & colVar() const { return _v.colVar(); }

    void print(std::ostream &os) const
    { os << "("; _u.print(os);os <<" + ";_v.print(os);os << ")"; }
};

/// Addition operator for expressions
template <typename E1, typename E2> EIGEN_STRONG_INLINE
add_expr<E1,E2> const operator+(_expr<E1> const& u, _expr<E2> const& v)
{ return add_expr<E1, E2>(u, v); }


/*// testing, |, ^, &, <<, >>, ||, &&,  unary ~
template <typename E1, typename E2> add_expr<E1,E2> const
operator|(_expr<E1> const& u, _expr<E2> const& v)
{ return add_expr<E1, E2>(u, v); }
template <typename E1, typename E2> add_expr<E1,E2> const
operator^(_expr<E1> const& u, _expr<E2> const& v)
{ return add_expr<E1, E2>(u, v); }
*/


/**
   Expression for (square) matrix summation operation
   
   M [r x r*k] is a list of matrices 
   Summation is done over k,
   [M1 M2 .. Mk]

   u [s x k] is a list of vectors

   Computed quantity is of size [r x r*s] and contains
   [ ... sum_k(Mk * u(s,k) ) ... ]_s
 */
template <typename E1, typename E2>
class summ_expr : public _expr<summ_expr<E1,E2> >
{
public:
    typedef typename E1::Scalar Scalar;

    enum {ScalarValued = 0, ColBlocks = E2::ColBlocks};
    
    summ_expr(E1 const& u, E2 const& M) : _u(u), _M(M) { }

    const gsMatrix<Scalar> & eval(const index_t k) const
    {
        MatExprType sl = _u.eval(k);
        MatExprType ml = _M.eval(k);
        const index_t sr  = sl.rows();
        const index_t mr  = ml.rows();
        const index_t mb  = ml.cols() / mr;

        res.setZero(mr, sr * mr);
        for (index_t i = 0; i!=sr; ++i)
            for (index_t t = 0; t!=mb; ++t)
                res.middleCols(i*mr,mr) += sl(i,t) * ml.middleCols(t*mr,mr);
        return res;
    }

    index_t rows() const { return _M.rows(); }    
    index_t cols() const { return _u.rows() * _M.rows(); }
    void setFlag() const { _u.setFlag(); _M.setFlag(); }
    void parse(gsSortedVector<const gsFunctionSet<Scalar>*> & evList) const
    { _u.parse(evList); _M.parse(evList); }

    bool rowSpan() const { return _u.rowSpan(); }
    bool colSpan() const { return false; }
    
    const gsFeVariable<Scalar> & rowVar() const { return _u.rowVar(); }
    const gsFeVariable<Scalar> & colVar() const { GISMO_ERROR("summ::colVar()"); }

    void print(std::ostream &os) const
    { os << "sum("; _M.print(os); os<<","; _u.print(os); os<<")"; }

private:
    typename E1::Nested_t _u;
    typename E2::Nested_t _M;
    
    mutable gsMatrix<Scalar> res;
};

/// Summation operator for expressions
template <typename E1, typename E2> EIGEN_STRONG_INLINE
summ_expr<E1,E2> const summ(E1 const & u, E2 const& M)
{ return summ_expr<E1,E2>(u, M); }

/**
   Expression for subtraction operation
 */
template <typename E1, typename E2>
class sub_expr : public _expr<sub_expr<E1, E2> >
{
    typename E1::Nested_t _u;
    typename E2::Nested_t _v;

public:
    enum {ScalarValued = E1::ScalarValued && E2::ScalarValued,
          ColBlocks = E1::ColBlocks && E2::ColBlocks };

    typedef typename E1::Scalar Scalar;
    
    sub_expr(_expr<E1> const& u, _expr<E2> const& v)
    : _u(u), _v(v)
    {
        GISMO_STATIC_ASSERT((int)E1::ColBlocks == (int)E2::ColBlocks, "Cannot subtract if the number of colums do not agree.");
        /* // Note: rows()/cols() might be still pending at construction time..
        GISMO_ASSERT(_u.rows() == _v.rows(),
                     "Wrong dimensions "<<_u.rows()<<"!="<<_v.rows()<<" in - operation");
        GISMO_ASSERT(_u.cols() == _v.cols(),
                     "Wrong dimensions "<<_u.cols()<<"!="<<_v.cols()<<" in - operation");
        */
    }

    AutoReturn_t eval(const index_t k) const
    {
        GISMO_ASSERT(_u.rows() == _v.rows(),
                     "Wrong dimensions "<<_u.rows()<<"!="<<_v.rows()<<" in - operation");
        GISMO_ASSERT(_u.cols() == _v.cols(),
                     "Wrong dimensions "<<_u.cols()<<"!="<<_v.cols()<<" in - operation");
        return (_u.eval(k) - _v.eval(k) );
    }
    
    index_t rows() const { return _u.rows(); }    
    index_t cols() const { return _u.cols(); }
    void setFlag() const { _u.setFlag(); _v.setFlag(); }
    void parse(gsSortedVector<const gsFunctionSet<Scalar>*> & evList) const
    { _u.parse(evList); _v.parse(evList); }

    bool rowSpan() const { return _u.rowSpan(); }
    bool colSpan() const { return _v.colSpan(); }
    
    const gsFeVariable<Scalar> & rowVar() const { return _u.rowVar(); }
    const gsFeVariable<Scalar> & colVar() const { return _v.colVar(); }

    void print(std::ostream &os) const
    { os << "("; _u.print(os); os<<" - ";_v.print(os); os << ")";}
};

/// Subtraction operator for expressions
template <typename E1, typename E2> EIGEN_STRONG_INLINE
sub_expr<E1,E2> const operator-(_expr<E1> const& u, _expr<E2> const& v)
{ return sub_expr<E1, E2>(u, v); }

template <typename E2> EIGEN_STRONG_INLINE
sub_expr<_expr<typename E2::Scalar>,E2> const
operator-(typename E2::Scalar const& s, _expr<E2> const& v)
{
    // assert E2::ScalarValued
    return sub_expr<_expr<typename E2::Scalar>, E2>(_expr<typename E2::Scalar>(s), v);
}

//template <typename E2> EIGEN_STRONG_INLINE
//sub_expr<E1,E2> const operator-(typename E2::Scalar const& u, _expr<E2> const& v)
//{ return sub_expr<E1, E2>(u, v); }

/**
   Expression for symmetrization operation
 */
template <typename E>
class symm_expr : public _expr<symm_expr<E> >
{
    typename E::Nested_t _u;

    mutable gsMatrix<typename E::Scalar> tmp;
public:
    typedef typename E::Scalar Scalar;
    
    symm_expr(_expr<E> const& u)
    : _u(u) { }

    MatExprType eval(const index_t k) const
    {
        //const MatExprType tmp = _u.eval(k);
        tmp = _u.eval(k);
        // todo: avoid temporary or not ?
        return tmp * tmp.transpose();
    }

    index_t rows() const { return _u.rows(); }    
    index_t cols() const { return _u.rows(); }
    void setFlag() const { _u.setFlag(); }
    void parse(gsSortedVector<const gsFunctionSet<Scalar>*> & evList) const
    { _u.parse(evList); }

    bool rowSpan() const { return _u.rowSpan(); }
    bool colSpan() const { return _u.rowSpan(); }
    
    const gsFeVariable<Scalar> & rowVar() const { return _u.rowVar(); }
    const gsFeVariable<Scalar> & colVar() const { return _u.rowVar(); }

    void print(std::ostream &os) const { os << "symm("; _u.print(os); os <<")"; }
};

/* Symmetrization operation
template <typename E> symm_expr<E> const
symm(_expr<E> const& u) { return symm_expr<E>(u);}
*/

#undef MatExprType
#undef AutoReturn_t

#if(__cplusplus >= 201402L)

// Shortcuts for common quantities, for instance function
// transformations by the geometry map \a G
#define GISMO_SHORTCUT_VAR_EXPRESSION(name,impl) template<class T> EIGEN_STRONG_INLINE \
auto name(const gsFeVariable<T> & u) { return impl; }
#define GISMO_SHORTCUT_MAP_EXPRESSION(name,impl) template<class T> EIGEN_STRONG_INLINE \
auto name(const gsGeometryMap<T> & G) { return impl; }
#define GISMO_SHORTCUT_PHY_EXPRESSION(name,impl) template<class T> EIGEN_STRONG_INLINE \
auto name(const gsFeVariable<T> & u, const gsGeometryMap<T> & G) { return impl; }

// Divergence
GISMO_SHORTCUT_VAR_EXPRESSION(  div, jac(u  ).trace()    )

// The unit (normalized) boundary (outer pointing) normal
GISMO_SHORTCUT_MAP_EXPRESSION(unv, nv(G).normalized()   ) //(!) bug + mem. leak

GISMO_SHORTCUT_PHY_EXPRESSION(igrad, grad(u)*jac(G).ginv() ) // transpose() problem ??
GISMO_SHORTCUT_PHY_EXPRESSION( ijac, jac(u) * jac(G).ginv())
GISMO_SHORTCUT_PHY_EXPRESSION( idiv, ijac(u,G).trace()    )
GISMO_SHORTCUT_PHY_EXPRESSION(ihess,
jac(G).ginv().tr()*( hess(u) - summ(igrad(u,G),hess(G)) ) * jac(G).ginv() )
GISMO_SHORTCUT_PHY_EXPRESSION(ilapl, ihess(u,G).trace()   )

#else

#define GISMO_SHORTCUT_VAR_EXPRESSION(name,impl) \
name(const gsFeVariable<T> & u) { return impl; }
#define GISMO_SHORTCUT_MAP_EXPRESSION(name,impl) \
name(const gsGeometryMap<T> & G) { return impl; }
#define GISMO_SHORTCUT_PHY_EXPRESSION(name,impl)  \
name(const gsFeVariable<T> & u, const gsGeometryMap<T> & G) { return impl; }

template<class T> EIGEN_STRONG_INLINE trace_expr<jac_expr<T> >
GISMO_SHORTCUT_VAR_EXPRESSION(div, jac(u).trace() )

template<class T> EIGEN_STRONG_INLINE normalized_expr<onormal_expr<T> >
GISMO_SHORTCUT_MAP_EXPRESSION(unv, nv(G).normalized() )

template<class T> EIGEN_STRONG_INLINE mult_expr<grad_expr<T>,jacGinv_expr<T>, 0>
GISMO_SHORTCUT_PHY_EXPRESSION(igrad, grad(u)*jac(G).ginv())

template<class T> EIGEN_STRONG_INLINE mult_expr<jac_expr<T>,jacGinv_expr<T>, 1>
GISMO_SHORTCUT_PHY_EXPRESSION(ijac, jac(u) * jac(G).ginv() )

template<class T> EIGEN_STRONG_INLINE trace_expr<mult_expr<jac_expr<T>,jacGinv_expr<T>, 1> >
GISMO_SHORTCUT_PHY_EXPRESSION(idiv, ijac(u,G).trace() )

template<class T> EIGEN_STRONG_INLINE mult_expr<mult_expr<tr_expr<jacGinv_expr<T> >,sub_expr<hess_expr<T>,summ_expr<mult_expr<grad_expr<T>, jacGinv_expr<T>, 0>, hess_expr<T> > >, 0>, jacGinv_expr<T>, 1>
GISMO_SHORTCUT_PHY_EXPRESSION(ihess, jac(G).ginv().tr()*(hess(u)-summ(igrad(u,G),hess(G)))*jac(G).ginv() )

//#define ilapl(u,G) ihess(u,G).trace()
template<class T> EIGEN_STRONG_INLINE trace_expr< mult_expr<mult_expr<tr_expr<jacGinv_expr<T> >, sub_expr<hess_expr<T>, summ_expr<mult_expr<grad_expr<T>, jacGinv_expr<T>, 0>, hess_expr<T> > >, 0>, jacGinv_expr<T>, 1> >
GISMO_SHORTCUT_PHY_EXPRESSION(ilapl, ihess(u,G).trace() )

#endif
#undef GISMO_SHORTCUT_PHY_EXPRESSION
#undef GISMO_SHORTCUT_VAR_EXPRESSION
#undef GISMO_SHORTCUT_MAP_EXPRESSION

} // namespace expr

} //namespace gismo
