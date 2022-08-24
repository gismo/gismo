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
#include <gsAssembler/gsDirichletValues.h>

namespace gismo
{

// Adaptor to compute Hessian
template<typename T>
void secDerToHessian(typename gsMatrix<T>::constRef & secDers,
                     const index_t dim,
                     Eigen::Matrix<T,Dynamic,Dynamic> & hessian)
{
    const index_t sz = dim*(dim+1)/2;
    const gsAsConstMatrix<T> ders(secDers.data(), sz, secDers.size() / sz );
    hessian.resize(dim*dim, ders.cols() );

    switch ( dim )
    {
    case 1:
        hessian = secDers; // ders
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

// Forward declaration in gismo namespace
template<class T> class gsExprHelper;

/** @namespace gismo::expr

    @brief
    This namespace contains expressions used for FE computations

    \ingroup Assembler
*/
namespace expr
{

#if __cplusplus >= 201402L || _MSVC_LANG >= 201402L // c++14
#  define MatExprType  auto
#  define AutoReturn_t auto
#  define GS_CONSTEXPR constexpr
//note: in c++11 auto-return requires -> decltype(.)
#else // 199711L, 201103L
#  define MatExprType typename gsMatrix<Scalar>::constRef
#  define AutoReturn_t typename util::conditional<ScalarValued,Scalar,MatExprType>::type
#  define GS_CONSTEXPR
#endif

template<class T> class gsFeSpace;
template<class T> class gsFeVariable;
template<class T> class gsFeSolution;
template<class E> class symm_expr;
template<class E> class symmetrize_expr;
template<class E> class trace_expr;
template<class E> class adjugate_expr;
template<class E> class norm_expr;
template<class E> class sqNorm_expr;
template<class E> class det_expr;
template<class E> class value_expr;
template<class E> class asdiag_expr;
template<class E> class rowsum_expr;
template<class E> class colsum_expr;
template<class E> class col_expr;
template<class T> class meas_expr;
template<class E> class inv_expr;
template<class E> class tr_expr;
template<class E> class sign_expr;
template<class T> class cdiam_expr;
template<class E> class temp_expr;
template<class E1, class E2, bool = E1::ColBlocks> class mult_expr
{using E1::GISMO_ERROR_mult_expr_has_invalid_template_arguments;};

/*
 Traits class for expressions
*/
template <typename E> struct expr_traits
{
public:
//    typedef typename E::Scalar Scalar;
  typedef real_t Scalar;//todo
  typedef const E Nested_t;
};

template <class E> struct is_arithmetic{enum{value=0};};
template <> struct is_arithmetic<real_t>{enum{value=1};};

template <typename E, bool = is_arithmetic<E>::value >
class _expr {using E::GISMO_ERROR_expr;};

/**
 \brief Base class for all expressions
*/
template <typename E>
class _expr<E, false>
{
protected://private:
    _expr(){}
    _expr(const _expr&) { }
    //friend E;//
public:

    enum {ScalarValued = 0, ColBlocks = 0};
    //todo: ValueType=0,1,2 (scalar,vector,matrix)
    //Space = 0,1,2 (instead of rowSpan ColSpan, RowCol ?)
    enum {Space = 0};

    //typedef typename E::Nested_t Nested_t;
    //typedef typename E::Scalar   Scalar;
    typedef typename expr_traits<E>::Nested_t Nested_t;
    typedef typename expr_traits<E>::Scalar   Scalar;

  /// Prints the expression as a string to \a os
  void print(std::ostream &os) const
  {
      static_cast<E const&>(*this).print(os);
  /*
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
  */
  }

  std::ostream & printDetail(std::ostream &os) const
  {
      os << (isScalar() ? "Scalar " :
            (isVector() ? "Vector " :
            (isMatrix() ? "Matrix " :
             "Unknown ") ) )
         <<"expression of size "<< rows() // bug: this is not fixed, or may not be known
         << " x "<<cols()<<"\n";
      print(os);
      return os;
  }

  /// Evaluates the expression at evaluation point indexed by \a k
  MatExprType eval(const index_t k) const
    { return static_cast<E const&>(*this).eval(k); }

    /// Returns the transpose of the expression
    tr_expr<E> tr() const
    { return tr_expr<E>(static_cast<E const&>(*this)); }

    /// Returns the sign of the expression
    sign_expr<E> sgn() const
    { return sign_expr<E>(static_cast<E const&>(*this)); }

    /// Returns an evaluation of the (sub-)expression in temporary memory
    temp_expr<E> temp() const
    { return temp_expr<E>(static_cast<E const&>(*this)); }

    /// Returns the inverse of the expression (for matrix-valued expressions)
    inv_expr<E> const inv() const
    { return inv_expr<E>(static_cast<E const&>(*this)); }

    /// Returns the trace of the expression (for matrix-valued expressions)
    trace_expr<E> trace() const
    { return trace_expr<E>(static_cast<E const&>(*this)); }

    /// Returns the adjugate of the expression (for matrix-valued expressions)
    adjugate_expr<E> adj() const
    { return adjugate_expr<E>(static_cast<E const&>(*this)); }

    /// Returns the Euclidean norm of the expression
    norm_expr<E> norm() const
    { return norm_expr<E>(static_cast<E const&>(*this)); }

    /// Returns the determinant of the expression
    det_expr<E> det() const
    { return det_expr<E>(static_cast<E const&>(*this)); }

    /// Returns the squared Euclidean norm of the expression
    sqNorm_expr<E> sqNorm() const
    { return sqNorm_expr<E>(static_cast<E const&>(*this)); }

    /// Returns the square root of the expression (component-wise)
    mult_expr<E,E,0> (sqr)() const { return (*this)*(*this); }

    symm_expr<E> symm() const
    { return symm_expr<E>(static_cast<E const&>(*this)); }

    symmetrize_expr<E> symmetrize() const
    { return symmetrize_expr<E>(static_cast<E const&>(*this)); }

    /// For matrix-valued expressions which are actually 1x1 matrix,
    /// returns a scalar valued expression
    value_expr<E> val() const
    { return value_expr<E>(static_cast<E const&>(*this)); }

    /// Returns a diagonal matrix expression of the vector expression
    asdiag_expr<E> asDiag() const
    { return asdiag_expr<E>(static_cast<E const&>(*this)); }

    /// Returns the rowSum of a matrix
    rowsum_expr<E> rowSum() const
    { return rowsum_expr<E>(static_cast<E const&>(*this)); }

    /// Returns the colSum of a matrix
    colsum_expr<E> colSum() const
    { return colsum_expr<E>(static_cast<E const&>(*this)); }

    col_expr<E> operator[](const index_t i) const
    { return col_expr<E>(static_cast<E const&>(*this),i); }

    /// Returns the row-size of the expression
    index_t rows() const
    { return static_cast<E const&>(*this).rows(); }

    /// Returns the column-size of the expression
    index_t cols() const
    { return static_cast<E const&>(*this).cols(); }

    index_t cardinality() const
    { return static_cast<E const&>(*this).cardinality_impl(); }

    static index_t cardinality_impl() { return 1; }
    //static index_t cardinality_impl() { return static_cast<E const&>(*this).rowVar().cardinality_impl(); }

    ///\brief Returns true iff the expression is scalar-valued.
    /// \note This is a runtime check, for compile-time check use E::ScalarValued
    bool isScalar() const { return rows()*cols()<=1; } //!rowSpan && !colSpan

    ///\brief Returns true iff the expression is vector-valued.
    static bool isVector() { return E::rowSpan && (!E::colSpan); }

    ///\brief Returns true iff the expression is matrix-valued.
    static bool isMatrix() { return E::rowSpan && E::colSpan; }

    //enum{rowSpan = 0, colSpan = 0};
    // static GS_CONSTEXPR bool rowSpan() { return E::rowSpan(); }
    // static bool colSpan { return E::colSpan; }

    ///\brief Sets the required evaluation flags
    void setFlag() const { static_cast<E const&>(*this).setFlag(); }

    ///\brief Parse the expression and discover the list of evaluation
    ///sources, also sets the required evaluation flags
    void parse(gsSortedVector<const gsFunctionSet<Scalar>*> & evList) const
    { static_cast<E const&>(*this).parse(evList); }

    /// Returns the space that is found on the left-most of the
    /// expression
    const gsFeSpace<Scalar> & rowVar() const
    {
        // assert ValueType!=0
        return static_cast<E const&>(*this).rowVar();
    }

    /// Returns the space that is found on the right-most of
    /// the expression
    const gsFeSpace<Scalar> & colVar() const
    {
        // assert ValueType==2
        return static_cast<E const&>(*this).colVar();
    }

    // Overload conversions, eg. converts _expr<mult_expr> to
    // mult_expr.
    operator E&()             { return static_cast<      E&>(*this); }
    operator E const&() const { return static_cast<const E&>(*this); }

    //E const & derived() const { return static_cast<const E&>(*this); }
};

/// Stream operator for expressions
template <typename E>
std::ostream &operator<<(std::ostream &os, const _expr<E> & b)
{b.print(os); return os; }

/*
   Null expression is a compatibility expression invalid at runtime
*/
template<class T>
class gsNullExpr : public _expr<gsNullExpr<T> >
{
public:

    operator const gsFeSpace<T> & () const
    {
        static gsFeSpace<T> vv;
        return vv;
    }

    typedef T Scalar;
    gsMatrix<T> eval(const index_t) const { GISMO_ERROR("gsNullExpr"); }
    inline index_t rows() const { GISMO_ERROR("gsNullExpr"); }
    inline index_t cols() const { GISMO_ERROR("gsNullExpr"); }
    inline void setFlag() const {/* gsInfo<<"gsNullExpr emtpy flag\n"; */ }
    void parse(gsSortedVector<const gsFunctionSet<T>*> &) const { }

    const gsFeSpace<T> & rowVar() const { GISMO_ERROR("gsNullExpr"); }
    const gsFeSpace<T> & colVar() const { GISMO_ERROR("gsNullExpr"); }

    enum{rowSpan = 0, colSpan = 0};

    void print(std::ostream &os) const { os << "NullExpr"; }

    static const gsNullExpr & get()
    {
        static gsNullExpr o;
        return o;
    }
//private:
    gsNullExpr() {}
};


template<class E>
class symbol_expr : public _expr<E>
{
public:
    typedef typename expr_traits<E>::Scalar Scalar;

    friend class gismo::gsExprHelper<Scalar>;
protected:
    //const gsFuncData<Scalar>    * m_fd2; // more data when needed
    const gsFunctionSet<Scalar> * m_fs; ///< Evaluation source for this FE variable
    const gsFuncData<Scalar>    * m_fd; ///< Temporary variable storing flags and evaluation data
    index_t m_d;                   ///< Dimension of this (scalar or vector) variable
    const gsMapData<Scalar>     * m_md; ///< If set, the variable is composed with a geometry map
    // comp(u,G)

public:

    enum {Space = 4};// remove!

    /// Returns the function source
    const gsFunctionSet<Scalar> & source() const {return *m_fs;}

    /// Returns the function data
    const gsFuncData<Scalar> & data() const {return *m_fd;}

    /// Returns the mapping data (precondition: composed()==true)
    const gsMapData<Scalar> & mapData() const {return *m_md;}

    /// Returns true if the variable is a composition
    bool composed() const {return NULL!=m_md;}

    index_t cardinality_impl() const { return m_d * m_fd->actives.rows(); }

private:

    void setSource(const gsFunctionSet<Scalar> & fs) { m_fs = &fs;}
    void setData(const gsFuncData<Scalar> & val) { m_fd = &val;}
    void clear() { m_fs = NULL; }
    // gsFuncData<Scalar> & data() {return *m_fd;}
    // gsMapData<Scalar> & mapData() {return *m_md;}

protected:

    explicit symbol_expr(index_t _d) : m_fs(NULL), m_fd(NULL), m_d(_d), m_md(NULL) { }

    void registerData(const gsFunctionSet<Scalar> & fs, const gsFuncData<Scalar> & val, index_t d)
    {
        GISMO_ASSERT(NULL==m_fs, "gsFeVariable: already registered");
        m_fs = &fs ;
        m_fd = &val;
        m_d  = d;
        m_md = NULL;
    }

    void registerData(const gsFunctionSet<Scalar> & fs, const gsFuncData<Scalar> & val, index_t d,
                      const gsMapData<Scalar> & md)
    {
        registerData(fs,val,d);
        m_md  = &md;
    }

public:

    bool isValid() const { return NULL!=m_fd && NULL!=m_fs; }

    // component
    // expr comp(const index_t i) const { return comp_expr<Scalar>(*this,i); }
    // eval(k).col(i)


    // The evaluation return rows for (basis) functions and columns
    // for (coordinate) components
    MatExprType eval(const index_t k) const
    { return m_fd->values[0].col(k).blockDiag(m_d); } //!!
    //{ return m_fd->values[0].col(k); }

    const gsFeSpace<Scalar> & rowVar() const {return gsNullExpr<Scalar>::get();}
    const gsFeSpace<Scalar> & colVar() const {return gsNullExpr<Scalar>::get();}

    enum{rowSpan = 0, colSpan = 0};

    index_t rows() const
    {
        GISMO_ASSERT(NULL!=m_fs, "FeVariable: Function member not registered");
        GISMO_ASSERT(NULL!=m_fd, "FeVariable: FuncData member not registered");
        return m_fs->targetDim();
        /*
        gsDebugVar(m_fd->flags & NEED_VALUE);
        gsDebugVar(m_fd->values[0].rows());
        gsDebugVar(m_fd->flags & NEED_ACTIVE);
        gsDebugVar(m_fd->actives.rows());
        gsDebugVar(m_fd->flags & NEED_DERIV);
        */

	/* HV: NOT NEDED SINCE CARDINALITY IMPLEMENTATION
        // note: precomputation is needed
        if (m_fd->flags & NEED_VALUE)
        {return m_d * m_fd->values[0].rows();}
        if (m_fd->flags & NEED_ACTIVE) // note: gsFunction coeff ??
        {return m_d * m_fd->actives.rows();}
        if (m_fd->flags & NEED_DERIV)
        {return m_d * m_fd->values[0].rows();}
        GISMO_ERROR("Cannot deduce row size.");
	*/
    }

    index_t cols() const { return m_d; }

    void setFlag() const
    {
        GISMO_ASSERT(NULL!=m_fd, "FeVariable: FuncData member not registered");
        m_fd->flags |= NEED_VALUE | NEED_ACTIVE;
        if (NULL!=m_md) m_md->flags |= NEED_VALUE;
    }

    void parse(gsSortedVector<const gsFunctionSet<Scalar>*> & evList) const
    {
        GISMO_ASSERT(NULL!=m_fd, "FeVariable: FuncData member not registered");
        evList.push_sorted_unique(m_fs);
        m_fd->flags |= NEED_VALUE;
        if (NULL!=m_md) m_md->flags |= NEED_VALUE;
    }

    void print(std::ostream &os) const { os << "u"; }

public:

    /// Returns the vector dimension of the FE variable
    index_t dim() const { return m_d;}

    /// Returns the vector dimension of the FE variable
    index_t targetDim() const { return m_d;}

    /// Returns the parameter domain dimension the FE variable
    index_t parDim() const
    {
        return m_fs->domainDim();
        //return m_fd->dim.first;
    }

    index_t cSize()  const
    {
        GISMO_ASSERT(0!=m_fd->values[0].size(),"Probable error.");
        return m_fd->values[0].rows();
    } // coordinate size
};

/*
    Column expression
*/
template<class E>
class col_expr : public _expr<col_expr<E> >
{
    typename E::Nested_t _c;
    const index_t _i;
public:
    typedef typename E::Scalar Scalar;
    typedef const col_expr<E> Nested_t;

    col_expr(const E & c, const index_t i) : _c(c), _i(i) { }

public:

    //ConstColXpr
    inline MatExprType eval(const index_t k) const { return _c.eval(k).col(_i); }

    index_t rows() const { return _c.rows(); }
    index_t cols() const { return 1; }
    void setFlag() const { _c.setFlag();}
    void parse(gsSortedVector<const gsFunctionSet<Scalar>*> & evList) const { _c.parse(evList); }

    enum{rowSpan = E::rowSpan, colSpan = 0};

    const gsFeSpace<Scalar> & rowVar() const { return _c.rowVar(); }
    const gsFeSpace<Scalar> & colVar() const { return _c.colVar(); }

    void print(std::ostream &os) const { os<<_c<<"["<<_i<<"]"; }
};

/*
    Expression for a constant value
 */
template<class T>
class _expr<T, true> : public _expr<_expr<T> >
{
    const T _c;
public:
    typedef T Scalar;
    typedef const _expr<T> Nested_t;

    _expr(const Scalar & c) : _c(c) { }

public:
    enum {ScalarValued = 1};
    enum {Space = 0};

    inline Scalar eval(const index_t ) const { return _c; }

    inline _expr val() const { return *this; }
    index_t rows() const { return 0; }
    index_t cols() const { return 0; }
    void setFlag() const { }
    void parse(gsSortedVector<const gsFunctionSet<Scalar>*> &) const { }

    enum{rowSpan = 0, colSpan = 0};
    const gsFeSpace<T> & rowVar() const { return gsNullExpr<T>(); }
    const gsFeSpace<T> & colVar() const { return gsNullExpr<T>(); }

    void print(std::ostream &os) const { os<<_c; }
};

/*
   Geometry map expression
 */
template<class T>
class gsGeometryMap : public _expr<gsGeometryMap<T> >
{
    const gsFunctionSet<T> * m_fs; ///< Evaluation source for this geometry map
    const gsMapData<T> *  m_fd;    ///< Temporary variable storing flags and evaluation data
    //index_t d, n;

public:

    enum {Space = 3};

    /// Returns the function source
    const gsFunctionSet<T> & source() const {return *m_fs;}

    /// Returns the function data
    const gsMapData<T> & data() const  { return *m_fd; }

    index_t targetDim() const { return m_fs->targetDim();}
public:
    typedef T Scalar;

    friend class gismo::gsExprHelper<T>;

    void print(std::ostream &os) const { os << "G"; }

    MatExprType eval(const index_t k) const { return m_fd->values[0].col(k); }

    void setFlag() const
    {
        GISMO_ASSERT(NULL!=m_fd, "GeometryMap not registered");
        m_fd->flags |= NEED_VALUE;
    }

protected:

    gsGeometryMap() : m_fs(NULL), m_fd(NULL) { }

    /// Registers the source function and evaluation data
    void registerData(const gsFunctionSet<T> & fs, const gsMapData<T> & val)
    {
        m_fs = &fs;
        m_fd = &val;
    }

public:

    /// Returns true iff the source function has been set
    bool isValid() const { return NULL!=m_fs; }

    index_t rows() const { return m_fd->dim.second; }
    index_t cols() const { return 1; }

    enum{rowSpan = 0, colSpan = 0};

    const gsFeSpace<T> & rowVar() const { return gsNullExpr<T>(); }
    const gsFeSpace<T> & colVar() const { return gsNullExpr<T>(); }

    void parse(gsSortedVector<const gsFunctionExpr<T>*> & evList) const
    {
        GISMO_ASSERT(NULL!=m_fd, "GeometryMap not registered");
        evList.push_unique(m_fs);
        m_fd->flags |= NEED_VALUE;
    }
};

// Traits for gsGeometryMap
template <typename T>  struct expr_traits<gsGeometryMap<T> >
{
    typedef T Scalar;
    typedef gsGeometryMap<T> const & Nested_t;
};

/*
   An element object collecting relevant expressions
*/
template<class T>
class gsFeElement
{
    friend class cdiam_expr<T>;

    const gsDomainIterator<T> * m_di; ///< Pointer to the domain iterator

    cdiam_expr<T> cd;
public:
    typedef T Scalar;

    gsFeElement() : m_di(NULL), cd(*this) { }

    void set(const gsDomainIterator<T> & di)
    { m_di = &di; }

    /// The diameter of the element
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
    const gsFeElement<T> & _e; ///<Reference to the element
public:
    typedef T Scalar;

    enum {ScalarValued = 1};

    explicit cdiam_expr(const gsFeElement<T> & el) : _e(el) { }

    T eval(const index_t ) const { return _e.m_di->getCellSize(); }

    inline cdiam_expr<T> val() const { return *this; }
    inline index_t rows() const { return 0; }
    inline index_t cols() const { return 0; }
    inline void setFlag() const { }
    void parse(gsSortedVector<const gsFunctionSet<Scalar>*> & ) const { }

    const gsFeSpace<T> & rowVar() const { return gsNullExpr<T>(); }
    const gsFeSpace<T> & colVar() const { return gsNullExpr<T>(); }

    enum{rowSpan = 0, colSpan = 0};

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
    void parse(gsSortedVector<const gsFunctionSet<Scalar>*> &) const { }
    const gsFeVariable<T> & rowVar() const { gsNullExpr<T>(); }
    const gsFeVariable<T> & colVar() const { gsNullExpr<T>(); }
    enum{rowSpan = 0, colSpan = 0};

    void print(std::ostream &os) const
    { os << "diam(e)"; }
};
*/

template <typename T>
struct expr_traits<gsFeVariable<T> >
{
    typedef T Scalar;
    typedef const gsFeVariable<T> & Nested_t;
};

/**
   Expression for finite element variables or PDE coefficient functionals

   This can be e.g. a diffusion coefficient, or a function that
   belongs to an isogeometric function space (derived class gsFeSpace)
*/
template<class T>
class gsFeVariable  : public symbol_expr< gsFeVariable<T> >
{
    friend class gismo::gsExprHelper<T>;
    typedef symbol_expr< gsFeVariable<T> > Base;
protected:
    explicit gsFeVariable(index_t _d = 1) : Base(_d) { }
public:
        enum{rowSpan = 0};
};

/**
   Expression for finite element variable in an isogeometric function
   space

   This corresponds to an FE variable that
   belongs to an isogeometric function space
*/
template<class T>
class gsFeSpace :public symbol_expr< gsFeSpace<T> >
//class gsFeSpace : public gsFeVariable<T>
{
    friend class gsNullExpr<T>;
    friend class gsFeSolution<T>;

protected:
    typedef symbol_expr< gsFeSpace<T> > Base;
    //typedef gsFeVariable<T> Base;

    index_t             m_id;

    // C^r coupling
    mutable index_t m_r; // iFace::handling

    typedef gsBoundaryConditions<T> bcList;
    typedef typename bcList::bcRefList bcRefList;
    mutable bcRefList m_bcs;
    mutable gsDofMapper m_mapper;
    gsMatrix<T> m_fixedDofs;

public:
        enum{rowSpan = 1};

public:

    typedef const gsFeSpace & Nested_t;

    enum {Space = 1};

    typedef T Scalar;

    const gsFeSpace<T> & rowVar() const {return *this;}

    gsDofMapper & mapper() {return m_mapper;}
    const gsDofMapper & mapper() const {return m_mapper;}

    inline const gsMatrix<T> & fixedPart() const {return m_fixedDofs;}
    gsMatrix<T> & fixedPart() {return m_fixedDofs;}

    const bcRefList & bc() const { return m_bcs; }
    void addBc(bcRefList bc) const { m_bcs = bc; }
    void clearBc() const { m_bcs.clear(); }
    size_t bcSize() const { return m_bcs.size(); }

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
        const index_t sz  = mb[p].size();
        result.resize(sz, dim); // (!)

        for (index_t c = 0; c!=dim; c++) // for all components
        {
            for (index_t i = 0; i < sz; ++i) // loop over all basis functions (even the eliminated ones)
            {
                const int ii = m_mapper.index(i, p, c);
                if ( m_mapper.is_free_index(ii) ) // DoF value is in the solVector
                        result(i,c) = solVector.at(ii);
                else // eliminated DoF: fill with Dirichlet data
                {
                    result(i,c) =  m_fixedDofs.at( m_mapper.global_to_bindex(ii) );
                }
            }
        }
    }

    // space restrictTo(boundaries);
    // space restrictTo(bcRefList domain);

    void reset()
    {
        m_mapper = gsDofMapper(); //reset ?

        if (const gsMultiBasis<T> * mb =
            dynamic_cast<const gsMultiBasis<T>*>(&this->source()) )
        {
            m_mapper = gsDofMapper(*mb, this->dim());
            //m_mapper.init(*mb); //bug
            if ( 0==this->interfaceCont() ) // Conforming boundaries ?
            {
                for ( gsBoxTopology::const_iiterator it = mb->topology().iBegin();
                      it != mb->topology().iEnd(); ++it )
                {
                    mb->matchInterface(*it, m_mapper);
                }
            }

            gsMatrix<index_t> bnd;
            for (typename bcRefList::const_iterator
                     it = this->bc().begin() ; it != this->bc().end(); ++it )
            {
                GISMO_ASSERT(static_cast<size_t>(it->get().ps.patch) < this->mapper().numPatches(),
                              "Problem: a boundary condition is set on a patch id which does not exist.");

                bnd = mb->basis(it->get().ps.patch).boundary( it->get().ps.side() );
                for (index_t c=0; c!= this->dim(); c++) // for all components
                    m_mapper.markBoundary(it->get().ps.patch, bnd, c);
            }
        }
        else if (const gsBasis<T> * b =
                 dynamic_cast<const gsBasis<T>*>(&this->source()) )
        {
            m_mapper = gsDofMapper(*b);
            gsMatrix<index_t> bnd;
            for (typename bcRefList::const_iterator
                     it = this->bc().begin() ; it != this->bc().end(); ++it )
            {
                GISMO_ASSERT( it->get().ps.patch == 0,
                              "Problem: a boundary condition is set on a patch id which does not exist.");

                bnd = b->boundary( it->get().ps.side() );
                m_mapper.markBoundary(0, bnd);
            }
        }
        else
        {
            GISMO_ASSERT( 0 == this->bc().size(), "Problem: BCs are ignored.");
            m_mapper.setIdentity(this->source().nPieces(), this->source().size());
        }

        m_mapper.finalize();
        //this->mapper().print();
    }

    void setup(const gsBoundaryConditions<T> & bc, const index_t dir_values, const index_t _icont = -1) const
    {
        this->setInterfaceCont(_icont);

        m_mapper = gsDofMapper(); //reset ?

        const gsMultiBasis<T> * smb = static_cast<const gsMultiBasis<T>*>(&this->source());

        if (const gsMultiBasis<T> * mb =
            dynamic_cast<const gsMultiBasis<T>*>(&this->source()) )
        {
            m_mapper = gsDofMapper(*mb, this->dim() );
            //m_mapper.init(*mb, this->dim()); //bug
            if ( 0==this->interfaceCont() ) // Conforming boundaries ?
            {
                for ( gsBoxTopology::const_iiterator it = mb->topology().iBegin();
                      it != mb->topology().iEnd(); ++it )
                {
                    mb->matchInterface(*it, m_mapper);
                }
            }

            // DIRICHLET
            gsMatrix<index_t> bnd;
            for (typename bcList::const_iterator
                     it = bc.begin("Dirichlet") ; it != bc.end("Dirichlet"); ++it )
            {
                const index_t cc = it->unkComponent();
                //if (it->unknown() == -1 || it->unkown() == this->id())

                GISMO_ASSERT(static_cast<size_t>(it->ps.patch) < this->mapper().numPatches(),
                              "Problem: a boundary condition is set on a patch id which does not exist.");

                bnd = mb->basis(it->ps.patch).boundary( it->ps.side() );
                if (cc==-1)
                    for (index_t c=0; c!= this->dim(); c++) // for all components
                        m_mapper.markBoundary(it->ps.patch, bnd, c);
                else
                    m_mapper.markBoundary(it->ps.patch, bnd, cc);
            }
            // CLAMPED
            gsMatrix<index_t> bnd1;
            for (typename bcList::const_iterator
                     it = bc.begin("Clamped") ; it != bc.end("Clamped"); ++it )
            {
                const index_t cc = it->unkComponent();
                //if (it->unknown() == -1 || it->unkown() == this->id())

                GISMO_ASSERT(static_cast<size_t>(it->ps.patch) < this->mapper().numPatches(),
                              "Problem: a boundary condition is set on a patch id which does not exist.");

                bnd  = mb->basis(it->ps.patch).boundaryOffset( it->ps.side(), 0);
                bnd1 = mb->basis(it->ps.patch).boundaryOffset( it->ps.side(), 1);
                // Cast to tensor b-spline basis
                if ( mb != NULL) // clamp adjacent dofs
                {
                    if ( ! it->ps.parameter() )
                        bnd.swap(bnd1);
                    if (cc==-1)
                        for (index_t c=0; c!= this->dim(); c++) // for all components
                            for ( index_t k=0; k<bnd.size(); ++k)
                                m_mapper.matchDof( it->ps.patch, (bnd)(k,0),
                                                   it->ps.patch, (bnd1)(k,0) , c);
                    else
                        for ( index_t k=0; k<bnd.size(); ++k)
                            m_mapper.matchDof( it->ps.patch, (bnd)(k,0),
                                           it->ps.patch, (bnd1)(k,0) , cc);
                }
                else
                    gsWarn<<"Unable to apply clamped condition.\n";
            }

            // COLLAPSED
            for (typename bcList::const_iterator
                     it = bc.begin("Collapsed") ; it != bc.end("Collapsed"); ++it )
            {
                const index_t cc = it->unkComponent();
                //if (it->unknown() == -1 || it->unkown() == this->id())

                GISMO_ASSERT(static_cast<size_t>(it->ps.patch) < this->mapper().numPatches(),
                              "Problem: a boundary condition is set on a patch id which does not exist.");

                bnd = mb->basis(it->ps.patch).boundary( it->ps.side() );

                // Cast to tensor b-spline basis
                if ( mb != NULL) // clamp adjacent dofs
                {
                    // match all DoFs to the first one of the side
                    if (cc==-1)
                        for (index_t c=0; c!= this->dim(); c++) // for all components
                            for ( index_t k=0; k<bnd.size()-1; ++k)
                                m_mapper.matchDof( it->ps.patch, (bnd)(0,0),
                                                   it->ps.patch, (bnd)(k+1,0) , c);
                    else
                        for ( index_t k=0; k<bnd.size()-1; ++k)
                                m_mapper.matchDof( it->ps.patch, (bnd)(0,0),
                                                   it->ps.patch, (bnd)(k+1,0) , cc);
                }
            }

            // corners
            for (typename gsBoundaryConditions<T>::const_citerator
                     it = bc.cornerBegin() ; it != bc.cornerEnd(); ++it )
            {
                //assumes (unk == -1 || it->unknown == unk)
                GISMO_ASSERT(static_cast<size_t>(it->patch) < smb->nBases(),
                             "Problem: a corner boundary condition is set on a patch id which does not exist.");
                m_mapper.eliminateDof(smb->basis(it->patch).functionAtCorner(it->corner), it->patch);
            }

        }
        else if (const gsBasis<T> * b =
                 dynamic_cast<const gsBasis<T>*>(&this->source()) )
        {
            m_mapper = gsDofMapper(*b);
            gsMatrix<index_t> bnd;
            for (typename bcList::const_iterator
                     it = bc.begin("Dirichlet") ; it != bc.end("Dirichlet"); ++it )
            {
                GISMO_ASSERT( it->ps.patch == 0,
                              "Problem: a boundary condition is set on a patch id which does not exist.");

                bnd = b->boundary( it->ps.side() );
                //const index_t cc = it->unkComponent();
                m_mapper.markBoundary(0, bnd, 0);
            }

            m_mapper = gsDofMapper(*b);
            for (typename bcList::const_iterator
                     it = bc.begin("Clamped") ; it != bc.end("Clamped"); ++it )
            {
                GISMO_ASSERT( it->ps.patch == 0,
                              "Problem: a boundary condition is set on a patch id which does not exist.");

                bnd = b->boundary( it->ps.side() );
                //const index_t cc = it->unkComponent();
                // m_mapper.markBoundary(0, bnd, 0);
            }

            m_mapper = gsDofMapper(*b);
            for (typename bcList::const_iterator
                     it = bc.begin("Collapsed") ; it != bc.end("Collapsed"); ++it )
            {
                GISMO_ASSERT( it->ps.patch == 0,
                              "Problem: a boundary condition is set on a patch id which does not exist.");

                bnd = b->boundary( it->ps.side() );
                //const index_t cc = it->unkComponent();
                // m_mapper.markBoundary(0, bnd, 0);
            }
        }
        else
        {
            GISMO_ASSERT( 0 == bc.size(), "Problem: BCs are ignored.");
            m_mapper.setIdentity(this->source().nPieces(), this->source().size());
        }

        m_mapper.finalize();

        // No more BCs
        //m_bcs = bc.get("Dirichlet");
        gsDirichletValues(bc, dir_values, *this);

        // corner values (overrides edge BCs)
        gsMatrix<T> & fixedDofs = const_cast<gsMatrix<T> &>(m_fixedDofs);
        for ( typename gsBoundaryConditions<T>::const_citerator
                  it = bc.cornerBegin(); it != bc.cornerEnd(); ++it )
        {
            const int i  = smb->basis(it->patch).functionAtCorner(it->corner);
            const int ii = m_mapper.bindex( i , it->patch, 0 );//component=0 for now! Todo.
            fixedDofs.at(ii) = it->value;
        }
    }

protected:
    friend class gismo::gsExprHelper<T>;
    explicit gsFeSpace(index_t _d = 1) : Base(_d), m_id(-1), m_r(-1)
    { }
};


/*
   Expression representing a function given by a vector of
   coefficients in a gsFeSpace.

   Typically it used for accessing the solution of a boundary-value
   problem.
*/
template<class T>
class gsFeSolution :public _expr<gsFeSolution<T> >
{
protected:
    const gsFeSpace<T> & _u;

    gsMatrix<T> * _Sv; ///< Pointer to a coefficient vector
public:
    typedef T Scalar;

    explicit gsFeSolution(const gsFeSpace<T> & u)
    : //Base(u.dim())
      //Base(u)
      _u(u), _Sv(NULL)
    { }

    gsFeSolution(const gsFeSpace<T> & u, gsMatrix<T> & Sv) : _u(u), _Sv(&Sv)
    { }

    mutable gsMatrix<T> res;
    const gsMatrix<T> & eval(index_t k) const
    {
        GISMO_ASSERT(1==_u.data().actives.cols(), "Single actives expected");

        res.setZero(_u.dim(), 1);
        const gsDofMapper & map = _u.mapper();

        for (index_t c = 0; c!=_u.dim(); c++) // for all components
        {
            for (index_t i = 0; i!=_u.data().actives.size(); ++i)
            {
                const index_t ii = map.index(_u.data().actives.at(i), _u.data().patchId, c);
                if ( map.is_free_index(ii) ) // DoF value is in the solVector
                {
                    // gsDebugVar(_Sv->at(ii));
                    // gsDebugVar(_u.data().values[0](i,k));
                    res.at(c) += _Sv->at(ii) * _u.data().values[0](i,k);
                }
                else
                    res.at(c) += _u.data().values[0](i,k) *
                        _u.fixedPart().at( map.global_to_bindex(ii) );
            }
        }
        return res;
    }

    enum{rowSpan = 0, colSpan = 0};

    index_t rows() const {return _u.dim(); }
    static index_t cols() {return 1; }

    void setFlag() const
    { _u.data().flags |= NEED_VALUE | NEED_ACTIVE; }

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

    void setSolutionVector(gsMatrix<T>& solVector)
    { _Sv = & solVector; }

    const gsMatrix<T> & coefs() const { return *_Sv; }

    //gsMatrix<T> & coefs() { return *_Sv; } // wd4702 ?

    /// Extract the coefficients of piece piece \a p
    void extract(gsMatrix<T> & result, const index_t p = 0) const
    {
        // GISMO_ASSERT(NULL!=_Sv, "Solution vector not set, call setSolutionVector"); // wd4702 ?
        _u.getCoeffs(*_Sv, result, p);
    }

    /// Extracts ALL the coefficients in a solution vector; including coupled and boundary DoFs
    void extractFull(gsMatrix<T> & result) const
    {
        index_t offset, ii, bi;
        result.resize(_u.mapper().mapSize(),1);
        for (index_t c=0; c!=_u.m_d; c++)
            for (size_t j=0; j!=_u.mapper().numPatches(); j++)
                for (size_t i=0; i!=_u.mapper().patchSize(j); i++) // loop over all DoFs (free and eliminated)
                {
                    offset = _u.mapper().offset(j);

                    ii = _u.mapper().index(i,j,c); // global index
                    if (_u.mapper().is_boundary(i,j,c))
                    {
                        bi = _u.mapper().global_to_bindex(ii); // boundary index
                        result(i+offset,0) = _u.fixedPart().at(bi);
                    }
                    else
                        result(i+offset,0) = _Sv->at(ii);
                }
    }

    /// Extract this variable as a multipatch object
    void extract(gsMultiPatch<T> & result) const
    {
        result.clear();

        if( const gsMultiBasis<T>* basis = dynamic_cast<const gsMultiBasis<T>* >(&_u.source()) )
            for (size_t i = 0; i != basis->nBases(); ++i)
            {
                memory::unique_ptr<gsGeometry<T> > p(this->extractPiece(i));
                result.addPatch(*p);
            }
    }

    /// Extract the piece \a p as a gsGeometry pointer
    memory::unique_ptr<gsGeometry<T> > extractPiece(const index_t p) const
    {
        if ( const gsBasis<T> * b = dynamic_cast<const gsBasis<T>*>(&_u.source().piece(p)) )
        {
            gsMatrix<T> cf;
            extract(cf, p);
            return b->makeGeometry(give(cf));
        }
        GISMO_ERROR("gsFeSolution: Extraction error");
    }

    // insert g-coefficients to the solution vector
    void insert(const gsGeometry<T> & g, const index_t p = 0) const
    {
        const index_t dim = _u.dim();
        const gsMatrix<T> & cf = g.coefs();
        gsMatrix<T> & sol = *_Sv;
        //gsMatrix<T> & fixedPart = _u.fixedPart();
        const gsDofMapper & mapper = _u.mapper();
        for (index_t c = 0; c!=_u.dim(); c++) // for all components
        {
            for (index_t i = 0; i != cf.rows(); ++i)
            {
                const index_t ii = mapper.index(i, p, c);
                if ( mapper.is_free_index(ii) ) // DoF value is in the solVector
                    sol.at(ii) = cf(i, c);
                /*
                  else
                  {
                  fixedPart.row(m_mapper.global_to_bindex(ii)) = cf.row(i);
                  }
                */
            }
        }
    }
};

/*
   Expression for the gradient of a gsFeSolution
 */
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
        for (index_t c = 0; c!= _u.dim(); c++)
        {
            for (index_t i = 0; i!=_u.data().actives.size(); ++i)
            {
                const index_t ii = map.index(_u.data().actives.at(i), _u.data().patchId,c);
                if ( map.is_free_index(ii) ) // DoF value is in the solVector
                {
                        res.row(c) += _u.coefs().at(ii) *
                            _u.data().values[1]
                            //.block(i*_u.parDim(),k,_u.parDim(),1).transpose();
                            .col(k).segment(i*_u.parDim(), _u.parDim()).transpose();
                }
                else
                {
                    res.noalias() +=
                        _u.fixedPart().row( map.global_to_bindex(ii) ).asDiagonal() *
                        _u.data().values[1].col(k).segment(i*_u.parDim(), _u.parDim())
                        .transpose().replicate(_u.dim(),1);
                }
            }
        }
        return res;
    }

    enum{rowSpan = 0, colSpan = 0};

    index_t rows() const {return _u.dim();}

    index_t cols() const {return _u.parDim(); }

    void setFlag() const
    {
        _u.data().flags |= NEED_GRAD|NEED_ACTIVE;
    }

    void parse(gsSortedVector<const gsFunctionSet<Scalar>*> & evList) const
    {
        //GISMO_ASSERT(NULL!=m_fd, "FeVariable: FuncData member not registered");
        evList.push_sorted_unique(&_u.source());
        _u.data().flags |= NEED_GRAD|NEED_ACTIVE;
    }

    void print(std::ostream &os) const { os << "grad(s)"; }
};


//template<class T>
//class gsFeMutVariable  : public _expr<gsFeMutVariable<T> >
// member gsFuncData

/*
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
    enum {Space = (E::Space==0?0:(E::Space==1?2:1))};

    mutable gsMatrix<Scalar> res;

    //MatExprType eval(const index_t k) const
    const gsMatrix<Scalar> & eval(const index_t k) const
    {
        //return _u.eval(k).transpose();
        // /*
        if (E::ColBlocks)
        {
            res = _u.eval(k).blockTranspose( _u.cardinality() ); return res;
            //return _u.eval(k).blockTranspose(_u.cols()/_u.rows());
        }
        else
        {
            res = _u.eval(k).transpose(); return res;
            //return _u.eval(k).blockTranspose(1); // buggy ?
        }
        //*/
    }

    index_t rows() const
    {
        // if (ColBlocks) // note: sq. blocks assumed
        //     return _u.rows();
        // else
            return _u.cols();
    }

    index_t cols() const
    {
        // if (ColBlocks) // note: sq. blocks assumed
        //     return _u.cols();
        // else
            return _u.rows();
    }
    void setFlag() const { _u.setFlag(); }

    void parse(gsSortedVector<const gsFunctionSet<Scalar>*> & evList) const
    { _u.parse(evList); }

    const gsFeSpace<Scalar> & rowVar() const { return _u.colVar(); }
    const gsFeSpace<Scalar> & colVar() const { return _u.rowVar(); }

    index_t cardinality_impl() const { return _u.cardinality_impl(); }

    enum{rowSpan = E::colSpan, colSpan = E::rowSpan};

    void print(std::ostream &os) const { os<<"("; _u.print(os); os <<")'"; }
private:
/*
    template<class U> EIGEN_STRONG_INLINE MatExprType
    eval_impl(const U k, typename util::enable_if<1==ColBlocks,U>::type* = nullptr)
    { return _u.eval(k).blockTranspose(_u.cols()/_u.rows()); }

    template<class U> EIGEN_STRONG_INLINE MatExprType
    eval_impl(const U k, typename util::enable_if<0==ColBlocks,U>::type* = nullptr)
    { return _u.eval(k).transpose(); }
*/
};


/*
   Expression for an evaluation of the (sub-)expression in temporary memory
*/
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

    const gsFeSpace<Scalar> & rowVar() const { return _u.rowVar(); }
    const gsFeSpace<Scalar> & colVar() const { return _u.colVar(); }

    enum{rowSpan = E::rowSpan, colSpan = E::colSpan};

    void print(std::ostream &os) const { _u.print(os); }
};

/*
   Expression for the trace of a (matrix) expression
 */
template<class E>
class trace_expr  : public _expr<trace_expr<E> >
{
public:
    typedef typename E::Scalar Scalar;
    enum {ScalarValued = 0};
    enum {Space = E::Space};
private:
    typename E::Nested_t _u;
    mutable gsMatrix<Scalar> res;

public:
    trace_expr(_expr<E> const& u) : _u(u)
    {
        // gcc 4.8.4: invalid read due to _u.rows() using gsFuncData
        //GISMO_ASSERT(0== _u.cols()%_u.rows(), "Expecting square-block expression, got " << _u.rows() <<" x "<< _u.cols() );
    }

    // choose if ColBlocks
    const gsMatrix<Scalar> & eval(const index_t k) const
    {
        // Assume mat ??
        MatExprType tmp = _u.eval(k);
        const index_t cb = _u.rows();
        const index_t r  = _u.cols() / cb;
        res.resize(r, 1);
        for (index_t i = 0; i!=r; ++i){
            res(i,0) = tmp.middleCols(i*cb,cb).trace();
        }
        return res;
    }

    // choose if !ColBlocks
    //todo: Scalar eval(const index_t k) const

    index_t rows() const { return _u.cols() / _u.rows(); }
    index_t cols() const { return 1; }
    void setFlag() const { _u.setFlag(); }

    void parse(gsSortedVector<const gsFunctionSet<Scalar>*> & evList) const
    { _u.parse(evList); }

    const gsFeSpace<Scalar> & rowVar() const { return _u.rowVar(); }
    const gsFeSpace<Scalar> & colVar() const { return _u.colVar(); }

    enum{rowSpan = E::rowSpan, colSpan = E::colSpan};

    void print(std::ostream &os) const { os << "trace("; _u.print(os); os<<")"; }
};

/*
   Expression for the adjugate of a (matrix) expression
 */
template<class E>
class adjugate_expr  : public _expr<adjugate_expr<E> >
{
public:
    typedef typename E::Scalar Scalar;
    enum {ScalarValued = 0, ColBlocks = E::ColBlocks};
    enum {Space = E::Space};
private:
    typename E::Nested_t _u;
    mutable gsMatrix<Scalar> res;

public:
    adjugate_expr(_expr<E> const& u) : _u(u)
    {
        // gcc 4.8.4: invalid read due to _u.rows() using gsFuncData
        //GISMO_ASSERT(0== _u.cols()%_u.rows(), "Expecting square-block expression, got " << _u.rows() <<" x "<< _u.cols() );
    }

    // choose if ColBlocks
    const gsMatrix<Scalar> & eval(const index_t k) const
    {
        // Assume mat ??
        MatExprType tmp = _u.eval(k);
        const index_t cb = _u.rows();
        const index_t r  = _u.cols() / cb;
        res.resize(_u.rows(),_u.cols());
        for (index_t i = 0; i!=r; ++i){
            res.middleCols(i*cb,cb) = tmp.middleCols(i*cb,cb).adjugate();
        }
        return res;
    }

    // choose if !ColBlocks
    //todo: Scalar eval(const index_t k) const

    index_t rows() const { return _u.rows(); }
    index_t cols() const { return _u.cols(); }
    void setFlag() const { _u.setFlag(); }

    void parse(gsSortedVector<const gsFunctionSet<Scalar>*> & evList) const
    { _u.parse(evList); }

    const gsFeSpace<Scalar> & rowVar() const { return _u.rowVar(); }
    const gsFeSpace<Scalar> & colVar() const { return _u.colVar(); }

    enum{rowSpan = E::rowSpan, colSpan = E::colSpan};

    void print(std::ostream &os) const { os << "adj("; _u.print(os); os<<")"; }
};

template<class E>
class reshape_expr  : public _expr<reshape_expr<E> >
{
public:
    typedef typename E::Scalar Scalar;
    enum {ScalarValued = 0};
private:
    typename E::Nested_t _u;
    index_t _n, _m;
    mutable gsMatrix<Scalar> tmp;

public:

    //the reshaping is done column-wise
    reshape_expr(_expr<E> const& u, index_t n, index_t m) : _u(u), _n(n), _m(m)
    {
        //GISMO_ASSERT( _u.rows()*_u.cols() == _n*_m, "Wrong dimension"); //
    }

    const gsAsConstMatrix<Scalar> eval(const index_t k) const
    {
        // Note: this assertion would fail in the constructore!
        GISMO_ASSERT( _u.rows()*_u.cols() == _n*_m, "Wrong dimension "
         << _u.rows() << " x "<<_u.cols() << "!=" << _n << " * "<< _m );
        tmp = _u.eval(k);
        return gsAsConstMatrix<Scalar>(tmp.data(),_n,_m);
    }

    index_t rows() const { return _n; }
    index_t cols() const { return _m; }
    void setFlag() const { _u.setFlag(); }

    void parse(gsSortedVector<const gsFunctionSet<Scalar>*> & evList) const
    { _u.parse(evList); }

    const gsFeSpace<Scalar> & rowVar() const { return _u.rowVar(); }
    const gsFeSpace<Scalar> & colVar() const { return _u.colVar(); }

    enum{rowSpan = E::rowSpan, colSpan = E::colSpan};

    void print(std::ostream &os) const { os << "reshape("; _u.print(os); os<<","<<_n<<","<<_m<<")"; }
};

/// Reshape an expression
template <typename E> EIGEN_STRONG_INLINE
reshape_expr<E> const reshape(E const & u, index_t n, index_t m)
{ return reshape_expr<E>(u, n, m); }

template<class E>
class replicate_expr  : public _expr<replicate_expr<E> >
{
public:
    typedef typename E::Scalar Scalar;
    enum {ScalarValued = 0};
private:
    typename E::Nested_t _u;
    index_t _n, _m;
    mutable gsMatrix<Scalar> tmp;

public:

    //the replicate is done nxm times
    replicate_expr(_expr<E> const& u, index_t n, index_t m) : _u(u), _n(n), _m(m)
    {
    }

    MatExprType eval(const index_t k) const
    {
        tmp = _u.eval(k);
        return tmp.replicate(_n,_m);
    }

    index_t rows() const { return _n*_u.rows(); }
    index_t cols() const { return _m*_u.cols(); }
    void setFlag() const { _u.setFlag(); }

    void parse(gsSortedVector<const gsFunctionSet<Scalar>*> & evList) const
    { _u.parse(evList); }

    const gsFeSpace<Scalar> & rowVar() const { return _u.rowVar(); }
    const gsFeSpace<Scalar> & colVar() const { return _u.colVar(); }
    index_t cardinality_impl() const { return _u.cardinality_impl(); }

    enum{rowSpan = E::rowSpan, colSpan = E::colSpan};

    void print(std::ostream &os) const { os << "replicate("; _u.print(os); os<<","<<_n<<","<<_m<<")"; }
};

/// Replicate an expression
template <typename E> EIGEN_STRONG_INLINE
replicate_expr<E> const replicate(E const & u, index_t n, index_t m = 1)
{ return replicate_expr<E>(u, n, m); }

/**
  Transforms a matrix expression into a vector expression by computing the vector
  [ a b c+d]^T
  for each matrix block
  [ a d ]
  [ c b ]
*/
template<class E>
class flat_expr  : public _expr<flat_expr<E> >
{
public:
    typedef typename E::Scalar Scalar;
    enum {ScalarValued = 0, Space = E::Space};
private:
    typename E::Nested_t _u;
    mutable gsMatrix<Scalar> tmp;

public:

    flat_expr(_expr<E> const& u) : _u(u)
    {
        //GISMO_ASSERT( _u.rows()*_u.cols() == _n*_m, "Wrong dimension"); //
    }

    const gsMatrix<Scalar> & eval(const index_t k) const
    {
        //_u.print(gsInfo); gsInfo<<std::endl;
        //GISMO_ASSERT( _u.rows()==2 , "Wrong dimension, got " << _u.rows() << ", " << _u.cols());

        tmp = _u.eval(k);
        const index_t numActives = _u.cardinality();

        for (index_t i = 0; i<numActives; ++i)
        {
            tmp(0,2*i+1) += tmp(1,2*i);
            // GISMO_ASSERT(tmp(0,2*i+1)==tmp(1,2*i),"ALERT!, tmp(0,2*i+1) = "<<tmp(0,2*i+1)<<", tmp(1,2*i) = "<<tmp(1,2*i));
            std::swap(tmp(1,2*i), tmp(1,2*i+1));
        }

        tmp.resize(4,numActives);
        tmp.conservativeResize(3,numActives);

        //tmp.resize(numActives,4);
        //tmp.conservativeResize(numActives,3);


        // gsDebugVar(rowSpan);
        // gsDebugVar(colSpan);
        if ( rowSpan )
            tmp.transposeInPlace();
        else if (!colSpan) // if not colSpan and not rowSpan
            tmp.transposeInPlace();

        return tmp;
    }

    index_t rows() const { return 1; }
    index_t cols() const { return 3; }
    void setFlag() const { _u.setFlag(); }

    void parse(gsSortedVector<const gsFunctionSet<Scalar>*> & evList) const
    { _u.parse(evList); }

    const gsFeSpace<Scalar> & rowVar() const { return _u.rowVar(); }
    const gsFeSpace<Scalar> & colVar() const { return _u.colVar(); }
    index_t cardinality_impl() const { return _u.cardinality_impl(); }

    enum{rowSpan = E::rowSpan, colSpan = E::colSpan};

    void print(std::ostream &os) const { os << "flat("; _u.print(os); os<<")"; }
};

/// Make a matrix 2x2 expression "flat"
template <typename E> EIGEN_STRONG_INLINE
flat_expr<E> const flat(E const & u)
{ return flat_expr<E>(u); }


/// Takes the row-sum of a matrix expression
template<class E>
class rowsum_expr  : public _expr<rowsum_expr<E> >
{
public:
    typedef typename E::Scalar Scalar;
    enum {ScalarValued = 0, Space = E::Space};
private:
    typename E::Nested_t _u;
    mutable gsMatrix<Scalar> tmp;

public:

    rowsum_expr(_expr<E> const& u) : _u(u)
    {
        //GISMO_ASSERT( _u.rows()*_u.cols() == _n*_m, "Wrong dimension"); //
    }

    const gsMatrix<Scalar> & eval(const index_t k) const
    {
        tmp = _u.eval(k).rowwise().sum();
        return tmp;
    }

    index_t rows() const { return _u.rows(); }
    index_t cols() const { return 1; }
    void setFlag() const { _u.setFlag(); }

    void parse(gsSortedVector<const gsFunctionSet<Scalar>*> & evList) const
    { _u.parse(evList); }

    const gsFeSpace<Scalar> & rowVar() const { return _u.rowVar(); }
    const gsFeSpace<Scalar> & colVar() const { return _u.colVar(); }
    index_t cardinality_impl() const { return _u.cardinality_impl(); }

    enum{rowSpan = E::rowSpan, colSpan = E::colSpan};

    void print(std::ostream &os) const { os << "rowsum("; _u.print(os); os<<")"; }
};

/// Takes the column-sum of a matrix or vector expression
template<class E>
class colsum_expr  : public _expr<colsum_expr<E> >
{
public:
    typedef typename E::Scalar Scalar;
    enum {ScalarValued = 0, Space = E::Space};
private:
    typename E::Nested_t _u;
    mutable gsMatrix<Scalar> tmp;

public:

    colsum_expr(_expr<E> const& u) : _u(u)
    {
        //GISMO_ASSERT( _u.rows()*_u.cols() == _n*_m, "Wrong dimension"); //
    }

    const gsMatrix<Scalar> & eval(const index_t k) const
    {
        tmp = _u.eval(k).colwise().sum();
        return tmp;
    }

    index_t rows() const { return _u.rows(); }
    index_t cols() const { return 1; }
    void setFlag() const { _u.setFlag(); }

    void parse(gsSortedVector<const gsFunctionSet<Scalar>*> & evList) const
    { _u.parse(evList); }

    const gsFeSpace<Scalar> & rowVar() const { return _u.rowVar(); }
    const gsFeSpace<Scalar> & colVar() const { return _u.colVar(); }
    index_t cardinality_impl() const { return _u.cardinality_impl(); }

    enum{rowSpan = E::rowSpan, colSpan = E::colSpan};

    void print(std::ostream &os) const { os << "colsum("; _u.print(os); os<<")"; }
};


/*
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

    const gsFeSpace<Scalar> & rowVar() const { return _u.rowVar(); }
    const gsFeSpace<Scalar> & colVar() const { return _u.colVar(); }

    enum{rowSpan = E::rowSpan, colSpan = E::colSpan};

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
    const gsFeSpace<Scalar> & rowVar() const {return gsNullExpr<Scalar>::get();} \
    const gsFeSpace<Scalar> & colVar() const {return gsNullExpr<Scalar>::get();} \
    void print(std::ostream &os) const                                  \
    { os << #name <<"("; _u.print(os); os <<")"; }                      \
    enum{rowSpan = E::rowSpan, colSpan = E::colSpan};  };               \
    // const gsFeSpace<Scalar> & rowVar() const { return _u.rowVar(); }
    // const gsFeSpace<Scalar> & colVar() const { return _u.colVar(); }

/// Eucledian Norm
GISMO_EXPR_VECTOR_EXPRESSION(norm,norm,1);
/// Squared Eucledian Norm
GISMO_EXPR_VECTOR_EXPRESSION(sqNorm,squaredNorm,1);
/// Normalization of a vector to unit measure
GISMO_EXPR_VECTOR_EXPRESSION(normalized,normalized,0); // (!) mem.
/// Inverse of a matrix expression
GISMO_EXPR_VECTOR_EXPRESSION(inv,inverse,0);
// GISMO_EXPR_VECTOR_EXPRESSION(cwSqr,array().square,0)
// GISMO_EXPR_VECTOR_EXPRESSION(sum,array().sum,1)
// GISMO_EXPR_VECTOR_EXPRESSION(sqrt,array().sqrt,0)
// GISMO_EXPR_VECTOR_EXPRESSION(abs,array().abs,0)

//Determinant
GISMO_EXPR_VECTOR_EXPRESSION(det,determinant,1);


//GISMO_EXPR_VECTOR_EXPRESSION(replicate,replicate,0);

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

    const gsFeSpace<Scalar> & rowVar() const { return _u.rowVar(); }
    const gsFeSpace<Scalar> & colVar() const { return _u.colVar(); }

    enum{rowSpan = E::rowSpan, colSpan = E::colSpan};

    index_t rows() const { return _u.rows(); }
    index_t cols() const { return _u.rows() * _u.cols(); }
    void setFlag() const { _u.setFlag(); }
    void parse(gsSortedVector<const gsFunctionSet<Scalar>*> & evList) const
    { _u.parse(evList); }

    void print(std::ostream &os) const { os << "diag("; _u.print(os); os <<")";}
};

/**
   Expression for the identity matrix
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

    gsMatrix<Scalar>::IdentityReturnType eval(const index_t) const
    {
        return gsMatrix<Scalar>::Identity(_dim,_dim);
    }

    index_t rows() const { return _dim; }
    index_t cols() const { return  _dim; }
    void setFlag() const { }
    void parse(gsSortedVector<const gsFunctionSet<Scalar>*> & ) const {  }

    enum{rowSpan = 0, colSpan = 0};

    const gsFeSpace<Scalar> & rowVar() const {return gsNullExpr<Scalar>::get();}
    const gsFeSpace<Scalar> & colVar() const {return gsNullExpr<Scalar>::get();}

    void print(std::ostream &os) const { os << "id("<<_dim <<")";}
};

/**
   Expression for the sign of another expression
 */
template<class E>
class sign_expr : public _expr<sign_expr<E> >
{
    typename E::Nested_t _u;
public:
    typedef typename E::Scalar Scalar;
    enum {ScalarValued = 1 };
    enum {Space = E::Space};

    sign_expr(_expr<E> const& u) : _u(u) { }

    Scalar eval(const index_t k) const
    {
        const Scalar v = _u.val().eval(k);
        return ( v>0 ? 1 : ( v<0 ? -1 : 0 ) );
    }

    static index_t rows() { return 0; }
    static index_t cols() { return 0; }

    void setFlag() const { _u.setFlag(); }

    void parse(gsSortedVector<const gsFunctionSet<Scalar>*> & ) const {  }

    static bool isScalar() { return true; }

    enum{rowSpan = 0, colSpan = 0};

    const gsFeSpace<Scalar> & rowVar() const {return gsNullExpr<Scalar>::get();}
    const gsFeSpace<Scalar> & colVar() const {return gsNullExpr<Scalar>::get();}

    void print(std::ostream &os) const { os<<"sgn("; _u.print(os); os <<")"; }
};


template<class E>
class pow_expr : public _expr<pow_expr<E> >
{
    typename E::Nested_t _u;

public:
    typedef typename E::Scalar Scalar;
    enum {ScalarValued = 1 };
    enum {Space = E::Space};

    Scalar _q;// power

    pow_expr(_expr<E> const& u, Scalar q) : _u(u), _q(q) { }

    Scalar eval(const index_t k) const
    {
        const Scalar v = _u.val().eval(k);
        return math::pow(v,_q);
    }

    static index_t rows() { return 0; }
    static index_t cols() { return 0; }

    void setFlag() const { _u.setFlag(); }

    void parse(gsSortedVector<const gsFunctionSet<Scalar>*> & ) const {  }

    static bool isScalar() { return true; }

    enum{rowSpan = 0, colSpan = 0};

    const gsFeSpace<Scalar> & rowVar() const {return gsNullExpr<Scalar>::get();}
    const gsFeSpace<Scalar> & colVar() const {return gsNullExpr<Scalar>::get();}

    void print(std::ostream &os) const { os<<"pow("; _u.print(os); os <<")"; }
};

// Call as pow(a,b)
template<class E>
pow_expr<E> pow(_expr<E> const& u, real_t q) { return pow_expr<E>(u,q); }


/**
computes outer products of a matrix by a space of dimension > 1
[Jg Jg Jg] * Jb ..
(d x d^2)  * (d^2 x N*d)  --> (d x N*d)
*/
template <typename E1, typename E2>
class matrix_by_space_expr  : public _expr<matrix_by_space_expr<E1,E2> >
{
public:
    typedef typename E1::Scalar Scalar;
    enum {ScalarValued = 0, ColBlocks = 1};
    enum {Space = E2::Space};
private:
    typename E1::Nested_t _u;
    typename E2::Nested_t _v;
    mutable gsMatrix<Scalar> res;

public:
    matrix_by_space_expr(E1 const& u, E2 const& v) : _u(u), _v(v) { }


    // choose if ColBlocks
    const gsMatrix<Scalar> & eval(const index_t k) const
    {
        const index_t r   = _u.rows();
        const index_t N  = _v.cols() / (r*r);

        const MatExprType uEv        = _u.eval(k);
        const MatExprType vEv  = _v.eval(k);

        res.resize(r, N*r*r);
        // gsDebugVar(res.cols());
        for (index_t s = 0; s!=r; ++s)
            for (index_t i = 0; i!=N; ++i)
            {
                res.middleCols((s*N + i)*r,r).noalias() =
                uEv.col(s) * vEv.middleCols((s*N + i)*r,r).row(s);
                //uEv*vEv.middleCols((s*N + i)*r,r);
            }
        //meaning: [Jg Jg Jg] * Jb ..
        return res;
    }

    index_t rows() const { return _u.cols(); }
    index_t cols() const { return _v.cols(); }
    void setFlag() const { _u.setFlag(); _v.setFlag(); }

    void parse(gsSortedVector<const gsFunctionSet<Scalar>*> & evList) const
    { _u.parse(evList); }

    const gsFeSpace<Scalar> & rowVar() const { return _v.rowVar(); }
    const gsFeSpace<Scalar> & colVar() const { return _v.colVar(); }

    enum{rowSpan = E2::rowSpan, colSpan = E2::colSpan};

    void print(std::ostream &os) const { os << "matrix_by_space("; _u.print(os); os<<")"; }
};

/**
computes outer products of a matrix by a space of dimension > 1
[Jg Jg Jg] * Jb ..
(d x d^2)  * (d^2 x N*d)  --> (d x N*d)
*/
template <typename E1, typename E2>
class matrix_by_space_expr_tr  : public _expr<matrix_by_space_expr_tr<E1,E2> >
{
public:
    typedef typename E1::Scalar Scalar;
    enum {ScalarValued = 0, ColBlocks = 1};
    enum {Space = E2::Space};
private:
    typename E1::Nested_t _u;
    typename E2::Nested_t _v;
    mutable gsMatrix<Scalar> res;

public:
    matrix_by_space_expr_tr(E1 const& u, E2 const& v) : _u(u), _v(v) { }


    // choose if ColBlocks
    const gsMatrix<Scalar> & eval(const index_t k) const
    {
        const index_t r   = _u.rows();
        const index_t N  = _v.cols() / (r*r);

        const MatExprType uEv        = _u.eval(k);
        const MatExprType vEv  = _v.eval(k);

        res.resize(r, N*r*r);
        // gsDebugVar(res.cols());
        for (index_t s = 0; s!=r; ++s)
            for (index_t i = 0; i!=N; ++i)
            {
                res.middleCols((s*N + i)*r,r).noalias() =
                //uEv.col(s) * vEv.middleCols((s*N + i)*r,r).row(s);
                uEv.transpose()*vEv.middleCols((s*N + i)*r,r).transpose();
            }
        //meaning: [Jg Jg Jg] * Jb ..
        return res;
    }

    index_t rows() const { return _u.cols(); }
    index_t cols() const { return _v.cols(); }
    void setFlag() const { _u.setFlag(); _v.setFlag(); }

    void parse(gsSortedVector<const gsFunctionSet<Scalar>*> & evList) const
    { _u.parse(evList); }

    const gsFeSpace<Scalar> & rowVar() const { return _v.rowVar(); }
    const gsFeSpace<Scalar> & colVar() const { return _v.colVar(); }

    enum{rowSpan = E2::rowSpan, colSpan = E2::colSpan};

    void print(std::ostream &os) const { os << "matrix_by_space_tr("; _u.print(os); os<<")"; }
};

/*
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

    enum{rowSpan = 0, colSpan = 0};

    const gsFeSpace<Scalar> & rowVar() const {return gsNullExpr<Scalar>::get();}
    const gsFeSpace<Scalar> & colVar() const {return gsNullExpr<Scalar>::get();}

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

/*
   Expression for the gradient of a finite element variable

   Transposed gradient vectors are returned as a matrix
 */
template<class E>
class grad_expr : public _expr<grad_expr<E> >
{
    typename E::Nested_t _u;
public:
    enum{ Space = E::Space };

    typedef typename E::Scalar Scalar;
    mutable gsMatrix<Scalar> tmp;

    grad_expr(const E & u) : _u(u)
    { GISMO_ASSERT(1==u.dim(),"grad(.) requires 1D variable, use jac(.) instead.");}

    const gsMatrix<Scalar> & eval(const index_t k) const
    {
        tmp = _u.data().values[1].reshapeCol(k, cols(), rows()).transpose();
        return tmp;
    }

    index_t rows() const
    {
        //return _u.data().values[0].rows();
        return _u.data().values[1].rows() / cols();
    }
    //index_t rows() const { return _u.data().actives.size(); }
    //index_t rows() const { return _u.rows(); }

    //index_t rows() const { return _u.source().targetDim() is wrong }
    index_t cols() const { return _u.source().domainDim(); }

    void setFlag() const
    {
        _u.data().flags |= NEED_GRAD;
        if (_u.composed() )
            _u.mapData().flags |= NEED_VALUE;
    }

    void parse(gsSortedVector<const gsFunctionSet<Scalar>*> & evList) const
    {
        //GISMO_ASSERT(NULL!=m_fd, "FeVariable: FuncData member not registered");
        evList.push_sorted_unique(&_u.source());
        _u.data().flags |= NEED_GRAD;
        if (_u.composed() )
            _u.mapData().flags |= NEED_VALUE;
    }

    const gsFeSpace<Scalar> & rowVar() const { return _u.rowVar(); }
    const gsFeSpace<Scalar> & colVar() const
    {return gsNullExpr<Scalar>::get();}

    enum{rowSpan = E::rowSpan, colSpan = 0};

    void print(std::ostream &os) const { os << "grad("; _u.print(os); os <<")"; }
};

/*
   Expression for the derivative of the jacobian of a spline geometry map,
   with respect to the coordinate c.

   It returns a matrix with the gradient of u in row d.
 */
template<class E>
class dJacdc_expr : public _expr<dJacdc_expr<E> >
{
    typename E::Nested_t _u;
public:
    enum{ Space = E::Space, ScalarValued = 0, ColBlocks = E::rowSpan};

    typedef typename E::Scalar Scalar;

    mutable gsMatrix<Scalar> res;
    index_t _c;

    dJacdc_expr(const E & u, index_t c) : _u(u), _c(c)
    { GISMO_ASSERT(1==u.dim(),"grad(.) requires 1D variable, use jac(.) instead.");}

    const gsMatrix<Scalar> & eval(const index_t k) const
    {
        index_t dd = _u.source().domainDim();
        index_t n = _u.rows();
        res.setZero(dd, dd*n);

        gsMatrix<Scalar> grad = _u.data().values[1].reshapeCol(k, dd, n);
        for(index_t i = 0; i < n; i++){
            res.row(_c).segment(i*dd,dd) = grad.col(i);
        }
        return res;
    }

    index_t rows() const
    {
        //return _u.data().values[0].rows();
        return _u.source().domainDim();
    }
    //index_t rows() const { return _u.data().actives.size(); }
    //index_t rows() const { return _u.rows(); }

    //index_t rows() const { return _u.source().targetDim() is wrong }
    index_t cols() const { return _u.source().domainDim()*_u.rows(); }

    void setFlag() const
    {
        _u.data().flags |= NEED_GRAD;
        if (_u.composed() )
            _u.mapData().flags |= NEED_VALUE;
    }

    void parse(gsSortedVector<const gsFunctionSet<Scalar>*> & evList) const
    {
        //GISMO_ASSERT(NULL!=m_fd, "FeVariable: FuncData member not registered");
        evList.push_sorted_unique(&_u.source());
        _u.data().flags |= NEED_GRAD;
        if (_u.composed() )
            _u.mapData().flags |= NEED_VALUE;
    }

    const gsFeSpace<Scalar> & rowVar() const { return _u.rowVar(); }
    const gsFeSpace<Scalar> & colVar() const
    {return gsNullExpr<Scalar>::get();}

    enum{rowSpan = E::rowSpan, colSpan = 0};

    void print(std::ostream &os) const { os << "dJacdc("; _u.print(os); os <<")"; }
};


/*
   Expression for the nabla (\f$\nabla\f$) of a finite element variable,
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

    const gsFeSpace<T> & rowVar() const { return u.rowVar(); }
    const gsFeSpace<T> & colVar() const
    {return gsNullExpr<T>();}

    enum{rowSpan = 1, colSpan = 0};

    void print(std::ostream &os) const { os << "nabla("; u.print(os); os <<")"; }
};

/*
   Expression for the nabla2 (\f$\nabla^2\f$ or Del2) of a finite element variable,
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

    const gsFeSpace<T> & rowVar() const { return u.rowVar(); }
    const gsFeSpace<T> & colVar() const
    {return gsNullExpr<T>::get();}

    enum{rowSpan = 1, colSpan = 0};
};

/// The nabla2 (\f$\nabla^2\f$) of a finite element variable
template<class T>
nabla2_expr<T> nabla2(const gsFeVariable<T> & u)
{ return nabla2_expr<T>(u); }
// #define lapl(x) nabla2(x).sum() // assume tarDim==1

/**
   Expression for the outer pointing normal of a geometry map. This
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

    const gsFeSpace<T> & rowVar() const {return gsNullExpr<T>::get();}
    const gsFeSpace<T> & colVar() const {return gsNullExpr<T>::get();}

    enum{rowSpan = 0, colSpan = 0};

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

/**
   Expression for the out of plane surface normal of a geometry map.
   The expression is valid for a surface or hypersurface.
 */
template<class T>
class normal_expr : public _expr<normal_expr<T> >
{
    typename gsGeometryMap<T>::Nested_t _G;

public:
    typedef T Scalar;

    normal_expr(const gsGeometryMap<T> & G) : _G(G) { }

    MatExprType eval(const index_t k) const
    {
        return _G.data().normals.col(k);
    }

    index_t rows() const { return _G.data().dim.second; }
    index_t cols() const { return 1; }

    enum{rowSpan = 0, colSpan = 0};

    const gsFeSpace<T> & rowVar() const {return gsNullExpr<T>::get();}
    const gsFeSpace<T> & colVar() const {return gsNullExpr<T>::get();}

    void setFlag() const { _G.data().flags |= NEED_NORMAL; }

    void parse(gsSortedVector<const gsFunctionSet<Scalar>*> & evList) const
    {
        //GISMO_ASSERT(NULL!=m_fd, "FeVariable: FuncData member not registered");
        evList.push_sorted_unique(&_G.source());
        _G.data().flags |= NEED_NORMAL;
    }

    // Normalized to unit length
    normalized_expr<normal_expr<T> > normalized()
    { return normalized_expr<normal_expr<T> >(*this); }

    void print(std::ostream &os) const { os << "sn("; _G.print(os); os <<")"; }
};

/**
   Expression for the tangent vector of a geometry map. This
   expression is valid only at the boundaries of a geometric patch and only for in-plane surfaces
 */
template<class T>
class tangent_expr : public _expr<tangent_expr<T> >
{
    typename gsGeometryMap<T>::Nested_t _G;

public:
    typedef T Scalar;

    tangent_expr(const gsGeometryMap<T> & G) : _G(G) { }

    mutable gsMatrix<Scalar> res;
    mutable gsVector<Scalar,3> normal,onormal,otangent;

    const gsMatrix<Scalar> eval(const index_t k) const
    {
        if (_G.targetDim()==2)
        {
            res = _G.data().outNormals.col(k);//2x1
            std::swap( res(0,0), res(1,0) );
            res(0,0) *= -1;
            return res;
        }
        else if (_G.targetDim()==3)
        {
            onormal = _G.data().outNormals.col(k);
            normal =  _G.data().normals.col(k);
            otangent = normal.template head<3>().cross(onormal.template head<3>());
            return otangent;
        }
        else
            GISMO_ERROR("Function not implemented for dimension"<<_G.targetDim());

    }

    index_t rows() const { return _G.data().dim.second; }
    index_t cols() const { return 1; }

    enum{rowSpan = 0, colSpan = 0};

    void setFlag() const
    {
        _G.data().flags |= NEED_NORMAL;
        _G.data().flags |= NEED_OUTER_NORMAL;
    }

    void parse(gsSortedVector<const gsFunctionSet<Scalar>*> & evList) const
    {
        //GISMO_ASSERT(NULL!=m_fd, "FeVariable: FuncData member not registered");
        evList.push_sorted_unique(&_G.source());
        _G.data().flags |= NEED_NORMAL;
        _G.data().flags |= NEED_OUTER_NORMAL;
    }

    // Normalized to unit length
    normalized_expr<tangent_expr<T> > normalized()
    { return normalized_expr<tangent_expr<T> >(*this); }

    void print(std::ostream &os) const { os << "tv("; _G.print(os); os <<")"; }
};

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

    index_t rows() const { return _u.data().laplacians.rows(); }
    index_t cols() const { return 1; }

    enum{rowSpan = 1, colSpan = 0};

    void setFlag() const { _u.data().flags |= NEED_LAPLACIAN; }

    void parse(gsSortedVector<const gsFunctionSet<Scalar>*> & evList) const
    {
        //GISMO_ASSERT(NULL!=m_fd, "FeVariable: FuncData member not registered");
        evList.push_sorted_unique(&_u.source());
        _u.data().flags |= NEED_LAPLACIAN;
    }

    void print(std::ostream &os) const { os << "lap("; _u.print(os); os <<")"; }
};

/*
   Expression for the Laplacian of a finite element solution
 */
template<class T>
class solLapl_expr : public _expr<solLapl_expr<T> >
{
protected:
    const gsFeSolution<T> & _u;

public:
    typedef T Scalar;

    solLapl_expr(const gsFeSolution<T> & u) : _u(u) { }

    mutable gsMatrix<T> res;
    const gsMatrix<T> eval(const index_t k) const
    {
        GISMO_ASSERT(1==_u.data().actives.cols(), "Single actives expected");

        res.setZero(_u.dim(), 1); //  scalar, but per component
        const gsDofMapper & map = _u.mapper();

        index_t numActs = _u.data().values[0].rows();
        index_t numDers = _u.parDim() * (_u.parDim() + 1) / 2;
        gsMatrix<T> deriv2;

        for (index_t c = 0; c!= _u.dim(); c++)
            for (index_t i = 0; i!=numActs; ++i)
            {
                const index_t ii = map.index(_u.data().actives.at(i), _u.data().patchId,c);
                deriv2 = _u.data().values[2].block(i*numDers,k,_u.parDim(),1); // this only takes d11, d22, d33 part. For all the derivatives [d11, d22, d33, d12, d13, d23]: col.block(i*numDers,k,numDers,1)
                if ( map.is_free_index(ii) ) // DoF value is in the solVector
                    res.at(c) += _u.coefs().at(ii) * deriv2.sum();
                else
                    res.at(c) +=_u.fixedPart().at( map.global_to_bindex(ii) ) * deriv2.sum();
            }
        return res;
    }

    index_t rows() const { return _u.dim(); }
    index_t cols() const { return 1; }

    enum{rowSpan = 1, colSpan = 0};

    void setFlag() const { _u.data().flags |= NEED_DERIV2; }

    void parse(gsSortedVector<const gsFunctionSet<Scalar>*> & evList) const
    {
        //GISMO_ASSERT(NULL!=m_fd, "FeVariable: FuncData member not registered");
        evList.push_sorted_unique(&_u.source());
        _u.data().flags |= NEED_DERIV2;
    }

    void print(std::ostream &os) const { os << "lap(s)"; }
};

template<class T>
class solHess_expr : public _expr<solHess_expr<T> >
{
protected:
    const gsFeSolution<T> & _u;

public:
    typedef T Scalar;

    solHess_expr(const gsFeSolution<T> & u) : _u(u) { }

    mutable gsMatrix<T> res;
    const gsMatrix<T> eval(const index_t k) const
    {
        GISMO_ASSERT(1==_u.data().actives.cols(), "Single actives expected");

        const gsDofMapper & map = _u.mapper();

        index_t numActs = _u.data().values[0].rows();
        index_t numDers = cols();
        gsMatrix<T> deriv2;

        res.setZero(rows(), cols());
        for (index_t c = 0; c!= _u.dim(); c++)
            for (index_t i = 0; i!=numActs; ++i)
            {
                const index_t ii = map.index(_u.data().actives.at(i), _u.data().patchId,c);
                deriv2 = _u.data().values[2].block(i*numDers,k,numDers,1).transpose(); // start row, start col, rows, cols
                if ( map.is_free_index(ii) ) // DoF value is in the solVector
                    res.row(c) += _u.coefs().at(ii) * deriv2;
                else
                    res.row(c) +=_u.fixedPart().at( map.global_to_bindex(ii) ) * deriv2;
            }
        return res;

    }

    index_t rows() const
    {
        return _u.dim(); //  number of components
    }
    index_t cols() const
    {// second derivatives in the columns; i.e. [d11, d22, d33, d12, d13, d23]
        return _u.parDim() * (_u.parDim() + 1) / 2;
    }

    enum{rowSpan = 1, colSpan = 0};

    void setFlag() const { _u.data().flags |= NEED_DERIV2; }

    void parse(gsSortedVector<const gsFunctionSet<Scalar>*> & evList) const
    {
        //GISMO_ASSERT(NULL!=m_fd, "FeVariable: FuncData member not registered");
        evList.push_sorted_unique(&_u.source());
        _u.data().flags |= NEED_DERIV2;
    }

    void print(std::ostream &os) const { os << "hess(s)"; }
};



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

    MatExprType eval(const index_t k) const
    {   // todo: fix funcdata
        //return _G.data().fundForm(k).transpose() * _G.data().fundForm(k) ;
        GISMO_NO_IMPLEMENTATION
    }

    index_t rows() const { return _G.data().dim.second; }
    index_t cols() const { return _G.data().dim.first ; }

    enum{rowSpan = 0, colSpan = 0};

    void setFlag() const { _G.data().flags |= NEED_GRAD_TRANSFORM; }

    void parse(gsSortedVector<const gsFunctionSet<Scalar>*> & evList) const
    {
        //GISMO_ASSERT(NULL!=m_fd, "FeVariable: FuncData member not registered");
        evList.push_sorted_unique(&_G.source());
        _G.data().flags |= NEED_GRAD_TRANSFORM;
    }

    const gsFeSpace<Scalar> & rowVar() const {return gsNullExpr<T>::get();}
    const gsFeSpace<Scalar> & colVar() const {return gsNullExpr<T>::get();}

    void print(std::ostream &os) const { os << "fform("; _G.print(os); os <<")"; }
};

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

    enum{rowSpan = 0, colSpan = 0};

    void setFlag() const { _G.data().flags |= NEED_GRAD_TRANSFORM; }

    void parse(gsSortedVector<const gsFunctionSet<Scalar>*> & evList) const
    {
        //GISMO_ASSERT(NULL!=m_fd, "FeVariable: FuncData member not registered");
        evList.push_sorted_unique(&_G.source());
        _G.data().flags |= NEED_GRAD_TRANSFORM;
    }

    const gsFeSpace<Scalar> & rowVar() const {return gsNullExpr<T>::get();}
    const gsFeSpace<Scalar> & colVar() const {return gsNullExpr<T>::get();}

    // todo mat_expr ?
    // tr() const --> _G.data().fundForm(k)

    void print(std::ostream &os) const { os << "jacInv("; _G.print(os); os <<")"; }
};

// #define GISMO_FUNCDATA_EXPRESSION( ...

/*
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

    index_t rows() const
    {
        return _G.source().targetDim();
        //return _G.data().dim.second;
    }

    index_t cols() const
    {
        return _G.source().domainDim();
        //return _G.data().dim.first;
    }

    enum{rowSpan = 0, colSpan = 0};

    static const gsFeSpace<Scalar> & rowVar() { return gsNullExpr<Scalar>::get(); }
    static const gsFeSpace<Scalar> & colVar() { return gsNullExpr<Scalar>::get(); }

    void setFlag() const { _G.data().flags |= NEED_DERIV; }

    void parse(gsSortedVector<const gsFunctionSet<Scalar>*> & evList) const
    {
        //GISMO_ASSERT(NULL!=m_fd, "FeVariable: FuncData member not registered");
        evList.push_sorted_unique(&_G.source());
        _G.data().flags |= NEED_DERIV;
    }

    meas_expr<T> absDet() const
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


/*
   Expression for the Jacobian matrix of a FE variable
 */
template<class E>
class jac_expr : public _expr<jac_expr<E> >
{
    typename E::Nested_t m_fev;
public:
    enum {ColBlocks = E::rowSpan };
    enum {Space = E::Space };

    typedef typename E::Scalar Scalar;

    mutable gsMatrix<Scalar> res;

    jac_expr(const E & _u)
    : m_fev(_u) { }

    MatExprType eval(const index_t k) const
    {
        // Dim x (numActive*Dim)
        res = m_fev.data().values[1].col(k).transpose().blockDiag(m_fev.dim());
        return res;
    }

    const gsFeSpace<Scalar> & rowVar() const { return m_fev; }
    const gsFeSpace<Scalar> & colVar() const { return gsNullExpr<Scalar>::get(); }

    index_t rows() const { return m_fev.dim(); }
    index_t cols() const
    {   //bug: Should return the column size of each BLOCK contained
        // return m_fev.dim() * m_fev.data().actives.rows() * m_fev.data().dim.first;
        return m_fev.data().dim.first;
    }

    index_t cardinality_impl() const
    {
        return m_fev.dim() * m_fev.data().actives.rows();
    }

    enum{rowSpan = 1, colSpan = 0};

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

    void print(std::ostream &os) const { os << "jac("; m_fev.print(os);os <<")"; }
};

/*
   Expression for the Jacobian matrix of a vector function
 */
template<class T>
class fjac_expr : public _expr<fjac_expr<T> >
{
    const gsFeVariable<T> & m_fev;
public:
    enum {ColBlocks = 0};// main difference from jac

    typedef T Scalar;

    fjac_expr(const gsFeVariable<T> & _u)
    : m_fev(_u)
    {
        gsInfo<<"(!) fjac(u) \n";
    }

    MatExprType eval(const index_t k) const
    {
        return m_fev.data().values[1].reshapeCol(k, cols(), rows()).transpose();
    }

    const gsFeSpace<T> & rowVar() const { return gsNullExpr<T>::get(); }
    const gsFeSpace<T> & colVar() const { return gsNullExpr<T>::get(); }

    index_t rows() const { return m_fev.data().dim.first; }
    index_t cols() const {return m_fev.data().dim.second; }

    enum{rowSpan = 1, colSpan = 0};

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

    void print(std::ostream &os) const { os << "fjac("; m_fev.print(os);os <<")"; }
};

/*
   Expression for the Hessian matrix of (a coordinate of) the geometry
   map or a finite element variable

template<class T>
class hess_expr : public _expr<hess_expr<T> >
{
public:
    typedef T Scalar;
    enum {ScalarValued = 0, ColBlocks = 1};

private:
    const gsFuncData<T> * m_data;
    mutable gsMatrix<Scalar> res;

    //hess_expr(const hess_expr & );
public:
    hess_expr(const gsGeometryMap<T> & G)
    : m_data(&G.data()) { } //ColBlocks=0 ?

    hess_expr(const gsFeVariable<T> & _u)
    : m_data(&_u.data())
    {
        GISMO_ASSERT(1==_u.dim(),"hess(.) requires 1D variable");
    }

    const gsMatrix<Scalar> & eval(const index_t k) const
    {
        const index_t sz = cols();
        res.resize(m_data->dim.first, sz*m_data->dim.first);
        secDerToHessian(m_data->values[2].col(k), m_data->dim.first, res);
        res.resize(m_data->dim.first, res.cols()*m_data->dim.first);
        // Note: auto returns by value here
        return res;
    }

    index_t rows() const
    {
        // gsDebugVar(m_data);
        // gsDebugVar(m_data->dim.first);
        // gsDebugVar(m_data->values[0].rows());
        return m_data->dim.first;
    }
    index_t cols() const
    {
        return 2*m_data->values[2].rows() / (1+m_data->dim.first);
    }

    void setFlag() const { m_data->flags |= NEED_2ND_DER; }

    void parse(gsSortedVector<const gsFunctionSet<Scalar>*> & evList) const
    {
        GISMO_ERROR("error 1712");
    }

    enum{rowSpan = 1, colSpan = 0};

    const gsFeVariable<T> & rowVar() const
    {
        GISMO_ERROR("error -- hess");
        return gsNullExpr<T>::get();
    }
    const gsFeSpace<T> & colVar() const { return gsNullExpr<T>::get(); }

    void print(std::ostream &os) const
    //    { os << "hess("; _u.print(os);os <<")"; }
    { os << "hess(U)"; }
};
 */

template<class E>
class hess_expr : public _expr<hess_expr<E> >
{
public:
    typedef typename E::Scalar Scalar;
private:
    typename E::Nested_t _u;
    mutable gsMatrix<Scalar> res;
public:
    enum {ScalarValued = 0, ColBlocks = E::rowSpan };
    enum {Space = E::Space };

public:
    hess_expr(const E & u) : _u(u)
    {
        //gsInfo << "expression is space ? "<<E::Space <<"\n"; _u.print(gsInfo);
        //GISMO_ASSERT(1==_u.dim(),"hess(.) requires 1D variable");
    }

    const gsMatrix<Scalar> & eval(const index_t k) const
    {
        const gsFuncData<Scalar> & dd = _u.data();
        const index_t sz = cols();
        res.resize(dd.dim.first, sz*dd.dim.first);
        secDerToHessian(dd.values[2].col(k), dd.dim.first, res);
        res.resize(dd.dim.first, res.cols()*dd.dim.first);
        // Note: auto returns by value here
        return res;
    }

    index_t rows() const
    {
        // gsDebugVar(m_data);
        // gsDebugVar(_u.data().dim.first);
        // gsDebugVar(_u.data().values[0].rows());
        return _u.data().dim.first;
    }
    index_t cols() const
    {
        return 2*_u.data().values[2].rows() / (1+_u.data().dim.first);
    }

    void setFlag() const { _u.data().flags |= NEED_2ND_DER; }

    void parse(gsSortedVector<const gsFunctionSet<Scalar>*> & evList) const
    {
        GISMO_ERROR("error 1712");
    }

    enum{rowSpan = E::rowSpan, colSpan = 0};

    const gsFeSpace<Scalar> & rowVar() const { return _u.rowVar(); }
    const gsFeSpace<Scalar> & colVar() const { return gsNullExpr<Scalar>::get(); }

    void print(std::ostream &os) const
    //    { os << "hess("; _u.print(os);os <<")"; }
    { os << "hess(U)"; }
};


/*
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

    enum{rowSpan = 0, colSpan = 0};
};


/*
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

    T eval(const index_t k) const
    {
        //gsDebugVar(_G.data().measures.at(k));
        return _G.data().measures.at(k);
    }

    index_t rows() const { return 0; }
    index_t cols() const { return 0; }
    void setFlag() const { _G.data().flags |= NEED_MEASURE; }
    void parse(gsSortedVector<const gsFunctionSet<Scalar>*> & evList) const
    {
        //GISMO_ASSERT(NULL!=m_fd, "FeVariable: FuncData member not registered");
        evList.push_sorted_unique(&_G.source());
        _G.data().flags |= NEED_MEASURE;
    }

    const gsFeSpace<T> & rowVar() const { return gsNullExpr<T>::get(); }
    const gsFeSpace<T> & colVar() const { return gsNullExpr<T>::get(); }

    enum{rowSpan = 0, colSpan = 0};

    void print(std::ostream &os) const { os << "meas("; _G.print(os); os <<")"; }
};

/*
   Expression for the curl
*/
template<class T>
class curl_expr : public _expr<curl_expr<T> >
{
public:
    typedef T Scalar;
private:
    typename gsFeVariable<T>::Nested_t _u;
    mutable gsMatrix<Scalar> res;
public:
    curl_expr(const gsFeVariable<T> & u) : _u(u)
    { GISMO_ASSERT(3==u.dim(),"curl(.) requires 3D variable."); }

    MatExprType eval(const index_t k) const
    {
        res.setZero( rows(), _u.dim());
        const index_t na = _u.data().values[0].rows();
        gsAsConstMatrix<T, Dynamic, Dynamic> pd =
            _u.data().values[1].reshapeCol(k, cols(), na);

        res.col(0).segment(na  ,na) = -pd.row(2);
        res.col(0).segment(2*na,na) =  pd.row(1);
        res.col(1).segment(0   ,na) =  pd.row(2);
        res.col(1).segment(2*na,na) = -pd.row(0);
        res.col(2).segment(0   ,na) = -pd.row(1);
        res.col(2).segment(na  ,na) =  pd.row(0);
        return res;
    }

    index_t rows() const { return _u.dim() * _u.data().values[0].rows(); }
    index_t cols() const { return _u.data().dim.first; }
    void setFlag() const
    {
        _u.data().flags |= NEED_GRAD;
        if (_u.composed() )
            _u.mapData().flags |= NEED_VALUE;
    }

    void parse(gsSortedVector<const gsFunctionSet<Scalar>*> & evList) const
    {
        //GISMO_ASSERT(NULL!=m_fd, "FeVariable: FuncData member not registered");
        evList.push_sorted_unique(&_u.source());
        _u.data().flags |= NEED_GRAD;
        if (_u.composed() )
            _u.mapData().flags |= NEED_VALUE;
    }

    const gsFeSpace<T> & rowVar() const { return _u.rowVar(); }
    const gsFeSpace<T> & colVar() const {return gsNullExpr<T>::get();}

    enum{rowSpan = 1, colSpan = 0};

    void print(std::ostream &os) const { os << "curl("; _u.print(os); os <<")"; }
};

/*
   Expression for multiplication operation (first version)

   First argument E1 has ColBlocks = false

   Partial specialization for (right) blockwise multiplication

   B * [A1 A2 A3] = [B*A1  B*A2  B*A3]

 */
template <typename E1, typename E2>
class mult_expr<E1,E2,false> : public _expr<mult_expr<E1, E2, false> >
{
    typename E1::Nested_t _u;
    typename E2::Nested_t _v;

public:
    enum {ScalarValued = E1::ScalarValued && E2::ScalarValued,
          ColBlocks = E2::ColBlocks};
    enum {Space = ( 0!=E1::Space ? int(E1::Space) : int(E2::Space) ) };

    typedef typename E1::Scalar Scalar;

    typedef typename
    util::conditional<ScalarValued,Scalar,typename gsMatrix<Scalar>::Base >::type Temporary_t;

    mutable Temporary_t tmp;

    mult_expr(_expr<E1> const& u,
              _expr<E2> const& v)
    : _u(u), _v(v) { }

    //EIGEN_STRONG_INLINE MatExprType
    const Temporary_t &
    eval(const index_t k) const
    {
        GISMO_ASSERT(0==_u.cols()*_v.rows() || _u.cols() == _v.rows(),
                     "Wrong dimensions "<<_u.cols()<<"!="<<_v.rows()<<" in * operation:\n"
                     << _u <<" times \n" << _v );
        // Note: a * b * c --> (a*b).eval()*c
        // gsDebugVar( _u.eval(k) );
        // gsDebugVar( _v.eval(k) );

        tmp = _u.eval(k) * _v.eval(k);
        // gsDebugVar(_u.eval(k));
        // gsDebugVar(_v.eval(k));
        // gsDebugVar(tmp);
        // gsDebugVar(E1::Space);
        // gsDebugVar(E2::Space);
        // gsDebugVar(Space);
        // gsDebugVar(E1::rowSpan);
        // gsDebugVar(E1::colSpan);
        // gsDebugVar(E2::rowSpan);
        // gsDebugVar(E2::colSpan);
        // gsDebugVar(rowSpan);
        // gsDebugVar(colSpan);
        // gsInfo<<"expression: "; _u.print(gsInfo); gsInfo<<"\n";
        // gsInfo<<"expression: "; _v.print(gsInfo); gsInfo<<"\n";
        return tmp; // assume result not scalarv
    }

    index_t rows() const { return E1::ScalarValued ? _v.rows()  : _u.rows(); }
    index_t cols() const { return E2::ScalarValued ? _u.cols()  : _v.cols(); }
    void setFlag() const { _u.setFlag(); _v.setFlag(); }
    void parse(gsSortedVector<const gsFunctionSet<Scalar>*> & evList) const
    { _u.parse(evList); _v.parse(evList); }

    enum{rowSpan = (0==E1::Space ? E2::rowSpan : int(E1::rowSpan)), colSpan = (0==E2::Space ? E1::colSpan : int(E2::colSpan)) };

    index_t cardinality_impl() const
    { return 0==E1::Space ? _v.cardinality(): _u.cardinality(); }

    const gsFeSpace<Scalar> & rowVar() const
    { return 0==E1::Space ? _v.rowVar() : _u.rowVar(); }
    const gsFeSpace<Scalar> & colVar() const
    {
        return 0==E2::Space ? _u.colVar() : _v.colVar();
        /*
        if ( _v.isScalar() )
            return _u.colVar();
        else
            return _v.colVar();
        */
    }

    void print(std::ostream &os) const { _u.print(os); os<<"*"; _v.print(os); }
};

/*
  Expression for multiplication operation (second version)

   Partial specialization for (right) blockwise multiplication
   [A1 A2 A3] * B = [A1*B  A2*B  A3*B]

   as well as

   [A1 A2 A3] * [B1 B2 B3] = [A1*B1  A2*B2  A3*B3]

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
    enum {Space = E1::Space};

    mult_expr(_expr<E1> const& u,
              _expr<E2> const& v)
    : _u(u), _v(v)
    {

    }

    const gsMatrix<Scalar> & eval(const index_t k) const
    {
        const index_t uc = _u.cols();
        const index_t ur = _u.rows();
        const index_t nb = _u.cardinality();
        const MatExprType tmpA = _u.eval(k);
        const MatExprType tmpB = _v.eval(k);

        // gsDebugVar(_u.eval(k));
        // gsDebugVar(_v.eval(k));
        // gsDebugVar(E1::rowSpan);
        // gsDebugVar(E1::colSpan);
        // gsDebugVar(E2::rowSpan);
        // gsDebugVar(E2::colSpan);
        // gsDebugVar(rowSpan);
        // gsDebugVar(colSpan);
        // gsDebugVar(rowVar());
        // gsDebugVar(colVar());
        // gsDebugVar(_u.rowVar());
        // gsDebugVar(_u.colVar());
        // gsDebugVar(_v.rowVar());
        // gsDebugVar(_v.colVar());
        // gsInfo<<"expression: "; _u.print(gsInfo); gsInfo<<"\n";
        // gsInfo<<"expression: "; _v.print(gsInfo); gsInfo<<"\n";

        // gsDebugVar(tmpA);
        // gsDebugVar(tmpB);


        // gsDebugVar(_u.ColBlocks);
        // gsDebugVar(_v.ColBlocks);

        //gsDebugVar(_v.cols());
        //gsDebugVar(ur);

        //GISMO_ASSERT(uc==vr, "Assumption for product failed");
        const index_t vc = _v.cols();
        // gsDebugVar( _u.cardinality() );
        // gsDebugVar( _v.cardinality() );

        // either _v.cardinality()==1 or _v.cardinality()==_u.cardinality()
        if (  1 == _v.cardinality() ) //second is not ColBlocks
        {
            res.resize(ur, vc*nb);
            GISMO_ASSERT(1==_v.cardinality(), "Dimension error");
            //gsInfo<<"cols = "<<res.cols()<<"; rows = "<<res.rows()<<"\n";
            for (index_t i = 0; i!=nb; ++i)
                res.middleCols(i*vc,vc).noalias()
                    = tmpA.middleCols(i*uc,uc) * tmpB;
        }
        // both are ColBlocks: [A1 A2 A3] * [B1 B2 B3] = [A1*B1  A2*B2  A3*B3]
        //                                               [A2*B1 ..           ]
        //                                               [                   ]
        else
        {
            const index_t nbv = _v.cardinality();
            res.resize(ur*nb, vc*nbv);
            for (index_t i = 0; i!=nb; ++i)
                for (index_t j = 0; j!=nbv; ++j)
                {
                    res.block(i*ur,j*vc,ur,vc).noalias() =
                        tmpA.middleCols(i*uc,uc) * tmpB.middleCols(j*vc,vc);

                // res.middleCols(i*vc,vc).noalias()
                //     = tmpA.middleCols(i*uc,uc) * tmpB.middleCols(i*vc,vc);
            }
        }


        return res;
    }

    index_t rows() const {
        return _u.rows();
    }
    index_t cols() const {
        // DEBUG changed by asgl, perhaps there was a bug here?
        //return _v.cols() * (_u.cols()/_u.rows());
        return _u.cols();
    }
    void setFlag() const { _u.setFlag(); _v.setFlag(); }
    void parse(gsSortedVector<const gsFunctionSet<Scalar>*> & evList) const
    { _u.parse(evList); _v.parse(evList); }

    enum{rowSpan = E1::rowSpan || E2::rowSpan,colSpan = E2::colSpan || E1::colSpan};

    index_t cardinality_impl() const { return  _u.cardinality(); }

    const gsFeSpace<Scalar> & rowVar() const { return _u.rowVar(); }
    const gsFeSpace<Scalar> & colVar() const
    {
        if ( 1 == _v.cardinality() )
            return _u.colVar();
        else
            return _v.colVar();
    }

    void print(std::ostream &os) const
    { os << "("; _u.print(os);os <<"*";_v.print(os);os << ")"; }
};

/*
  Expression for multiplication operation (third version)

   Scalar multiplication
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
    enum {Space = E2::Space};

    mult_expr(Scalar const & c, _expr<E2> const& v)
    : _c(c), _v(v) { }

    EIGEN_STRONG_INLINE AutoReturn_t eval(const index_t k) const
    {
        return ( _c * _v.eval(k) );

    }

    index_t rows() const { return _v.rows(); }
    index_t cols() const { return _v.cols(); }
    void setFlag() const { _v.setFlag(); }
    void parse(gsSortedVector<const gsFunctionSet<Scalar>*> & evList) const
    { _v.parse(evList); }

    enum{rowSpan = E2::rowSpan, colSpan = E2::colSpan};

    const gsFeSpace<Scalar> & rowVar() const { return _v.rowVar(); }
    const gsFeSpace<Scalar> & colVar() const { return _v.colVar(); }

    void print(std::ostream &os) const { os << _c <<"*";_v.print(os); }
};


template <typename E1, typename E2>
class collapse_expr : public _expr<collapse_expr<E1, E2> >
{
    typename E1::Nested_t _u;
    typename E2::Nested_t _v;

public:
    enum {ScalarValued = 0, ColBlocks = 0};
    enum { Space = E1::Space || E2::Space};

    typedef typename E1::Scalar Scalar;

    mutable gsMatrix<Scalar> res;

    collapse_expr(_expr<E1> const& u,
              _expr<E2> const& v)
    : _u(u), _v(v) { }

    //EIGEN_STRONG_INLINE MatExprType
    const gsMatrix<Scalar> &
    eval(const index_t k) const
    {
        const index_t nb = rows();
        const MatExprType tmpA = _u.eval(k);
        const MatExprType tmpB = _v.eval(k);

        if (E1::ColBlocks)
        {
            const index_t ur = _v.rows();
            res.resize(nb, ur);
            for (index_t i = 0; i!=nb; ++i)
            {
                res.row(i).transpose().noalias() = tmpA.middleCols(i*ur,ur) * tmpB;
            }
        }
        else if (E2::ColBlocks)
        {
            const index_t ur = _u.cols();
            res.resize(nb, ur);
            for (index_t i = 0; i!=nb; ++i)
            {
                res.row(i).noalias() = tmpA * tmpB.middleCols(i*ur,ur);
            }
        }

        return res;
    }

    index_t rows() const { return E1::ColBlocks ? _u.cols() / _v.rows() : _v.cols() / _u.cols() ; }
    index_t cols() const { return E1::ColBlocks ? _v.rows()  : _u.cols(); }
    void setFlag() const { _u.setFlag(); _v.setFlag(); }
    void parse(gsSortedVector<const gsFunctionSet<Scalar>*> & evList) const
    { _u.parse(evList); _v.parse(evList); }

    enum{rowSpan = E1::Space ? E1::rowSpan : E2::rowSpan, colSpan = E2::Space ? E2::colSpan : E1::colSpan};

    const gsFeSpace<Scalar> & rowVar() const
    { return E1::ColBlocks ? _u.rowVar() : _v.rowVar(); }
    const gsFeSpace<Scalar> & colVar() const
    {
        GISMO_ERROR("none");
    }

    void print(std::ostream &os) const { _u.print(os); os<<"~"; _v.print(os); }
};

// Multi-matrix collapsed by a vector
template <typename E1, typename E2> //EIGEN_STRONG_INLINE
//collapse_expr<E1,E2> const  operator&(<E1> const& u, _expr<E2> const& v)
collapse_expr<E1,E2> collapse( _expr<E1> const& u, _expr<E2> const& v)
{ return collapse_expr<E1, E2>(u, v); }


/*
   Expression for the Frobenius matrix (or double dot) product (first
   version) Also block-wise

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
        //todo: add check() functions, which will evaluate expressions on an empty matrix (no points) to setup initial dimensions ???
        //GISMO_ASSERT(_u.rows() == _v.rows(),
        //             "Wrong dimensions "<<_u.rows()<<"!="<<_v.rows()<<" in % operation");
        //GISMO_ASSERT(_u.cols() == _v.cols(),
        //             "Wrong dimensions "<<_u.cols()<<"!="<<_v.cols()<<" in % operation");
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

    enum{rowSpan = E1::rowSpan, colSpan = E2::rowSpan};

    const gsFeSpace<Scalar> & rowVar() const { return _u.rowVar(); }
    const gsFeSpace<Scalar> & colVar() const { return _v.rowVar(); }

    void print(std::ostream &os) const
    { os << "("; _u.print(os); os<<" % "; _v.print(os); os<<")";}
};

/*

   Expression for the Frobenius matrix (or double dot) product (second
   version), When left hand only side is block-wise

   [A1 A2 A3] : B = [A1:B  A2:B  A3:B]
*/
template <typename E1, typename E2>
class frprod_expr<E1,E2,false> : public _expr<frprod_expr<E1, E2,false> >
{
public:
    typedef typename E1::Scalar Scalar;
    enum {ScalarValued = 0, Space = E1::Space};

private:
    typename E1::Nested_t _u;
    typename E2::Nested_t _v;

    mutable gsMatrix<Scalar> res;

public:

    frprod_expr(_expr<E1> const& u, _expr<E2> const& v)
    : _u(u), _v(v)
    {
        // gsInfo << "expression is space ? "<<E1::Space <<"\n"; _u.print(gsInfo);
        // GISMO_ASSERT(_u.rows() == _v.rows(),
        //              "Wrong dimensions "<<_u.rows()<<"!="<<_v.rows()<<" in % operation");
        // GISMO_ASSERT(_u.cols() == _v.cols(),
        //              "Wrong dimensions "<<_u.cols()<<"!="<<_v.cols()<<" in % operation");
    }

    const gsMatrix<Scalar> & eval(const index_t k) const //todo: specialize for nb==1
    {
        // assert _u.size()==_v.size()
        MatExprType A = _u.eval(k);
        MatExprType B = _v.eval(k);
        const index_t rb = A.rows(); //==cb
        const index_t nb = A.cols() / rb;
        res.resize(nb, 1);
        for (index_t i = 0; i!=nb; ++i) // all with all
                res(i,0) =
                (A.middleCols(i*rb,rb).array() * B.array()).sum();
        return res;
    }

    index_t rows() const { return _u.cols() / _u.rows(); }
    index_t cols() const { return 1; }
    void setFlag() const { _u.setFlag(); _v.setFlag(); }
    void parse(gsSortedVector<const gsFunctionSet<Scalar>*> & evList) const
    { _u.parse(evList); _v.parse(evList); }

    enum{rowSpan = E1::rowSpan, colSpan = E2::rowSpan};

    const gsFeSpace<Scalar> & rowVar() const { return _u.rowVar(); }
    const gsFeSpace<Scalar> & colVar() const { return _v.rowVar(); }

    void print(std::ostream &os) const
    { os << "("; _u.print(os); os<<" % "; _v.print(os); os<<")";}
};

/*
   Expression for scalar division operation (first version)
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
    enum {Space = E1::Space}; // The denominator E2 has to be scalar.

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

    enum{rowSpan = E1::rowSpan, colSpan = E1::colSpan};

    const gsFeSpace<Scalar> & rowVar() const { return _u.rowVar(); }
    const gsFeSpace<Scalar> & colVar() const { return _u.colVar(); }

    void print(std::ostream &os) const
    { os << "("; _u.print(os);os <<" / ";_v.print(os);os << ")"; }
};

/*
   Division specialization (second version) for constant value denominator
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

    enum{rowSpan = E1::rowSpan, colSpan = E1::colSpan};

    const gsFeSpace<Scalar> & rowVar() const { return _u.rowVar(); }
    const gsFeSpace<Scalar> & colVar() const { return _u.colVar(); }

    void print(std::ostream &os) const
    { os << "("; _u.print(os);os <<"/"<< _c << ")"; }
};

/*
   Division specialization (third version) for constant value
   numerator
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

    enum{rowSpan = 0, colSpan = 0};

    const gsFeSpace<Scalar> & rowVar() const { return false; }
    const gsFeSpace<Scalar> & colVar() const { return false; }

    void print(std::ostream &os) const
    { os << "("<< _c <<"/";_u.print(os);os << ")";}
};

/*
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
    enum {Space = E1::Space}; // == E2::Space

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

    enum{rowSpan = E1::rowSpan, colSpan = E1::colSpan};

    index_t cardinality_impl() const { return _u.cardinality_impl(); }

    const gsFeSpace<Scalar> & rowVar() const { return _u.rowVar(); }
    const gsFeSpace<Scalar> & colVar() const { return _v.colVar(); }

    void print(std::ostream &os) const
    { os << "("; _u.print(os);os <<" + ";_v.print(os);os << ")"; }
};

/*// testing, |, ^, &, <<, >>, ||, &&,  unary ~
template <typename E1, typename E2> add_expr<E1,E2> const
operator|(_expr<E1> const& u, _expr<E2> const& v)
{ return add_expr<E1, E2>(u, v); }
template <typename E1, typename E2> add_expr<E1,E2> const
operator^(_expr<E1> const& u, _expr<E2> const& v)
{ return add_expr<E1, E2>(u, v); }
*/


/*
lincom_expr (lc) ?
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

    enum {ScalarValued = 0, ColBlocks = E1::rowSpan};
    enum {Space = E1::Space };

    summ_expr(E1 const& u, E2 const& M) : _u(u), _M(M) { }

    const gsMatrix<Scalar> & eval(const index_t k) const
    {
        MatExprType sl = _u.eval(k);
        const index_t sr  = sl.rows();
        MatExprType ml = _M.eval(k);
        const index_t mr  = ml.rows();
        const index_t mb  = ml.cols() / mr;

        res.setZero(mr, sr * mr);
        for (index_t i = 0; i!=sr; ++i)
            for (index_t t = 0; t!=mb; ++t) // lc
                res.middleCols(i*mr,mr) += sl(i,t) * ml.middleCols(t*mr,mr);
        return res;
    }

    index_t rows() const { return _M.rows(); }
    index_t cols() const { return _u.rows() * _M.rows(); }
    void setFlag() const { _u.setFlag(); _M.setFlag(); }
    void parse(gsSortedVector<const gsFunctionSet<Scalar>*> & evList) const
    { _u.parse(evList); _M.parse(evList); }

    enum{rowSpan = E1::rowSpan, colSpan = 0};

    const gsFeSpace<Scalar> & rowVar() const { return _u.rowVar(); }
    const gsFeSpace<Scalar> & colVar() const { return gsNullExpr<Scalar>::get(); }

    void print(std::ostream &os) const
    { os << "sum("; _M.print(os); os<<","; _u.print(os); os<<")"; }

private:
    typename E1::Nested_t _u;
    typename E2::Nested_t _M;

    mutable gsMatrix<Scalar> res;
};


/*
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
    enum {Space = E1::Space}; // == E2::Space

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
                     "Wrong dimensions "<<_u.rows()<<"!="<<_v.rows()<<" in - operation:\n" << _u <<" minus \n" << _v );
        GISMO_ASSERT(_u.cols() == _v.cols(),
                     "Wrong dimensions "<<_u.cols()<<"!="<<_v.cols()<<" in - operation:\n" << _u <<" minus \n" << _v );
        //gsDebugVar( (_u.eval(k) - _v.eval(k)) );
        //gsDebugVar( (_u.eval(k) - _v.eval(k)).squaredNorm() );
        return (_u.eval(k) - _v.eval(k) );
    }

    index_t rows() const { return _u.rows(); }
    index_t cols() const { return _u.cols(); }
    void setFlag() const { _u.setFlag(); _v.setFlag(); }
    void parse(gsSortedVector<const gsFunctionSet<Scalar>*> & evList) const
    { _u.parse(evList); _v.parse(evList); }

    enum{rowSpan = E1::rowSpan, colSpan = E2::colSpan};

    const gsFeSpace<Scalar> & rowVar() const { return _u.rowVar(); }
    const gsFeSpace<Scalar> & colVar() const { return _v.colVar(); }

    void print(std::ostream &os) const
    { os << "("; _u.print(os); os<<" - ";_v.print(os); os << ")";}
};

//template <typename E2> EIGEN_STRONG_INLINE
//sub_expr<E1,E2> const operator-(typename E2::Scalar const& u, _expr<E2> const& v)
//{ return sub_expr<E1, E2>(u, v); }

/*
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

    enum{rowSpan = E::rowSpan, colSpan = E::rowSpan};

    const gsFeSpace<Scalar> & rowVar() const { return _u.rowVar(); }
    const gsFeSpace<Scalar> & colVar() const { return _u.rowVar(); }

    void print(std::ostream &os) const { os << "symm("; _u.print(os); os <<")"; }
};

template <typename E>
class symmetrize_expr : public _expr<symmetrize_expr<E> >
{
    typename E::Nested_t _u;

    mutable gsMatrix<typename E::Scalar> tmp;
public:
    enum {Space = E::Space};
    typedef typename E::Scalar Scalar;

    symmetrize_expr(_expr<E> const& u)
    : _u(u) { }

    MatExprType eval(const index_t k) const
    {
        //const MatExprType tmp = _u.eval(k);
        tmp = _u.eval(k);
        // todo: avoid temporary or not ?
        return tmp + tmp.transpose();
    }

    index_t rows() const { return _u.rows(); }
    index_t cols() const { return _u.rows(); }
    void setFlag() const { _u.setFlag(); }
    void parse(gsSortedVector<const gsFunctionSet<Scalar>*> & evList) const
    { _u.parse(evList); }

    enum{rowSpan = E::rowSpan, colSpan = E::rowSpan};

    const gsFeSpace<Scalar> & rowVar() const { return _u.rowVar(); }
    const gsFeSpace<Scalar> & colVar() const { return _u.rowVar(); }
    index_t cardinality_impl() const { return _u.cardinality_impl(); }

    void print(std::ostream &os) const { os << "symmetrize("; _u.print(os); os <<")"; }
};

/* Symmetrization operation
template <typename E> symm_expr<E> const
symm(_expr<E> const& u) { return symm_expr<E>(u);}
*/

#undef MatExprType
#undef AutoReturn_t
//----------------------------------------------------------------------------------

// Returns the unit as an expression
//EIGEN_STRONG_INLINE _expr<real_t> one() { return _expr<real_t>(1); }

/// The identity matrix of dimension \a dim
EIGEN_STRONG_INLINE idMat_expr id(const index_t dim) { return idMat_expr(dim); }

/// The gradient of a solution variable
template<class T> EIGEN_STRONG_INLINE
solGrad_expr<T> grad(const gsFeSolution<T> & u) { return solGrad_expr<T>(u); }

/// The gradient of a variable
template<class E> EIGEN_STRONG_INLINE
grad_expr<E> grad(const E & u) { return grad_expr<E>(u); }

/// The derivative of the jacobian of a geometry map with respect to a coordinate.
template<class E> EIGEN_STRONG_INLINE
dJacdc_expr<E> dJacdc(const E & u, index_t c) { return dJacdc_expr<E>(u,c); }

/// The curl of a finite element variable
template<class T> EIGEN_STRONG_INLINE
curl_expr<T> curl(const gsFeVariable<T> & u) { return curl_expr<T>(u); }

/// The nabla (\f$\nabla\f$) of a finite element variable
template<class T> EIGEN_STRONG_INLINE
nabla_expr<T> nabla(const gsFeVariable<T> & u) { return nabla_expr<T>(u); }

/// The (outer pointing) boundary normal of a geometry map
template<class T> EIGEN_STRONG_INLINE
onormal_expr<T> nv(const gsGeometryMap<T> & u) { return onormal_expr<T>(u); }

template<class T> EIGEN_STRONG_INLINE
normal_expr<T> sn(const gsGeometryMap<T> & u) { return normal_expr<T>(u); }

/// The tangent boundary vector of a geometry map in 2D
template<class T> EIGEN_STRONG_INLINE
tangent_expr<T> tv(const gsGeometryMap<T> & u) { return tangent_expr<T>(u); }

/// The laplacian of a finite element variable
template<class T> EIGEN_STRONG_INLINE
lapl_expr<T> lapl(const gsFeVariable<T> & u) { return lapl_expr<T>(u); }

template<class T> EIGEN_STRONG_INLINE
solLapl_expr<T> slapl(const gsFeSolution<T> & u) { return solLapl_expr<T>(u); }

/// The hessian of a solution variable
template<class T> EIGEN_STRONG_INLINE
solHess_expr<T> shess(const gsFeSolution<T> & u) { return solHess_expr<T>(u); }

/// The first fundamental form of \a G //fform is buggy ?
// template<class T> EIGEN_STRONG_INLINE fform_expr<T> fform(const gsGeometryMap<T> & G) { return fform_expr<T>(G); }

/// The Jacobian matrix of a geometry map
template<class T> EIGEN_STRONG_INLINE
jacG_expr<T> jac(const gsGeometryMap<T> & G) { return jacG_expr<T>(G); }

/// The Jacobian matrix of a FE variable
template<class E> EIGEN_STRONG_INLINE
jac_expr<E> jac(const E & u) { return jac_expr<E>(u); }

/// The Jacobian matrix of a vector function
template<class T> EIGEN_STRONG_INLINE
fjac_expr<T> fjac(const gsFeVariable<T> & u) { return fjac_expr<T>(u); }

/// The Hessian matrix of (each coordianate of) the geometry map \a G
//template<class T> EIGEN_STRONG_INLINE hess_expr<T> hess(const gsGeometryMap<T> & G) { return hess_expr<T >(G); }

/// The Hessian matrix of a finite element variable
//template<class T> EIGEN_STRONG_INLINE hess_expr<T> hess(const gsFeVariable<T> & u) { return hess_expr<T>(u); }

template<class E> EIGEN_STRONG_INLINE
hess_expr<E> hess(const E & u) { return hess_expr<E>(u); }

/// The partial derivatives of the Jacobian matrix of a geometry map
template<class T> EIGEN_STRONG_INLINE
dJacG_expr<T> dJac(const gsGeometryMap<T> & G) { return dJacG_expr<T>(G); }

/// The measure of a geometry map
template<class T> EIGEN_STRONG_INLINE
meas_expr<T> meas(const gsGeometryMap<T> & G) { return meas_expr<T>(G); }

/// Multiplication operator for expressions
template <typename E1, typename E2> EIGEN_STRONG_INLINE
mult_expr<E1,E2> const operator*(_expr<E1> const& u, _expr<E2> const& v)
{ return mult_expr<E1, E2>(u, v); }

template <typename E2> EIGEN_STRONG_INLINE
mult_expr<typename E2::Scalar,E2,false> const
operator*(typename E2::Scalar const& u, _expr<E2> const& v)
{ return mult_expr<typename E2::Scalar, E2, false>(u, v); }
//{ return mult_expr<_expr<typename E2::Scalar>,E2, false>(u, v); }

template <typename E1> EIGEN_STRONG_INLINE
mult_expr<typename E1::Scalar,E1,false> const
operator*(_expr<E1> const& v, typename E1::Scalar const& u)
{ return mult_expr<typename E1::Scalar,E1, false>(u, v); }
//{ return mult_expr<_expr<typename E1::Scalar>,E1, false>(u, v); }

template <typename E1> EIGEN_STRONG_INLINE
mult_expr<typename E1::Scalar,E1,false> const
operator-(_expr<E1> const& u)
{ return mult_expr<typename E1::Scalar,E1, false>(-1, u); }
//{ return mult_expr<_expr<typename E1::Scalar>,E1, false>(-1, u); }

/*
template <typename E1> mult_expr<gsMatrix<typename E1::Scalar>,E1,false> const
operator*(gsMatrix<typename E1::Scalar> const& u, _expr<E1> const& v)
{ return mult_expr<gsMatrix<typename E1::Scalar>,E1, false>(u, v); }
*/

/// Frobenious product (also known as double dot product) operator for expressions
template <typename E1, typename E2> EIGEN_STRONG_INLINE
frprod_expr<E1,E2> const  operator%(_expr<E1> const& u, _expr<E2> const& v)
{ return frprod_expr<E1, E2>(u, v); }

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

/// Addition operator for expressions
template <typename E1, typename E2> EIGEN_STRONG_INLINE
add_expr<E1,E2> const operator+(_expr<E1> const& u, _expr<E2> const& v)
{ return add_expr<E1, E2>(u, v); }

/// Matrix-summation operator for expressions
template <typename E1, typename E2> EIGEN_STRONG_INLINE
summ_expr<E1,E2> const summ(E1 const & u, E2 const& M)
{ return summ_expr<E1,E2>(u, M); }

/// Matrix by space TODO: find better name and/or description? And is this the best place?
/// [Jg Jg Jg] * Jb ..
template <typename E1, typename E2> EIGEN_STRONG_INLINE
matrix_by_space_expr<E1,E2> const matrix_by_space(E1 const & u, E2 const& v)
{ return matrix_by_space_expr<E1,E2>(u, v); }

/// Matrix by space TODO: find better name and/or description? And is this the best place?
/// [Jg Jg Jg] * Jb ..
template <typename E1, typename E2> EIGEN_STRONG_INLINE
matrix_by_space_expr_tr<E1,E2> const matrix_by_space_tr(E1 const & u, E2 const& v)
{ return matrix_by_space_expr_tr<E1,E2>(u, v); }

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


//----------------------------------------------------------------------------------
#if __cplusplus >= 201402L || _MSVC_LANG >= 201402L

// Shortcuts for common quantities, for instance function
// transformations by the geometry map \a G
#define GISMO_SHORTCUT_VAR_EXPRESSION(name,impl) template<class E> EIGEN_STRONG_INLINE \
auto name(const E & u) { return impl; }
#define GISMO_SHORTCUT_MAP_EXPRESSION(name,impl) template<class T> EIGEN_STRONG_INLINE \
auto name(const gsGeometryMap<T> & G) { return impl; }
#define GISMO_SHORTCUT_PHY_EXPRESSION(name,impl) template<class E> EIGEN_STRONG_INLINE \
    auto name(const E & u, const gsGeometryMap<typename E::Scalar> & G) { return impl; }

// Divergence
GISMO_SHORTCUT_VAR_EXPRESSION(  div, jac(u).trace() )
GISMO_SHORTCUT_PHY_EXPRESSION( idiv, ijac(u,G).trace()    )

// The unit (normalized) boundary (outer pointing) normal
GISMO_SHORTCUT_MAP_EXPRESSION(unv, nv(G).normalized()   ) //(!) bug + mem. leak

GISMO_SHORTCUT_PHY_EXPRESSION(igrad, grad(u)*jac(G).ginv() ) // transpose() problem ??
GISMO_SHORTCUT_VAR_EXPRESSION(igrad, grad(u) ) // u is presumed to be defined over G

GISMO_SHORTCUT_PHY_EXPRESSION( ijac, jac(u) * jac(G).ginv())

GISMO_SHORTCUT_PHY_EXPRESSION(ihess,
jac(G).ginv().tr()*( hess(u) - summ(igrad(u,G),hess(G)) ) * jac(G).ginv() )
GISMO_SHORTCUT_VAR_EXPRESSION(ihess, hess(u) )

GISMO_SHORTCUT_PHY_EXPRESSION(ilapl, ihess(u,G).trace()   )
GISMO_SHORTCUT_VAR_EXPRESSION(ilapl, hess(u).trace() )

GISMO_SHORTCUT_VAR_EXPRESSION(fform, jac(u).tr()*jac(u) )

#else
// note for C++11 and older: to fix
//error: could not convert .. from ‘A' to 'B'
// we need to update the return type from B to A

#define GISMO_SHORTCUT_VAR_EXPRESSION(name,impl) \
name(const E & u) { return impl; }
#define GISMO_SHORTCUT_MAP_EXPRESSION(name,impl) \
name(const gsGeometryMap<T> & G) { return impl; }

//#define GISMO_SHORTCUT_PHY_EXPRESSION(name,impl) name(const gsFeVariable<T> & u, const gsGeometryMap<T> & G) { return impl; }
#define GISMO_SHORTCUT_PHY_EXPRESSION(name,impl) \
name(const E & u, const gsGeometryMap<typename E::Scalar> & G) { return impl; }

// Divergence
template<class E> EIGEN_STRONG_INLINE trace_expr<jac_expr<E> >
GISMO_SHORTCUT_VAR_EXPRESSION(div, jac(u).trace() )

template<class T> EIGEN_STRONG_INLINE normalized_expr<onormal_expr<T> >
GISMO_SHORTCUT_MAP_EXPRESSION(unv, nv(G).normalized() )

//template<class T> EIGEN_STRONG_INLINE
template<class E> EIGEN_STRONG_INLINE
mult_expr<grad_expr<E>,jacGinv_expr<typename E::Scalar>, 0>
GISMO_SHORTCUT_PHY_EXPRESSION(igrad, grad(u)*jac(G).ginv())

template<class E> EIGEN_STRONG_INLINE grad_expr<E> // u is presumed to be defined over G
GISMO_SHORTCUT_VAR_EXPRESSION(igrad, grad(u))

template<class E> EIGEN_STRONG_INLINE
mult_expr<jac_expr<E>,jacGinv_expr<typename E::Scalar>, 1>
GISMO_SHORTCUT_PHY_EXPRESSION(ijac, jac(u) * jac(G).ginv() )

template<class E> EIGEN_STRONG_INLINE
trace_expr<mult_expr<jac_expr<E>,jacGinv_expr<typename E::Scalar>, 1> >
GISMO_SHORTCUT_PHY_EXPRESSION(idiv, ijac(u,G).trace() )

template<class E> EIGEN_STRONG_INLINE
mult_expr<mult_expr<tr_expr<jacGinv_expr<typename E::Scalar> >,sub_expr<hess_expr<E>,summ_expr<mult_expr<grad_expr<E>, jacGinv_expr<typename E::Scalar>, 0>, hess_expr<gsGeometryMap<typename E::Scalar> > > >, 0>, jacGinv_expr<typename E::Scalar>, 0>
GISMO_SHORTCUT_PHY_EXPRESSION(ihess, jac(G).ginv().tr()*(hess(u)-summ(igrad(u,G),hess(G)))*jac(G).ginv() )

template<class E> EIGEN_STRONG_INLINE trace_expr<
mult_expr<mult_expr<tr_expr<jacGinv_expr<typename E::Scalar> >,sub_expr<hess_expr<E>,summ_expr<mult_expr<grad_expr<E>, jacGinv_expr<typename E::Scalar>, 0>, hess_expr<gsGeometryMap<typename E::Scalar> > > >, 0>, jacGinv_expr<typename E::Scalar>, 0>
    >
GISMO_SHORTCUT_PHY_EXPRESSION(ilapl, ihess(u,G).trace() )

template<class E> EIGEN_STRONG_INLINE hess_expr<E>
GISMO_SHORTCUT_VAR_EXPRESSION(ihess, hess(u) )

template<class E> EIGEN_STRONG_INLINE trace_expr<hess_expr<E> >
GISMO_SHORTCUT_VAR_EXPRESSION(ilapl, hess(u).trace() )

#endif
#undef GISMO_SHORTCUT_PHY_EXPRESSION
#undef GISMO_SHORTCUT_VAR_EXPRESSION
#undef GISMO_SHORTCUT_MAP_EXPRESSION

} // namespace expr

} //namespace gismo
