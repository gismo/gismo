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
#include <gsAssembler/gsDirichletValues.h>
#include <gsMSplines/gsMappedBasis.h>


namespace gismo
{

// Adaptor to compute Hessian
template <typename Derived>
void secDerToHessian(const Eigen::DenseBase<Derived> &  secDers,
                     const index_t dim,
                     gsMatrix<typename Derived::Scalar> & hessian)
{
    const index_t sz = dim*(dim+1)/2;
    const gsAsConstMatrix<typename Derived::Scalar>
        ders(secDers.derived().data(), sz, secDers.size() / sz );
    hessian.resize(dim*dim, ders.cols() );

    switch ( dim )
    {
    case 1:
        hessian = secDers.transpose(); //==ders
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

/// Struct containing information for matrix assembly
template<class T> struct gsFeSpaceData
{
    gsFeSpaceData(const gsFunctionSet<T> & _fs, index_t _dim, index_t _id):
    fs(&_fs), dim(give(_dim)), id(give(_id)) { }

    const gsFunctionSet<T> * fs;
    index_t dim, id;
    gsDofMapper mapper;
    gsMatrix<T> fixedDofs;
    index_t cont; //int. coupling

    bool valid() const
    {
        GISMO_ASSERT(nullptr!=fs, "Invalid pointer.");
        return static_cast<size_t>(fs->size()*dim)==mapper.mapSize();
    }

    void init()
    {
        GISMO_ASSERT(nullptr!=fs, "Invalid pointer.");
        if (const gsMultiBasis<T> * mb =
            dynamic_cast<const gsMultiBasis<T>*>(fs) )
            mapper = gsDofMapper(*mb, dim );
        else if (const gsBasis<T> * b =
                 dynamic_cast<const gsBasis<T>*>(fs) )
            mapper = gsDofMapper(*b, dim );
        mapper.finalize();
        fixedDofs.clear();
        cont = -1;
    }
};

// Forward declaration in gismo namespace
template<class T> class gsExprHelper;

/** @namespace gismo::expr

    @brief
    This namespace contains expressions used for FE computations

    \ingroup Assembler
*/
namespace expr
{

template <class E> struct is_arithmetic{enum{value=0};};
template <> struct is_arithmetic<real_t>{enum{value=1};};
template <typename E, bool = is_arithmetic<E>::value >
class _expr {using E::GISMO_ERROR_expr;};

template<class T> class gsFeSpace;
template<class T> class gsFeVariable;
template<class T> class gsFeSolution;
template<class E> class symm_expr;
template<class E> class symmetrize_expr;
template<class E> class normalized_expr;
template<class E> class trace_expr;
template<class E> class integral_expr;
template<class E> class adjugate_expr;
template<class E> class norm_expr;
template<class E> class sqNorm_expr;
template<class E> class det_expr;
template<class E> class value_expr;
template<class E> class asdiag_expr;
template<class E> class max_expr;
template<class E> class rowsum_expr;
template<class E> class colsum_expr;
template<class E> class col_expr;
template<class T> class meas_expr;
template<class E> class inv_expr;
template<class E, bool cw = false> class tr_expr;
template<class E> class cb_expr;
template<class E> class abs_expr;
template<class E> class pow_expr;
template<class E> class sign_expr;
template<class E> class ppart_expr;
template<class T> class cdiam_expr;
template<class E> class temp_expr;
template<class E1, class E2, bool = E1::ColBlocks && !E1::ScalarValued && !E2::ScalarValued> class mult_expr
{using E1::GISMO_ERROR_mult_expr_has_invalid_template_arguments;};

// Call as pow(a,b)
template<class E> pow_expr<E>
pow(_expr<E> const& u, real_t q) { return pow_expr<E>(u,q); }

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

#  define Temporary_t typename util::conditional<ScalarValued,Scalar,   \
        typename gsMatrix<Scalar>::Base >::type
#if __cplusplus >= 201402L || _MSVC_LANG >= 201402L // c++14
#  define MatExprType  auto
#  define AutoReturn_t auto
//note: in c++11 auto-return requires -> decltype(.)
#else // 199711L, 201103L
#  define MatExprType typename gsMatrix<real_t>::constRef
#  define AutoReturn_t typename util::conditional<ScalarValued,real_t,MatExprType>::type
#endif

/**
   \brief Base class for all expressions
*/
template <typename E>
class _expr<E, false>
{
protected://private:
    _expr(){}
    _expr(const _expr&) { }
public:
    // Defined in derived classes: enum { Space, ScalarValued, ColBlocks }
    // - ScalarValued: 0 is a scalar (must have Space=0),1 one denotes gsMatrix
    // - ColBlocks: the expression stacks matrices per basis function
    // - Space: 0: not a trial nor a test object (eg. normal vector, force function)
    //          1: a test object  (essentially a right-hand side vector expression)
    //          2: a trial object
    //          3: a trial+trial object (essentially a matrix expression)

    typedef typename expr_traits<E>::Nested_t Nested_t;
    typedef typename expr_traits<E>::Scalar   Scalar;

    /// Prints the expression as a string to \a os
    void print(std::ostream &os) const
    {
        //gsInfo<<"\n Space="<<E::Space<<", ScV="<<E::ScalarValued<<", ColBlocks="<<E::ColBlocks<<"\n";
        static_cast<E const&>(*this).print(os);
        os<<"\n";
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
        os << (isVectorTr() ? "VectorTr " :
               (isVector() ? "Vector " :
                (isMatrix() ? "Matrix " :
                 "Scalar ") ) )
           <<"expression of size "<< rows() // bug: this might be invalid if unparsed
           << " x "<<cols()<<"\n";
        print(os);
        return os;
    }

    /// Evaluates the expression at evaluation point indexed by \a k
    MatExprType eval(const index_t k) const
    { return static_cast<E const&>(*this).eval(k); }

    /// Returns the transpose of the expression
    tr_expr<E> tr() const
    { return tr_expr<E,false>(static_cast<E const&>(*this)); }

    /// Returns the coordinate-wise transpose of the expression
    tr_expr<E,true> cwisetr() const
    { return tr_expr<E,true>(static_cast<E const&>(*this)); }

    /// Returns the puts the expression to colBlocks
    cb_expr<E> cb() const
    { return cb_expr<E>(static_cast<E const&>(*this)); }

    /// Returns the sign of the expression
    sign_expr<E> sgn(Scalar tolerance=0) const
    { return sign_expr<E>(static_cast<E const&>(*this), tolerance); }

    /// Returns the expression's positive part
    ppart_expr<E> ppart() const
    { return ppart_expr<E>(static_cast<E const&>(*this)); }

    /// Returns the expression's negative part
    mult_expr<real_t, ppart_expr<mult_expr<double,E,false>> , false> 
    npart() const { return -1* ( -(*this) ).ppart() ; }

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

    /// Returns the vector normalized to unit length
    normalized_expr<E> normalized() const
    { return normalized_expr<E>(static_cast<E const&>(*this)); }

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
    max_expr<E> max() const
    { return max_expr<E>(static_cast<E const&>(*this)); }

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

    ///\brief Returns true iff the expression is scalar-valued.
    /// \note This is a runtime check, for compile-time check use E::ScalarValued
    bool isScalar() const { return rows()*cols()<=1; } //!rowSpan && !colSpan

    static bool isVector  () { return 1==E::Space; }
    static bool isVectorTr() { return 2==E::Space; }
    static bool isMatrix  () { return 3==E::Space; }

    ///\brief Parse the expression and discover the list of evaluation
    ///sources, also sets the required evaluation flags
    void parse(gsExprHelper<Scalar> & evList) const
    { static_cast<E const&>(*this).parse(evList); }

    template<class op> void apply(op & _op) const
    { static_cast<E const&>(*this).apply(_op); }

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

    E const & derived() const { return static_cast<const E&>(*this); }
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
        static gsFeSpace<T> vv(-1);
        return vv;
    }

    typedef T Scalar;
    gsMatrix<T> eval(const index_t) const { GISMO_ERROR("gsNullExpr"); }
    inline index_t rows() const { GISMO_ERROR("gsNullExpr"); }
    inline index_t cols() const { GISMO_ERROR("gsNullExpr"); }
    void parse(gsExprHelper<T> &) const { }

    const gsFeSpace<T> & rowVar() const { GISMO_ERROR("gsNullExpr"); }
    const gsFeSpace<T> & colVar() const { GISMO_ERROR("gsNullExpr"); }

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
    const gsFunctionSet<Scalar> * m_fs; ///< Evaluation source for this FE variable
    const gsFuncData<Scalar>    * m_fd; ///< Temporary variable storing flags and evaluation data
    index_t m_d;                   ///< Dimension of this (scalar or vector) variable
    bool m_isAcross; ///< true when this expression is evaluated across an interface

public:

    /// Returns whether this expression is evaluated across an interface
    bool isAcross() const { return m_isAcross; }

    E right() const
    {
        E ac(this->derived());
        ac.m_fs = m_fs;//needed?
        ac.m_isAcross = true;
        return ac;
    }

    E left() const
    {
        E ac(this->derived());
        ac.m_fs = m_fs;
        ac.m_isAcross = false;
        return ac;
    }

    /// Returns the function source
    const gsFunctionSet<Scalar> & source() const {return *m_fs;}

    /// Returns the function data
    const gsFuncData<Scalar> & data() const
    {
        GISMO_ASSERT(NULL!=m_fd, "FuncData member not registered "<<this<<"/"<< m_fs);
        return *m_fd;
    }

    // used by FeSpace, FeVariable, ..
    void parse(gsExprHelper<Scalar> & evList) const
    {
        evList.add(*this);
        this->m_fd->flags |= NEED_VALUE | NEED_ACTIVE;
    }

    index_t cardinality_impl() const
    {
        GISMO_ASSERT(this->data().actives.rows()!=0,"Cardinality depends on the NEED_ACTIVE flag");
        return m_d * this->data().actives.rows();
    }

    //public for now due to Bc
    void setSource(const gsFunctionSet<Scalar> & fs) { m_fs = &fs;}

private:
    void setData(const gsFuncData<Scalar> & val) { m_fd = &val;}
    void setDim(index_t _d) { m_d = _d; }
    void clear() { m_fs = NULL; }

protected:
    explicit symbol_expr(index_t _d)
    : m_fs(NULL), m_fd(NULL), m_d(_d), m_isAcross(false) { }

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

    index_t rows() const
    {
        GISMO_ASSERT(NULL!=m_fs, "FeVariable: Function source not registered");
        return m_fs->targetDim();
    }

    index_t cols() const { return m_d; }

    void print(std::ostream &os) const { os << "u"; }

public:

    /// Returns the vector dimension of the FE variable
    index_t dim() const { return m_d;}

    /// Returns the target dimension of the FE variable
    /// before vector-replication
    index_t targetDim() const { return m_fs->targetDim(); }

    /// Returns the parameter domain dimension the FE variable
    index_t parDim() const { return m_fs->domainDim(); }

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

    enum { Space = E::Space, ScalarValued = 0, ColBlocks = 0 };

    col_expr(const E & c, const index_t i) : _c(c), _i(i) { }

public:

    //ConstColXpr
    inline MatExprType eval(const index_t k) const { return _c.eval(k).col(_i); }

    index_t rows() const { return _c.rows(); }
    index_t cols() const { return 1; }
    void parse(gsExprHelper<Scalar> & evList) const { _c.parse(evList); }

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

    _expr(Scalar c) : _c(give(c)) { }

public:
    enum {Space = 0, ScalarValued = 1, ColBlocks= 0};

    inline Scalar eval(const index_t ) const { return _c; }

    inline _expr val() const { return *this; }
    index_t rows() const { return 0; }
    index_t cols() const { return 0; }
    void parse(gsExprHelper<Scalar> &) const { }

    const gsFeSpace<T> & rowVar() const { return gsNullExpr<T>::get(); }
    const gsFeSpace<T> & colVar() const { return gsNullExpr<T>::get(); }

    void print(std::ostream &os) const { os<<_c; }
};

/*
  Geometry map expression
*/
template<class T>
class gsGeometryMap : public _expr<gsGeometryMap<T> >
{
    const gsFunctionSet<T> * m_fs; ///< Evaluation source for this geometry map
    const gsMapData<T>     * m_fd; ///< Temporary variable storing flags and evaluation data
    //index_t d, n;

    bool m_isAcross; ///< true when the patch evaluated is across an interface

public:
    enum {Space = 0, ScalarValued= 0, ColBlocks= 0};

    bool isAcross() const { return m_isAcross; }

    gsGeometryMap right() const
    {
        gsGeometryMap ac;
        ac.m_fs = m_fs;
        ac.m_isAcross = true;
        return ac;
    }

    gsGeometryMap left() const
    {
        gsGeometryMap ac;
        ac.m_fs = m_fs;
        ac.m_isAcross = false;
        return ac;
    }

    /// Returns the function source
    const gsFunctionSet<T> & source() const {return *m_fs;}

    /// Returns the function data
    const gsMapData<T> & data() const
    {
        GISMO_ASSERT(NULL!=m_fd, "gsGeometryMap: invalid data "<< m_fs <<","<<m_fd);
        return *m_fd;
    }

    index_t targetDim() const { return m_fs->targetDim();}

    void deformBy( const gsFeSolution<T> & deformation) const
    {
        const gsMatrix<T> &defVector = deformation.coefs();
        const index_t dim = m_fs->domainDim();

        const gsMultiBasis<T> & mb = static_cast<const gsMultiBasis<T>&>(deformation.space().source());
        const gsMultiPatch<T> & mp = static_cast<const gsMultiPatch<T>&>(*this->m_fs );
        GISMO_ASSERT( dynamic_cast<const gsMultiBasis<T>*>(&deformation.space().source()), "error");
        GISMO_ASSERT( dynamic_cast<const gsMultiPatch<T>*>( this->m_fs), "error");

        // For every patch of the MultiPatch
        for ( index_t p=0; p < mp.nPatches(); p++ )
        {
            // Get the patch's coefficients
            gsMatrix<T> &result = mp.patch(p).coefs();

            // Number of basis functions of patch with index p
            const index_t sz  = mb[p].size();

            // For all components
            for (index_t c = 0; c!=dim; c++)
            {
                // loop over all basis functions (even the eliminated ones)
                for (index_t i = 0; i < sz; ++i)
                {
                    const int ii = deformation.mapper().index(i, p, c);
                    if ( deformation.mapper().is_free_index(ii) ) // DoF value is in the defVector
                    {
                        result(i,c) += defVector.at(ii);
                    }

                }
            }
        }
    }
public:
    typedef T Scalar;

    friend class gismo::gsExprHelper<Scalar>;

    void print(std::ostream &os) const { os << "G"; }

    auto eval(const index_t k) const -> decltype(m_fd->values[0].col(k))
    { return m_fd->values[0].col(k); }

protected:

    gsGeometryMap() : m_fs(NULL), m_fd(NULL), m_isAcross(false) { }

    void setSource(const gsFunctionSet<Scalar> & fs) { m_fs = &fs;}
    void setData(const gsMapData<Scalar> & val) { m_fd = &val;}

public:

    index_t rows() const { return m_fs->targetDim(); }
    index_t cols() const { return 1; }

    const gsFeSpace<T> & rowVar() const { return gsNullExpr<T>::get(); }
    const gsFeSpace<T> & colVar() const { return gsNullExpr<T>::get(); }

    void parse(gsExprHelper<Scalar> & evList) const
    {
        evList.add(*this);
        m_fd->flags |= NEED_VALUE;
    }
};

// Traits for gsGeometryMap
template <typename T>  struct expr_traits<gsGeometryMap<T> >
{
    typedef T Scalar;
    typedef const gsGeometryMap<T> Nested_t; // nesting without ref!
};

/*
  Expression for the measure of a geometry map
*/
template<class T>
class meas_expr : public _expr<meas_expr<T> >
{
    typename gsGeometryMap<T>::Nested_t _G;

public:
    enum {Space = 0, ScalarValued = 1, ColBlocks = 0};

    typedef T Scalar;

    meas_expr(const gsGeometryMap<T> & G) : _G(G) { }

    T eval(const index_t k) const
    {
        return _G.data().measures.at(k);
    }

    index_t rows() const { return 0; }
    index_t cols() const { return 0; }
    void parse(gsExprHelper<Scalar> & evList) const
    {
        evList.add(_G);
        _G.data().flags |= NEED_MEASURE;
    }

    const gsFeSpace<T> & rowVar() const { return gsNullExpr<T>::get(); }
    const gsFeSpace<T> & colVar() const { return gsNullExpr<T>::get(); }

    void print(std::ostream &os) const { os << "meas("; _G.print(os); os <<")"; }
};


/*
  An element object collecting relevant expressions
*/
template<class T>
class gsFeElement
{
    friend class cdiam_expr<T>;

    const gsDomainIterator<T> * m_di; ///< Pointer to the domain iterator

    const gsVector<T> * m_weights;
    //const gsMatrix<T> * m_points;

    gsFeElement(const gsFeElement &);
public:
    typedef T Scalar;

    gsFeElement() : m_di(NULL), m_weights(nullptr) { }

    void set(const gsDomainIterator<T> & di, const gsVector<T> & weights)
    { m_di = &di, m_weights = &weights; }

    bool isValid() const { return nullptr!=m_weights; }

    const gsVector<T> & weights() const {return *m_weights;}

    template<class E>
    integral_expr<E> integral(const _expr<E>& ff) const
    { return integral_expr<E>(*this,ff); }

    typedef integral_expr<T> AreaRetType;
    AreaRetType area() const
    { return integral(_expr<T,true>(1)); }

    typedef integral_expr<meas_expr<T> > PHAreaRetType;
    /// The diameter of the element on the physical space
    PHAreaRetType area(const gsGeometryMap<Scalar> & _G) const
    { return integral(meas_expr<T>(_G)); }

    typedef pow_expr<integral_expr<T> > DiamRetType;
    /// The diameter of the element (on parameter space)
    DiamRetType diam() const //-> int(1)^(1/d)
    { return pow(integral(_expr<T,true>(1)),(T)(1)/(T)(2)); }

    typedef pow_expr<integral_expr<meas_expr<T> > > PHDiamRetType;
    /// The diameter of the element on the physical space
    PHDiamRetType diam(const gsGeometryMap<Scalar> & _G) const
    { return pow(integral(meas_expr<T>(_G)),(T)(1)/(T)(2)); }

    //const gsMatrix<T> points() const {return pts;}

    //index_t dim() { return di->
    
    void print(std::ostream &os) const { os << "e"; }

    void parse(gsExprHelper<T> & evList) const
    {
        GISMO_ERROR("EL");
        evList.add(*this);
        this->data().flags |= NEED_VALUE;
    }
};

/**
   An expression of the element diameter
*/
template<class E>
class integral_expr : public _expr<integral_expr<E> >
{
public:
    //typedef typename E::Scalar Scalar;
    typedef real_t Scalar;
    mutable Scalar m_val;
private:
    const gsFeElement<Scalar> & _e; ///<Reference to the element
    typename _expr<E>::Nested_t _ff;
public:
    enum {Space= 0, ScalarValued= 1, ColBlocks = 0};

    integral_expr(const gsFeElement<Scalar> & el, const _expr<E> & u)
    : m_val(-1), _e(el), _ff(u) { }

    const Scalar & eval(const index_t k) const
    {
        GISMO_ENSURE(_e.isValid(), "Element is valid within integrals only.");
        // if (0==k)
        {
            const Scalar * w = _e.weights().data();
            m_val = (*w) * _ff.val().eval(0);
            for (index_t j = 1; j != _e.weights().rows(); ++j)
                m_val += (*(++w)) * _ff.val().eval(j);
        }
        return m_val;
    }

    inline const integral_expr<E> & val() const { return *this; }
    inline index_t rows() const { return 0; }
    inline index_t cols() const { return 0; }
    void parse(gsExprHelper<Scalar> & evList) const
    {
        _ff.parse(evList);
    }

    const gsFeSpace<Scalar> & rowVar() const { return gsNullExpr<Scalar>::get(); }
    const gsFeSpace<Scalar> & colVar() const { return gsNullExpr<Scalar>::get(); }

    void print(std::ostream &os) const
    {
        os << "integral(";
        _ff.print(os);
        os <<")";
    }
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
  void parse(gsExprHelper<Scalar> &) const { }
  const gsFeVariable<T> & rowVar() const { gsNullExpr<Scalar>::get(); }
  const gsFeVariable<T> & colVar() const { gsNullExpr<Scalar>::get(); }

  void print(std::ostream &os) const
  { os << "diam(e)"; }
  };
*/

template <typename T>
struct expr_traits<gsFeVariable<T> >
{
    typedef T Scalar;
    typedef const gsFeVariable<T> Nested_t;
};

/**
   Expression for finite element variables or PDE coefficient functionals.
   This can be e.g. a diffusion coefficient, or an isogeometric function.
*/
template<class T>
class gsFeVariable  : public symbol_expr< gsFeVariable<T> >
{
    friend class gismo::gsExprHelper<T>;
    typedef symbol_expr< gsFeVariable<T> > Base;
protected:
    explicit gsFeVariable(index_t _d = 1) : Base(_d) { }
public:
    enum {Space = 0, ScalarValued = 0, ColBlocks = 0};
};

template<class T>
class gsComposition : public symbol_expr< gsComposition<T> >
{ //comp(f,G)
    friend class gismo::gsExprHelper<T>;
    typedef symbol_expr< gsComposition<T> > Base;
    typename gsGeometryMap<T>::Nested_t _G;
protected:
    explicit gsComposition(const gsGeometryMap<T> & G, index_t _d = 1)
    : Base(_d), _G(G) { }
public:
    enum {Space = 0, ScalarValued= 0, ColBlocks= 0};

    typename gsMatrix<T>::constColumn
    eval(const index_t k) const { return this->m_fd->values[0].col(k); }

    const gsGeometryMap<T> & inner() const { return _G;};

    void parse(gsExprHelper<T> & evList) const
    {
        //evList.add(_G); //done in gsExprHelper
        evList.add(*this);
        this->data().flags |= NEED_VALUE|NEED_ACTIVE;
        //_G.data().flags  |= NEED_VALUE; //done in gsExprHelper
    }
};


/**
   Expression for finite element variable in an isogeometric function
   space

   This corresponds to an FE variable that
   belongs to an isogeometric function space
*/
template<class T>
class gsFeSpace :public symbol_expr< gsFeSpace<T> >
{
    friend class gsNullExpr<T>;

protected:
    typedef symbol_expr< gsFeSpace<T> > Base;

    // contains id, mapper, fixedDofs, etc
    gsFeSpaceData<T> * m_sd;

public:
    enum{Space = 1, ScalarValued=0, ColBlocks=0};// test space

    typedef const gsFeSpace Nested_t; //no ref

    typedef T Scalar;

    const gsFeSpace<T> & rowVar() const {return *this;}

    gsDofMapper & mapper()
    {
        GISMO_ASSERT(NULL!=m_sd, "Space/mapper not properly initialized.");
        return m_sd->mapper;
    }

    const gsDofMapper & mapper() const
    {return const_cast<gsFeSpace*>(this)->mapper();}

    inline const gsMatrix<T> & fixedPart() const {return m_sd->fixedDofs;}
    gsMatrix<T> & fixedPart() {return m_sd->fixedDofs;}

    index_t   id() const { return (m_sd ? m_sd->id : -101); }
    void setSpaceData(gsFeSpaceData<T>& sd) {m_sd = &sd;}

    index_t   interfaceCont() const {return m_sd->cont;}
    index_t & setInterfaceCont(const index_t _r) const
    {
        GISMO_ASSERT(_r>-2 && _r<1, "Invalid or not implemented (r="<<_r<<").");
        return m_sd->cont = _r;
    }

    gsFeSolution<T> function(const gsMatrix<T>& solVector) const
    { return gsFeSolution<T>(*this); }

    void getCoeffs(const gsMatrix<T>& solVector, gsMatrix<T> & result,
                   const index_t p = 0) const
    {
        const index_t dim = this->dim();

        const gsMultiBasis<T> & mb = static_cast<const gsMultiBasis<T>&>(this->source());
        GISMO_ASSERT( dynamic_cast<const gsMultiBasis<T>*>(&this->source()), "error");

        // Reconstruct solution coefficients on patch p
        const index_t sz  = mb[p].size();
        result.resize(sz, dim!=1 ? dim : solVector.cols()); // (!)

        for (index_t c = 0; c!=dim; c++) // for all components
        {
            // loop over all basis functions (even the eliminated ones)
            for (index_t i = 0; i < sz; ++i)
            {
                const int ii = m_sd->mapper.index(i, p, c);
                if ( m_sd->mapper.is_free_index(ii) ) // DoF value is in the solVector
                    // result.row(i) = solVector.row(ii);
                    result(i,c) = solVector.at(ii);
                else // eliminated DoF: fill with Dirichlet data
                {
                    result(i,c) =  m_sd->fixedDofs.at( m_sd->mapper.global_to_bindex(ii) );
                }
            }
        }
    }

    // space restrictTo(boundaries);
    // space restrictTo(bcRefList domain);

    void setupMapper(gsDofMapper dofsMapper) const
    {
        GISMO_ASSERT( dofsMapper.isFinalized(), "The provided dof-mapper is not finalized.");
        GISMO_ASSERT( dofsMapper.mapSize()==static_cast<size_t>(this->source().size()), "The dof-mapper is not consistent.");
        m_sd->mapper = give(dofsMapper);
    }

    void setup(const index_t _icont = -1) const
    {
        this->setInterfaceCont(_icont);
        m_sd->mapper = gsDofMapper();

        if (const gsMultiBasis<T> * mb =
            dynamic_cast<const gsMultiBasis<T>*>(&this->source()) )
        {
            m_sd->mapper = gsDofMapper(*mb, this->dim() );
            //m_mapper.init(*mb, this->dim()); //bug
            if ( 0==this->interfaceCont() ) // Conforming boundaries ?
            {
                for ( gsBoxTopology::const_iiterator it = mb->topology().iBegin();
                      it != mb->topology().iEnd(); ++it )
                {
                    mb->matchInterface(*it, m_sd->mapper);
                }
            }
        }

        if (const gsMappedBasis<2,T> * mb =
            dynamic_cast<const gsMappedBasis<2,T>*>(&this->source()) )
        {
            m_sd->mapper.setIdentity(mb->nPatches(), mb->size() , this->dim());
        }

        m_sd->mapper.finalize();
    }

    void setup(const gsBoundaryConditions<T> & bc, const index_t dir_values,
               const index_t _icont = -1) const
    {
        this->setInterfaceCont(_icont);
        m_sd->mapper = gsDofMapper();
        const gsMultiBasis<T> *mb = dynamic_cast<const gsMultiBasis<T> *>(&this->source());
        if (mb != nullptr)
        {
            m_sd->mapper = gsDofMapper(*mb, this->dim());
            //m_mapper.init(*mb, this->dim()); //bug
            if (0 == this->interfaceCont()) // Conforming boundaries ?
            {
                for (gsBoxTopology::const_iiterator it = mb->topology().iBegin();
                     it != mb->topology().iEnd(); ++it) {
                    mb->matchInterface(*it, m_sd->mapper);
                }
            }

            // Strong Dirichlet conditions
            gsMatrix<index_t> bnd;
            for (typename gsBoundaryConditions<T>::const_iterator
                     it = bc.begin("Dirichlet"); it != bc.end("Dirichlet"); ++it)
            {
                const index_t cc = it->unkComponent();
                GISMO_ASSERT(static_cast<size_t>(it->ps.patch) < this->mapper().numPatches(),
                             "Problem: a boundary condition is set on a patch id which does not exist.");

                bnd = mb->basis(it->ps.patch).boundary(it->ps.side());
                m_sd->mapper.markBoundary(it->ps.patch, bnd, cc);
            }
            // Clamped boundary condition (per DoF)
            gsMatrix<index_t> bnd1;
            for (typename gsBoundaryConditions<T>::const_iterator
                     it = bc.begin("Clamped"); it != bc.end("Clamped"); ++it)
            {
                const index_t cc = it->unkComponent();

                GISMO_ASSERT(static_cast<size_t>(it->ps.patch) < this->mapper().numPatches(),
                             "Problem: a boundary condition is set on a patch id which does not exist.");

                bnd = mb->basis(it->ps.patch).boundaryOffset(it->ps.side(), 0);
                bnd1 = mb->basis(it->ps.patch).boundaryOffset(it->ps.side(), 1);
                // Cast to tensor b-spline basis
                if (!it->ps.parameter())
                        bnd.swap(bnd1);
                for (index_t k = 0; k < bnd.size(); ++k)
                    m_sd->mapper.matchDof(it->ps.patch, (bnd)(k, 0),
                                          it->ps.patch, (bnd1)(k, 0), cc);
            }

            // Collapsed
            for (typename gsBoundaryConditions<T>::const_iterator
                     it = bc.begin("Collapsed"); it != bc.end("Collapsed"); ++it)
            {
                const index_t cc = it->unkComponent();

                GISMO_ASSERT(static_cast<size_t>(it->ps.patch) < this->mapper().numPatches(),
                             "Problem: a boundary condition is set on a patch id which does not exist.");

                bnd = mb->basis(it->ps.patch).boundary(it->ps.side());

                // match all DoFs to the first one of the side
                for (index_t k = 0; k < bnd.size() - 1; ++k)
                    m_sd->mapper.matchDof(it->ps.patch, (bnd)(0, 0),
                                          it->ps.patch, (bnd)(k + 1, 0), cc);
            }

            // corners
            for (typename gsBoundaryConditions<T>::const_citerator
                     it = bc.cornerBegin(); it != bc.cornerEnd(); ++it)
            {
                for (index_t r = 0; r!=this->dim(); ++r)
                {
                    if (it->component!=-1 && r!=it->component) continue;

                    //assumes (unk == -1 || it->unknown == unk)
                    GISMO_ASSERT(static_cast<size_t>(it->patch) < mb->nBases(),
                                 "Problem: a corner boundary condition is set on a patch id which does not exist.");
                    m_sd->mapper.eliminateDof(mb->basis(it->patch).functionAtCorner(it->corner),
                                              it->patch, it->component);
                }
            }

        } else if (const gsBasis<T> *b =
                   dynamic_cast<const gsBasis<T> *>(&this->source()))
        {
            m_sd->mapper = gsDofMapper(*b, this->dim() );
            gsMatrix<index_t> bnd;
            for (typename gsBoundaryConditions<T>::const_iterator
                     it = bc.begin("Dirichlet"); it != bc.end("Dirichlet"); ++it) {
                GISMO_ASSERT(it->ps.patch == 0,
                             "Problem: a boundary condition is set on a patch id which does not exist.");

                bnd = b->boundary(it->ps.side());
                m_sd->mapper.markBoundary(0, bnd, it->unkComponent());
            }

            for (typename gsBoundaryConditions<T>::const_iterator
                     it = bc.begin("Clamped"); it != bc.end("Clamped"); ++it) {
                GISMO_ASSERT(it->ps.patch == 0,
                             "Problem: a boundary condition is set on a patch id which does not exist.");

                bnd = b->boundary(it->ps.side());
                //const index_t cc = it->unkComponent();
                // m_sd->mapper.markBoundary(0, bnd, 0);
            }

            m_sd->mapper = gsDofMapper(*b);
            for (typename gsBoundaryConditions<T>::const_iterator
                     it = bc.begin("Collapsed"); it != bc.end("Collapsed"); ++it) {
                GISMO_ASSERT(it->ps.patch == 0,
                             "Problem: a boundary condition is set on a patch id which does not exist.");

                bnd = b->boundary(it->ps.side());
                //const index_t cc = it->unkComponent();
                // m_sd->mapper.markBoundary(0, bnd, 0);
            }
        } else if (const gsMappedBasis<2, T> *mapb =
                   dynamic_cast<const gsMappedBasis<2, T> *>(&this->source())) {
            m_sd->mapper.setIdentity(mapb->nPatches(), mapb->size(), this->dim());

            if (0 == this->interfaceCont()) // C^0 matching interface
            {
                gsMatrix<index_t> int1, int2;
                for (gsBoxTopology::const_iiterator it = mapb->getTopol().iBegin();
                     it != mapb->getTopol().iEnd(); ++it) {
                    int1 = mapb->basis(it->first().patch).boundaryOffset(it->first().side(), 0);
                    int2 = mapb->basis(it->second().patch).boundaryOffset(it->second().side(), 0);

                    m_sd->mapper.matchDofs(it->first().patch, int1, it->second().patch, int2);
                }
            }
            if (1 == this->interfaceCont()) // C^1 matching interface
            {
                GISMO_ERROR("Boundary offset function is not implemented for gsMappedBasis in general.");
            }

            gsMatrix<index_t> bnd;
            for (typename gsBoundaryConditions<T>::const_iterator
                     it = bc.begin("Dirichlet"); it != bc.end("Dirichlet"); ++it) {
                const index_t cc = it->unkComponent();
                GISMO_ASSERT(static_cast<size_t>(it->ps.patch) < this->mapper().numPatches(),
                             "Problem: a boundary condition is set on a patch id which does not exist.");

                bnd = mapb->basis(it->ps.patch).boundary(it->ps.side());
                m_sd->mapper.markBoundary(it->ps.patch, bnd, cc);
            }

            // Clamped boundary condition (per DoF)
            gsMatrix<index_t> bnd1;
            for (typename gsBoundaryConditions<T>::const_iterator
                     it = bc.begin("Clamped"); it != bc.end("Clamped"); ++it) {
                const index_t cc = it->unkComponent();

                GISMO_ASSERT(static_cast<size_t>(it->ps.patch) < this->mapper().numPatches(),
                             "Problem: a boundary condition is set on a patch id which does not exist.");

                bnd = mapb->basis(it->ps.patch).boundaryOffset(it->ps.side(), 0);
                bnd1 = mapb->basis(it->ps.patch).boundaryOffset(it->ps.side(), 1);

                // Cast to tensor b-spline basis
                if (mapb != NULL) // clamp adjacent dofs
                {
                    if (!it->ps.parameter())
                        bnd.swap(bnd1);
                    for (index_t k = 0; k < bnd.size(); ++k)
                        m_sd->mapper.matchDof(it->ps.patch, (bnd)(k, 0),
                                              it->ps.patch, (bnd1)(k, 0), cc);
                } else
                    gsWarn << "Unable to apply clamped condition.\n";
            }

            // COLLAPSED
            for (typename gsBoundaryConditions<T>::const_iterator
                     it = bc.begin("Collapsed"); it != bc.end("Collapsed"); ++it) {
                const index_t cc = it->unkComponent();

                GISMO_ASSERT(static_cast<size_t>(it->ps.patch) < this->mapper().numPatches(),
                             "Problem: a boundary condition is set on a patch id which does not exist.");

                bnd = mapb->basis(it->ps.patch).boundary(it->ps.side());

                // Cast to tensor b-spline basis
                if (mapb != NULL) // clamp adjacent dofs
                {
                    // match all DoFs to the first one of the side
                    for (index_t k = 0; k < bnd.size() - 1; ++k)
                        m_sd->mapper.matchDof(it->ps.patch, (bnd)(0, 0),
                                              it->ps.patch, (bnd)(k + 1, 0), cc);
                }
            }

            // corners
            for (typename gsBoundaryConditions<T>::const_citerator
                     it = bc.cornerBegin(); it != bc.cornerEnd(); ++it) {
                //assumes (unk == -1 || it->unknown == unk)
                GISMO_ASSERT(static_cast<size_t>(it->patch) < mb->nBases(),
                             "Problem: a corner boundary condition is set on a patch id which does not exist.");
                m_sd->mapper.eliminateDof(mapb->basis(it->patch).functionAtCorner(it->corner), it->patch, it->component);
            }
        } else
        {
            GISMO_ASSERT(0 == bc.size(), "Problem: BCs are ignored.");
            m_sd->mapper.setIdentity(this->source().nPieces(), this->source().size());
        }

        m_sd->mapper.finalize();

        // Compute Dirichlet node values
        gsDirichletValues(bc, dir_values, *this);
    }

    void print(std::ostream &os) const { os << "u"; }

protected:
    friend class gismo::gsExprHelper<Scalar>;
    friend class symbol_expr<gsFeSpace>;
    explicit gsFeSpace(index_t _d = 1) : Base(_d), m_sd(nullptr) { }
};

template<class T> inline bool
operator== (const gsFeSpace<T> & a, const gsFeSpace<T> & b)
{ return a.id()== b.id() && a.isAcross()==b.isAcross(); }

/*
  Expression representing a function given by a vector of
  coefficients in a gsFeSpace.

  Typically it used for accessing the solution of a boundary-value
  problem.
*/
template<class T>
class gsFeSolution : public _expr<gsFeSolution<T> >
{
protected:
    const gsFeSpace<T> _u;
    gsMatrix<T> * _Sv; ///< Pointer to a coefficient vector
    bool m_isAcross; ///< true when this expression is evaluated across an interface

public:
    typedef T Scalar;
    enum {Space = 0, ScalarValued= 0, ColBlocks= 0};

    bool isAcross() const { return m_isAcross; }

    gsFeSolution right() const
    {
        gsFeSolution ac(*this);
        ac.m_isAcross = true;
        return ac;
    }

    gsFeSolution left() const { return gsFeSolution(*this); }

    explicit gsFeSolution(const gsFeSpace<T> & u) : _u(u), _Sv(NULL) { }

    gsFeSolution(const gsFeSpace<T> & u, gsMatrix<T> & Sv) : _u(u), _Sv(&Sv) { }

    const gsFeSpace<T> & space() const {return _u;};

    mutable gsMatrix<T> res;
    const gsMatrix<T> & eval(index_t k) const
    {
        GISMO_ASSERT(1==_u.data().actives.cols(), "Single actives expected");

        res.setZero(_u.dim(), 1);
        const gsDofMapper & map = _u.mapper();
        GISMO_ASSERT(_Sv->size()==map.freeSize(), "The solution vector has wrong dimensions: "<<_Sv->size()<<" != "<<map.freeSize());

        for (index_t c = 0; c!=_u.dim(); c++) // for all components
        {
            for (index_t i = 0; i!=_u.data().actives.size(); ++i)
            {
                const index_t ii = map.index(_u.data().actives.at(i), _u.data().patchId, c);
                if ( map.is_free_index(ii) ) // DoF value is in the solVector
                    res.at(c) += _Sv->at(ii) * _u.data().values[0](i,k);
                else
                    res.at(c) += _u.data().values[0](i,k) *
                        _u.fixedPart().at( map.global_to_bindex(ii) );
            }
        }
        return res;
    }

    //template<class U>
    //linearComb(U & ie){ sum up ie[_u] times the _Sv  }
    // ie.eval(k), _u.data().actives(), fixedPart() - see lapl_expr

    const gsFeSpace<Scalar> & rowVar() const {return gsNullExpr<Scalar>::get();}
    const gsFeSpace<Scalar> & colVar() const {return gsNullExpr<Scalar>::get();}

    index_t rows() const {return _u.dim(); }

    static index_t cols() {return 1; }

    void parse(gsExprHelper<Scalar> & evList) const
    {
        evList.add(_u);
        _u.data().flags |= NEED_VALUE | NEED_ACTIVE;
    }

    void print(std::ostream &os) const { os << "s"; }

public:
    index_t dim() const { return _u.dim();}

    index_t parDim() const
    { return _u.source().domainDim(); }

    gsDofMapper & mapper() {return _u.mapper();}
    const gsDofMapper & mapper() const {return _u.mapper();}

    inline const gsMatrix<T> & fixedPart() const {return _u.fixedPart();}
    gsMatrix<T> & fixedPart() {return _u.fixedPart();}

    gsFuncData<T> & data() {return *_u.data();}
    const gsFuncData<T> & data() const {return _u.data();}

    void setSolutionVector(gsMatrix<T>& solVector)
    { _Sv = & solVector; }

    const gsMatrix<T> & coefs() const { return *_Sv; }
    //gsMatrix<T> & coefs() { return *_Sv; } // wd4702 ?

    /// val: perturbation value, j: local bf index, p: patch
    void perturbLocal(T val, index_t j, index_t p = 0)
    {
        GISMO_ASSERT(1==_u.data().actives.cols(), "Single actives expected");

        auto qr = std::div(j, _u.data().actives.size() );
        _u.mapper().print();
        gsDebugVar(_u.mapper().asVector());
        const index_t ii = _u.mapper().index(qr.rem, p, qr.quot);
        gsDebugVar(ii);
        gsDebugVar(j);
        gsDebugVar(qr.rem );
        gsDebugVar(qr.quot);
        gsDebugVar( _Sv);
        if (_u.mapper().is_free_index(ii) )
        {
            GISMO_ASSERT(ii<_Sv->size(), "Solution vector is not initialized/allocated, sz="<<_Sv->size() );
            gsDebugVar( _Sv);
            _Sv->at(ii) += val;
        }
        //else
        //    _u.fixedPart().at( _u.mapper().global_to_bindex(ii) ) += val;
    }

    /// Extract the coefficients of piece \a p
    void extract(gsMatrix<T> & result, const index_t p = 0) const
    { _u.getCoeffs(*_Sv, result, p); }

    /// Extracts ALL the coefficients in a solution vector; including
    /// coupled and boundary DoFs
    void extractFull(gsMatrix<T> & result) const
    {
        index_t offset;
        const index_t dim = _u.dim();
        const size_t totalSz = _u.mapper().mapSize();
        result.resize(totalSz, 1);
        for (size_t p=0; p!=_u.mapper().numPatches(); ++p)
        {
            offset = _u.mapper().offset(p);
            // Reconstruct solution coefficients on patch p

            for (index_t c = 0; c!=dim; c++) // for all components
            {
                const index_t sz  = _u.mapper().patchSize(p,c);

                // loop over all basis functions (even the eliminated ones)
                for (index_t i = 0; i < sz; ++i)
                {
                    //gsDebugVar(i);
                    const int ii = _u.mapper().index(i, p, c);
                    //gsDebugVar(ii);
                    if ( _u.mapper().is_free_index(ii) ) // DoF value is in the solVector
                    {
                        result(i+offset,0) = _Sv->at(ii);
                    }
                    else // eliminated DoF: fill with Dirichlet data
                        result(i+offset,0) =  _u.fixedPart().at( _u.mapper().global_to_bindex(ii) );
                }
                offset += sz;
            }
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
                  fixedPart.row(m_sd->mapper.global_to_bindex(ii)) = cf.row(i);
                  }
                */
            }
        }
    }
};

/*
  Expression for the transpose of an expression
*/
template<class E, bool cw>
class tr_expr : public _expr<tr_expr<E,cw> >
{
    typename E::Nested_t _u;

public:

    typedef typename E::Scalar Scalar;

    tr_expr(_expr<E> const& u)
    : _u(u) { }

public:
    enum {ColBlocks = E::ColBlocks, ScalarValued=E::ScalarValued};
    enum {Space = cw?E::Space:(E::Space==1?2:(E::Space==2?1:E::Space))};

    mutable Temporary_t res;
    const Temporary_t & eval(const index_t k) const
    {
        if (E::ColBlocks)
            res = _u.eval(k).blockTranspose( _u.cardinality() );
        else
            res = _u.eval(k).transpose();
        return res;
    }

    index_t rows() const { return _u.cols(); }

    index_t cols() const { return _u.rows(); }

    void parse(gsExprHelper<Scalar> & evList) const
    { _u.parse(evList); }

    const gsFeSpace<Scalar> & rowVar() const { return cw?_u.rowVar():_u.colVar(); }
    const gsFeSpace<Scalar> & colVar() const { return cw?_u.colVar():_u.rowVar(); }

    index_t cardinality_impl() const { return _u.cardinality_impl(); }

    void print(std::ostream &os) const { os<<"("; _u.print(os); os <<")\u1D40"; }
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
  Expression to make an expression colblocks
*/
template<class E>
class cb_expr : public _expr<cb_expr<E> >
{
    typename E::Nested_t _u;

public:

    typedef typename E::Scalar Scalar;

    cb_expr(_expr<E> const& u)
    : _u(u) { }

public:
    enum {ColBlocks = 1, ScalarValued=E::ScalarValued};
    enum {Space = E::Space};

    mutable gsMatrix<Scalar> ev, res;

    const gsMatrix<Scalar> & eval(const index_t k) const
    {
        //return _u.eval(k).transpose();
        // /*
        ev = _u.eval(k);
        if (E::ColBlocks)
        {
            return ev;
        }
        else
        {
            res = _u.eval(k);

            GISMO_ASSERT(res.rows() % _u.rows() == 0 && res.cols() % _u.cols() == 0,"Result is not a multiple of the space dimensions?");

            index_t cardinality;
            if ( (cardinality = res.rows() / _u.rows()) >= 1 && res.cols() / _u.cols() == 1 ) // stored in rows
            {
                res.resize(_u.rows(), cardinality * _u.cols());
                for (index_t r = 0; r!=cardinality; r++)
                    res.block(0 , r * _u.cols(), _u.rows(), _u.cols()) = ev.block( r * _u.rows(), 0, _u.rows(), _u.cols() );
            }
            else if ( (cardinality = res.rows() / _u.rows()) == 1 && res.cols() / _u.cols() >= 1 ) // stored in cols ----->>>> This is already colBlocks???
            {
                res.resize(_u.rows(), cardinality * _u.cols());
                for (index_t r = 0; r!=cardinality; r++)
                    res.block(0 , r * _u.cols(), _u.rows(), _u.cols()) = ev.block( 0, r * _u.cols(), _u.rows(), _u.cols() );
            }
        }
        return res;
    }

    index_t rows() const { return _u.rows(); }

    index_t cols() const { return _u.cols(); }

    void parse(gsExprHelper<Scalar> & evList) const
    { _u.parse(evList); }

    const gsFeSpace<Scalar> & rowVar() const { return _u.rowVar(); }
    const gsFeSpace<Scalar> & colVar() const { return _u.colVar(); }

    index_t cardinality_impl() const
    {
        res = _u.eval(0);
        index_t cardinality;
        if ( res.rows() / _u.rows() >= 1 && res.cols() / _u.cols() == 1 ) // stored in rows
            cardinality = res.rows() / _u.rows();
        else if ( res.rows() / _u.rows() == 1 && res.cols() / _u.cols() >= 1 )
            cardinality = res.cols() / _u.cols();
        else
            GISMO_ERROR("Cardinality for cb_expr cannot be determined.");

        return cardinality;
    }

    void print(std::ostream &os) const { os<<"{"; _u.print(os); os <<"}"; }
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
    enum {Space = E::Space, ScalarValued = E::ScalarValued,
        ColBlocks = E::ColBlocks};

    // template<bool S  = ColBlocks>
    // typename util::enable_if<S,MatExprType>::type
    const gsMatrix<Scalar> & eval(const index_t k) const
    {
        tmp = _u.eval(k);
        return tmp;
    }

    index_t rows() const { return _u.rows(); }
    index_t cols() const { return _u.cols(); }
    void parse(gsExprHelper<Scalar> & evList) const { _u.parse(evList); }
    const gsFeSpace<Scalar> & rowVar() const { return _u.rowVar(); }
    const gsFeSpace<Scalar> & colVar() const { return _u.colVar(); }

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
    enum {ScalarValued = 0, Space = E::Space, ColBlocks= 0};

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
        auto tmp = _u.eval(k);
        const index_t cb = _u.rows();
        const index_t r  = _u.cardinality();
        if (Space==1)
            res.resize(r, 1);
        else
            res.resize(1, r);

        for (index_t i = 0; i!=r; ++i)
            res.at(i) = tmp.middleCols(i*cb,cb).trace();
        return res;
    }

    // choose if !ColBlocks
    //todo: Scalar eval(const index_t k) const

    index_t rows() const { return _u.cols() / _u.rows(); } //_u.cardinality()?
    index_t cols() const { return 1; }

    index_t cardinality_impl() const { return _u.cardinality(); }

    void parse(gsExprHelper<Scalar> & evList) const
    { _u.parse(evList); }

    const gsFeSpace<Scalar> & rowVar() const { return _u.rowVar(); }
    const gsFeSpace<Scalar> & colVar() const { return _u.colVar(); }

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
        auto tmp = _u.eval(k);
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

    void parse(gsExprHelper<Scalar> & evList) const
    { _u.parse(evList); }

    const gsFeSpace<Scalar> & rowVar() const { return _u.rowVar(); }
    const gsFeSpace<Scalar> & colVar() const { return _u.colVar(); }

    void print(std::ostream &os) const { os << "adj("; _u.print(os); os<<")"; }
};

template<class E>
class reshape_expr  : public _expr<reshape_expr<E> >
{
public:
    typedef typename E::Scalar Scalar;
    enum {ScalarValued = 0, ColBlocks = E::ColBlocks};
    enum {Space = E::Space};
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
        // Note: this assertion would fail in the constructor!
        GISMO_ASSERT( _u.rows()*_u.cols() == _n*_m, "Wrong dimension "
                      << _u.rows() << " x "<<_u.cols() << "!=" << _n << " * "<< _m );
        tmp = _u.eval(k);
        return gsAsConstMatrix<Scalar>(tmp.data(),_n,_m);
    }

    index_t rows() const { return _n; }
    index_t cols() const { return _m; }

    void parse(gsExprHelper<Scalar> & evList) const
    { _u.parse(evList); }

    const gsFeSpace<Scalar> & rowVar() const { return _u.rowVar(); }
    const gsFeSpace<Scalar> & colVar() const { return _u.colVar(); }

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
    enum {ScalarValued = 0, Space = E::Space, ColBlocks= E::ColBlocks};
private:
    typename E::Nested_t _u;
    index_t _n, _m;
    mutable gsMatrix<Scalar> tmp;

public:

    //the replicate is done nxm times
    replicate_expr(_expr<E> const& u, index_t n, index_t m) : _u(u), _n(n), _m(m)
    {
    }

    auto eval(const index_t k) const -> decltype(tmp.replicate(_n,_m))
    {
        tmp = _u.eval(k);
        return tmp.replicate(_n,_m);
    }

    index_t rows() const { return _n*_u.rows(); }
    index_t cols() const { return _m*_u.cols(); }

    void parse(gsExprHelper<Scalar> & evList) const
    { _u.parse(evList); }

    const gsFeSpace<Scalar> & rowVar() const { return _u.rowVar(); }
    const gsFeSpace<Scalar> & colVar() const { return _u.colVar(); }
    index_t cardinality_impl() const { return _u.cardinality_impl(); }

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
    enum {ScalarValued = 0, Space = E::Space, ColBlocks= 0}; // to do: ColBlocks
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
        tmp = _u.eval(k);
        const index_t numActives = _u.cardinality();

        for (index_t i = 0; i<numActives; ++i)
        {
            tmp(0,2*i+1) += tmp(1,2*i);
            std::swap(tmp(1,2*i), tmp(1,2*i+1));
        }

        tmp.resize(4,numActives);
        tmp.conservativeResize(3,numActives);

        if ( 1==Space )
            tmp.transposeInPlace();
        else if (2!=Space) // if not colSpan and not rowSpan
            tmp.transposeInPlace();

        return tmp;
    }

    index_t rows() const { return 1; }
    index_t cols() const { return 3; }

    void parse(gsExprHelper<Scalar> & evList) const
    { _u.parse(evList); }

    const gsFeSpace<Scalar> & rowVar() const { return _u.rowVar(); }
    const gsFeSpace<Scalar> & colVar() const { return _u.colVar(); }
    index_t cardinality_impl() const { return _u.cardinality_impl(); }

    void print(std::ostream &os) const { os << "flat("; _u.print(os); os<<")"; }
};

/// Make a matrix 2x2 expression "flat"
template <typename E> EIGEN_STRONG_INLINE
flat_expr<E> const flat(E const & u)
{ return flat_expr<E>(u); }


/*
  Expression for the diagonal(s) of a (matrix) expression
*/
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
        { 
            GISMO_ASSERT(0== _u.cols()%_u.rows(), "Expecting square-block expression, got "
            << _u.rows() <<" x "<< _u.cols() ); 
        }

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

        void parse(gsExprHelper<Scalar> & evList) const
        { _u.parse(evList); }

        const gsFeSpace<Scalar> & rowVar() const { return _u.rowVar(); }
        const gsFeSpace<Scalar> & colVar() const { return _u.colVar(); }


  void print(std::ostream &os) const { os << "diag("; _u.print(os); os<<")"; }
  };

/// Get diagonal elements of matrix as a vector
template <typename E> EIGEN_STRONG_INLINE
diag_expr<E> const diagonal(E const & u)
{ return diag_expr<E>(u); }

#define GISMO_EXPR_VECTOR_EXPRESSION(name, mname, isSv)                 \
    template<class E> class name##_##expr  : public _expr<name##_##expr<E> > { \
        typename E::Nested_t _u;                                        \
    public:                                                             \
    typedef typename E::Scalar Scalar;                                  \
    enum {Space= E::Space, ScalarValued= isSv, ColBlocks= E::ColBlocks}; \
    name##_##expr(_expr<E> const& u) : _u(u) { }                        \
    mutable Temporary_t tmp;                                            \
    const Temporary_t & eval(const index_t k) const {                   \
        tmp = _u.eval(k).mname(); return tmp; }                         \
    index_t rows() const { return isSv ? 0 : _u.rows(); }               \
    index_t cols() const { return isSv ? 0 : _u.cols(); }               \
    void parse(gsExprHelper<Scalar> & evList) const { _u.parse(evList); } \
    const gsFeSpace<Scalar> & rowVar() const {return gsNullExpr<Scalar>::get();} \
    const gsFeSpace<Scalar> & colVar() const {return gsNullExpr<Scalar>::get();} \
    void print(std::ostream &os) const                                  \
        { os << #name <<"("; _u.print(os); os <<")"; }                  \
    };

/// Eucledian Norm
GISMO_EXPR_VECTOR_EXPRESSION(norm,norm,1);
/// Squared Eucledian Norm
GISMO_EXPR_VECTOR_EXPRESSION(sqNorm,squaredNorm,1);
/// Normalization of a vector to unit measure
GISMO_EXPR_VECTOR_EXPRESSION(normalized,normalized,0);
/// Inverse of a matrix expression
GISMO_EXPR_VECTOR_EXPRESSION(inv,cramerInverse,0);
// GISMO_EXPR_VECTOR_EXPRESSION(cwSqr,array().square,0)
// GISMO_EXPR_VECTOR_EXPRESSION(sum,array().sum,1)
// GISMO_EXPR_VECTOR_EXPRESSION(sqrt,array().sqrt,0)
//GISMO_EXPR_VECTOR_EXPRESSION(abs,array().abs,0)

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
    enum{Space = E::Space, ScalarValued= 0, ColBlocks= E::ColBlocks};

    asdiag_expr(_expr<E> const& u) : _u(u) { }

public:

    const gsMatrix<Scalar> & eval(const index_t k) const
    {
        auto m = _u.eval(k);
        const index_t r = m.rows();
        const index_t c = m.cols();
        res.resize(r,r*c);
        for (index_t i = 0; i!=c; ++i)
            res.middleCols(i*r,r) = m.col(i).asDiagonal();
        return res;
    }

    const gsFeSpace<Scalar> & rowVar() const { return _u.rowVar(); }
    const gsFeSpace<Scalar> & colVar() const { return _u.colVar(); }

    index_t cardinality_impl() const { return _u.cardinality_impl(); }

    index_t rows() const { return _u.rows(); }
    index_t cols() const { return _u.rows() * _u.cols(); }
    void parse(gsExprHelper<Scalar> & evList) const
    { _u.parse(evList); }

    void print(std::ostream &os) const { os << "diag("; _u.print(os); os <<")";}
};

// Takes the max of a vector
template<class E>
class max_expr  : public _expr<max_expr<E> >
{
public:
    typedef typename E::Scalar Scalar;
    enum {ScalarValued = 0, Space = E::Space, ColBlocks = 1};
private:
    typename E::Nested_t _u;
    mutable gsMatrix<Scalar> tmp;
    mutable gsMatrix<Scalar> res;

public:

    max_expr(_expr<E> const& u) : _u(u)
    {
        //GISMO_ASSERT( _u.rows()*_u.cols() == _n*_m, "Wrong dimension"); //
    }

    const gsMatrix<Scalar> & eval(const index_t k) const {return eval_impl(_u,k); }

    index_t rows() const { return 1; }
    index_t cols() const { return 1; }
    void setFlag() const { _u.setFlag(); }

    void parse(gsExprHelper<Scalar> & evList) const
    { _u.parse(evList); }

    const gsFeSpace<Scalar> & rowVar() const { return _u.rowVar(); }
    const gsFeSpace<Scalar> & colVar() const { return _u.colVar(); }
    index_t cardinality_impl() const { return _u.cardinality_impl(); }

    void print(std::ostream &os) const { os << "max("; _u.print(os); os<<")"; }
private:
    template<class U> inline
    typename util::enable_if< util::is_same<U,gsFeSpace<Scalar> >::value, const gsMatrix<Scalar> & >::type
    eval_impl(const U & u, const index_t k)  const
    {
        tmp = u.eval(k);

        res.resize(1,u.cardinality());
        if (E::ColBlocks)
            for (index_t c=0; c!=_u.cardinality(); c++)
                res(0,c) = tmp.block(0,c*u.cols(),u.rows(),u.cols()).maxCoeff();
        else
            for (index_t c=0; c!=_u.rows(); c++)
                res(0,c) = tmp.block(c*u.rows(),0,u.rows(),u.cols()).maxCoeff();
        return res;
    }


    template<class U> inline
    typename util::enable_if< !util::is_same<U,gsFeSpace<Scalar> >::value, const gsMatrix<Scalar> & >::type
    eval_impl(const U & u, const index_t k)  const
    {
        res = u.eval(k).colwise().maxCoeff();
        return res;
    }
};

template<class E>
class rowsum_expr  : public _expr<rowsum_expr<E> >
{
public:
    typedef typename E::Scalar Scalar;
    enum {ScalarValued = 0, Space = E::Space, ColBlocks = E::ColBlocks};
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

    void parse(gsExprHelper<Scalar> & evList) const
    { _u.parse(evList); }

    const gsFeSpace<Scalar> & rowVar() const { return _u.rowVar(); }
    const gsFeSpace<Scalar> & colVar() const { return _u.colVar(); }
    index_t cardinality_impl() const { return _u.cardinality_impl(); }

    void print(std::ostream &os) const { os << "rowsum("; _u.print(os); os<<")"; }
};

template<class E>
class colsum_expr  : public _expr<colsum_expr<E> >
{
public:
    typedef typename E::Scalar Scalar;
    enum {ScalarValued = 0, Space = E::Space, ColBlocks = E::ColBlocks};
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

    void parse(gsExprHelper<Scalar> & evList) const
    { _u.parse(evList); }

    const gsFeSpace<Scalar> & rowVar() const { return _u.rowVar(); }
    const gsFeSpace<Scalar> & colVar() const { return _u.colVar(); }
    index_t cardinality_impl() const { return _u.cardinality_impl(); }

    void print(std::ostream &os) const { os << "colsum("; _u.print(os); os<<")"; }
};

/**
   Expression for the identity matrix
*/
class idMat_expr : public _expr<idMat_expr >
{
public:
    typedef real_t Scalar;
    enum {Space = 0, ScalarValued = 0, ColBlocks = 0};
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
    void parse(gsExprHelper<Scalar> & ) const {  }

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
    typename E::Scalar _tol;
public:
    typedef typename E::Scalar Scalar;
    enum {ScalarValued = 1, Space = E::Space, ColBlocks= 0};

    sign_expr(_expr<E> const& u, Scalar tolerance = 0.0) : _u(u),_tol(tolerance){ 
        GISMO_ASSERT( _tol >= 0, "Tolerance for sign_expr should be a positive number.");
    }

    Scalar eval(const index_t k) const
    {
        const Scalar v = _u.val().eval(k);
        return ( v>_tol ? 1 : ( v<-_tol ? -1 : 0 ) );
    }

    static index_t rows() { return 0; }
    static index_t cols() { return 0; }

    void parse(gsExprHelper<Scalar> & el) const
    { _u.parse(el); }

    static bool isScalar() { return true; }

    const gsFeSpace<Scalar> & rowVar() const {return gsNullExpr<Scalar>::get();}
    const gsFeSpace<Scalar> & colVar() const {return gsNullExpr<Scalar>::get();}

    void print(std::ostream &os) const { os<<"sgn("; _u.print(os); os <<")"; }
};

/**
   Expression for the component-wise positive part
*/
template<class E>
class ppart_expr : public _expr<ppart_expr<E> >
{
public:
    typedef typename E::Scalar Scalar;
    enum {ScalarValued = E::ScalarValued, Space = E::Space, ColBlocks= E::ColBlocks};
private:
    typename E::Nested_t _u;
    mutable gsMatrix<Scalar> res;
public:

    ppart_expr(_expr<E> const& u) : _u(u) { }

    const gsMatrix<Scalar> & eval(index_t k) const
    {
        res = _u.eval(k).cwiseMax(0.0); // component-wise maximum with zero
        return res;
    }


    const index_t rows() const { return _u.rows(); }
    const index_t cols() const { return _u.cols(); }

    void parse(gsExprHelper<Scalar> & el) const
    { _u.parse(el); }

    const gsFeSpace<Scalar> & rowVar() const {return _u.rowVar();}
    const gsFeSpace<Scalar> & colVar() const {return _u.colVar();}

    void print(std::ostream &os) const { os<<"posPart("; _u.print(os); os <<")"; }
};


template<class E>
class pow_expr : public _expr<pow_expr<E> >
{
    typename E::Nested_t _u;

public:
    typedef typename E::Scalar Scalar;
    enum {ScalarValued = 1, Space = E::Space, ColBlocks= E::ColBlocks};

    Scalar _q;// power

    pow_expr(_expr<E> const& u, Scalar q) : _u(u), _q(q) { }

    Scalar eval(const index_t k) const
    {
        const Scalar v = _u.val().eval(k);
        return math::pow(v,_q);
    }

    static index_t rows() { return 0; }
    static index_t cols() { return 0; }

    void parse(gsExprHelper<Scalar> & el) const
    { _u.parse(el); }

    static bool isScalar() { return true; }

    const gsFeSpace<Scalar> & rowVar() const {return gsNullExpr<Scalar>::get();}
    const gsFeSpace<Scalar> & colVar() const {return gsNullExpr<Scalar>::get();}

    void print(std::ostream &os) const { os<<"pow("; _u.print(os); os <<")"; }
};

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

        const auto uEv  = _u.eval(k);
        const auto vEv  = _v.eval(k);

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

    void parse(gsExprHelper<Scalar> & evList) const
    { _u.parse(evList); _v.parse(evList); }

    const gsFeSpace<Scalar> & rowVar() const { return _v.rowVar(); }
    const gsFeSpace<Scalar> & colVar() const { return _v.colVar(); }

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
        const index_t r  = _u.rows();
        const index_t N  = _v.cols() / (r*r);

        const auto uEv  = _u.eval(k);
        const auto vEv  = _v.eval(k);

        res.resize(r, N*r*r);
        for (index_t s = 0; s!=r; ++s)
            for (index_t i = 0; i!=N; ++i)
            {
                res.middleCols((s*N + i)*r,r).noalias() =
                    uEv.transpose()*vEv.middleCols((s*N + i)*r,r).transpose();
            }
        //meaning: [Jg Jg Jg] * Jb ..
        return res;
    }

    index_t rows() const { return _u.cols(); }
    index_t cols() const { return _v.cols(); }

    void parse(gsExprHelper<Scalar> & evList) const
    { _u.parse(evList); _v.parse(evList); }

    const gsFeSpace<Scalar> & rowVar() const { return _v.rowVar(); }
    const gsFeSpace<Scalar> & colVar() const { return _v.colVar(); }

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
    enum {Space= 0, ScalarValued= 1, ColBlocks= 0};

    Scalar eval(const index_t k) const { return eval_impl(_u,k); }

    // enables expr.val().val()
    inline value_expr<E> val() const { return *this; }
    index_t rows() const { return 0; }
    index_t cols() const { return 0; }
    void parse(gsExprHelper<Scalar> & evList) const
    { _u.parse(evList); }

    static bool isScalar() { return true; }

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

template<class E>
class abs_expr  : public _expr<abs_expr<E> >
{
    typename E::Nested_t _u;

public:
    typedef typename E::Scalar Scalar;
    explicit abs_expr(_expr<E> const& u) : _u(u) { }

public:
    enum {Space= 0, ScalarValued= 1, ColBlocks= 0};

    Scalar eval(const index_t k) const { return abs_expr::eval_impl(_u,k); }

    index_t rows() const { return _u.rows(); }
    index_t cols() const { return _u.cols(); }
    void parse(gsExprHelper<Scalar> & evList) const
    { _u.parse(evList); }

    static bool isScalar() { return true; }

    const gsFeSpace<Scalar> & rowVar() const {return gsNullExpr<Scalar>::get();}
    const gsFeSpace<Scalar> & colVar() const {return gsNullExpr<Scalar>::get();}

    void print(std::ostream &os) const { _u.print(os); }

    // Math functions. eg
    // sqrt_mexpr<T> sqrt() { return sqrt_mexpr<T>(*this); }
private:
    template<class U> static inline
    typename util::enable_if<U::ScalarValued,Scalar>::type
    eval_impl(const U & u, const index_t k) {return math::abs(u.eval(k)); }
    template<class U> static inline
    typename util::enable_if<!U::ScalarValued,gsMatrix<Scalar> >::type
    eval_impl(const U & u, const index_t k) { return u.eval(k).cwiseAbs(); }
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
    enum {Space = E::Space, ScalarValued= 0, ColBlocks= 0};

    typedef typename E::Scalar Scalar;
    mutable gsMatrix<Scalar> tmp;

    grad_expr(const E & u) : _u(u)
    { GISMO_ASSERT(1==u.dim(),"grad(.) requires 1D variable, use jac(.) instead.");}

    const gsMatrix<Scalar> & eval(const index_t k) const
    {
        // Assumes: derivatives are in _u.data().values[1]
        // gsExprHelper acounts for compositions/physical expressions
        // so that derivs are directly computed
        tmp = _u.data().values[1].reshapeCol(k, cols(), cardinality_impl()).transpose();
        return tmp;
    }

    index_t rows() const { return 1 /*==u.dim()*/; }

    index_t cols() const { return _u.source().domainDim(); }

    index_t cardinality_impl() const
    { return _u.data().values[1].rows() / cols(); }

    void parse(gsExprHelper<Scalar> & evList) const
    {
        evList.add(_u);
        _u.data().flags |= NEED_GRAD;
    }

    const gsFeSpace<Scalar> & rowVar() const { return _u.rowVar(); }
    const gsFeSpace<Scalar> & colVar() const
    {return gsNullExpr<Scalar>::get();}

    void print(std::ostream &os) const { os << "\u2207("; _u.print(os); os <<")"; }
private:

    template<class U> static inline
    typename util::enable_if<util::is_same<U,gsComposition<Scalar> >::value,
                             const gsMatrix<Scalar> &>::type
    eval_impl(const U & u, const index_t k)
    {
        return u.eval(k);
    }
};

/*
  \brief Expression for the gradient of a finite element variable

  Transposed gradient vectors are returned as a matrix.
  This specialization is for a gsFeSolution object
*/

template<class T>
class grad_expr<gsFeSolution<T> > : public _expr<grad_expr<gsFeSolution<T> > >
{
protected:
    const gsFeSolution<T> _u;

public:
    typedef T Scalar;
    enum {Space = 0, ScalarValued= 0, ColBlocks= 0};

    explicit grad_expr(const gsFeSolution<T> & u) : _u(u) { }

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
                        _u.data().values[1].col(k).segment(i*_u.parDim(), _u.parDim()).transpose();
                }
                else
                {
                    res.row(c) +=
                        _u.fixedPart().at( map.global_to_bindex(ii) ) *
                        _u.data().values[1].col(k).segment(i*_u.parDim(), _u.parDim()).transpose();
                }
            }
        }
        return res;
    }

    index_t rows() const {return _u.dim();}
    index_t cols() const {return _u.parDim(); }

    const gsFeSpace<Scalar> & rowVar() const
    {return gsNullExpr<Scalar>::get();}
    const gsFeSpace<Scalar> & colVar() const
    {return gsNullExpr<Scalar>::get();}

    void parse(gsExprHelper<Scalar> & evList) const
    {
        _u.parse(evList);                         // add symbol
        evList.add(_u.space());
        _u.data().flags |= NEED_GRAD|NEED_ACTIVE; // define flags
    }

    void print(std::ostream &os) const { os << "\u2207(s)"; }
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
    enum{ Space = E::Space, ScalarValued = 0, ColBlocks = (1==E::Space?1:0)};

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

    index_t rows() const { return _u.source().domainDim(); }

    index_t cols() const { return _u.source().domainDim()*_u.rows(); }

    void parse(gsExprHelper<Scalar> & evList) const
    {
        evList.add(_u);
        _u.data().flags |= NEED_GRAD;
    }

    const gsFeSpace<Scalar> & rowVar() const { return _u.rowVar(); }
    const gsFeSpace<Scalar> & colVar() const
    {return gsNullExpr<Scalar>::get();}

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
    enum{Space = 1};

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

    void parse(gsExprHelper<Scalar> & evList) const
    {
        evList.add(u);
        u.data().flags |= NEED_GRAD;
    }

    const gsFeSpace<T> & rowVar() const { return u.rowVar(); }
    const gsFeSpace<T> & colVar() const
    {return gsNullExpr<Scalar>::get();}

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
    enum{Space = 1};

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

    void parse(gsExprHelper<Scalar> & evList) const
    {
        evList.add(u);
        u.data().flags |= NEED_DERIV2;
    }

    const gsFeSpace<T> & rowVar() const { return u.rowVar(); }
    const gsFeSpace<T> & colVar() const
    {return gsNullExpr<T>::get();}
};

/// The nabla2 (\f$\nabla^2\f$) of a finite element variable
template<class T>
nabla2_expr<T> nabla2(const gsFeVariable<T> & u) { return nabla2_expr<T>(u); }
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
    enum {Space = 0, ScalarValued= 0, ColBlocks= 0};

    onormal_expr(const gsGeometryMap<T> & G) : _G(G) { }

    auto eval(const index_t k) const -> decltype(_G.data().outNormals.col(k))
    { return _G.data().outNormals.col(k); }

    index_t rows() const { return _G.data().dim.second; }
    index_t cols() const { return 1; }

    const gsFeSpace<T> & rowVar() const {return gsNullExpr<T>::get();}
    const gsFeSpace<T> & colVar() const {return gsNullExpr<T>::get();}

    void parse(gsExprHelper<Scalar> & evList) const
    {
        evList.add(_G);
        _G.data().flags |= NEED_OUTER_NORMAL;
    }

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
    enum {Space = 0, ScalarValued= 0, ColBlocks= 0};

    normal_expr(const gsGeometryMap<T> & G) : _G(G) { }

    auto eval(const index_t k) const -> decltype(_G.data().normals.col(k))
    { return _G.data().normals.col(k); }

    index_t rows() const { return _G.data().dim.second; }
    index_t cols() const { return 1; }

    const gsFeSpace<T> & rowVar() const {return gsNullExpr<T>::get();}
    const gsFeSpace<T> & colVar() const {return gsNullExpr<T>::get();}

    void parse(gsExprHelper<Scalar> & evList) const
    {
        evList.add(_G);
        _G.data().flags |= NEED_NORMAL;
    }

    void print(std::ostream &os) const { os << "sn("; _G.print(os); os <<")"; }
};

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
    enum {Space = 0, ScalarValued= 0, ColBlocks= 0};

    tangent_expr(const gsGeometryMap<T> & G) : _G(G) { }

    mutable gsVector<Scalar> res;
    const gsVector<Scalar> & eval(const index_t k) const
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
            res.resize(3);
            res.col3d(0) = _G.data().normals.col3d(k)
                .cross( _G.data().outNormals.col3d(k) );
            return res;
        }
        else
            GISMO_ERROR("Function not implemented for dimension"<<_G.targetDim());

    }

    index_t rows() const { return _G.data().dim.second; }
    index_t cols() const { return 1; }

    static const gsFeSpace<Scalar> & rowVar() {return gsNullExpr<Scalar>::get();}
    static const gsFeSpace<Scalar> & colVar() {return gsNullExpr<Scalar>::get();}

    void parse(gsExprHelper<Scalar> & evList) const
    {
        evList.add(_G);
        _G.data().flags |= NEED_NORMAL;
        _G.data().flags |= NEED_OUTER_NORMAL;
    }

    void print(std::ostream &os) const { os << "tv("; _G.print(os); os <<")"; }
};

/**
   Expression for the Laplacian of a finite element variable
*/
template<class E>
class lapl_expr : public _expr<lapl_expr<E> >
{
    typename E::Nested_t _u;

public:
    typedef typename E::Scalar Scalar;
    enum {Space = E::Space, ScalarValued= 0, ColBlocks= 0};

    lapl_expr(const E & u) : _u(u) { }

    auto eval(const index_t k) const -> decltype(_u.data().laplacians.col(k))
    {
        // numActive x 1
        return _u.data().laplacians.col(k);
        //todo: replace by
        // NEED_DERIV2
        // ..nabla2.sum()
    }

    index_t rows() const { return _u.data().laplacians.rows(); }
    index_t cols() const { return 1; }

    index_t cardinality_impl() const { return _u.cardinality_impl(); }

    void parse(gsExprHelper<Scalar> & evList) const
    {
        evList.add(_u);
        _u.data().flags |= NEED_LAPLACIAN;
    }

    static const gsFeSpace<Scalar> & rowVar() {return E::rowVar();}
    static const gsFeSpace<Scalar> & colVar() {return gsNullExpr<Scalar>::get();}

    void print(std::ostream &os) const { os << "\u2206("; _u.print(os); os <<")"; } //or \u0394
};

/*
  Expression for the Laplacian of a finite element solution
*/
template<class T>
class lapl_expr<gsFeSolution<T> > : public _expr<lapl_expr<gsFeSolution<T> > >
{
protected:
    const gsFeSolution<T> _u;

public:
    typedef T Scalar;
    enum {Space = 0, ScalarValued= 0, ColBlocks= 0};

    lapl_expr(const gsFeSolution<T> & u) : _u(u) { }

    mutable gsMatrix<T> res;
    const gsMatrix<T> & eval(const index_t k) const
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

    void parse(gsExprHelper<Scalar> & evList) const
    {
        evList.add(_u.space());
        _u.data().flags |= NEED_ACTIVE | NEED_DERIV2;
    }

    const gsFeSpace<Scalar> & rowVar() const {return gsNullExpr<T>::get();}
    const gsFeSpace<Scalar> & colVar() const {return gsNullExpr<T>::get();}

    void print(std::ostream &os) const { os << "\u2206(s)"; }
};

/*
  Expression for the (precomputed) second fundamental form of a surface
*/
template<class T>
class fform2nd_expr  : public _expr<fform2nd_expr<T> >
{
    typename gsGeometryMap<T>::Nested_t _G;
public:
    typedef T Scalar;
    enum {Space = 0, ScalarValued = 0, ColBlocks = 0};

    fform2nd_expr(const gsGeometryMap<T> & G) : _G(G) { }

    const gsAsConstMatrix<Scalar> eval(const index_t k) const
    {
        return gsAsConstMatrix<Scalar>(_G.data().fundForms.col(k).data(),rows(),cols());
    }

    index_t rows() const { return _G.data().dim.first ; }
    index_t cols() const { return _G.data().dim.first ; }

    void parse(gsExprHelper<Scalar> & evList) const
    {
        evList.add(_G);
        _G.data().flags |= NEED_2ND_FFORM;
    }

    const gsFeSpace<Scalar> & rowVar() const {return gsNullExpr<T>::get();}
    const gsFeSpace<Scalar> & colVar() const {return gsNullExpr<T>::get();}

    void print(std::ostream &os) const { os << "fform2nd("; _G.print(os); os <<")"; }
};

/*
  Expression for the (precomputed) inverse of the Jacobian matrix of
  a geometry map
*/
template<class T>
class jacInv_expr  : public _expr<jacInv_expr<T> >
{
    typename gsGeometryMap<T>::Nested_t _G;
public:
    typedef T Scalar;
    enum {Space = 0, ScalarValued = 0, ColBlocks = 0};

    jacInv_expr(const gsGeometryMap<T> & G) : _G(G)
    {
        // Note: for non-square Jacobian matrices, generalized inverse, i.e.: (J^t J)^{-t} J^t
        //GISMO_ASSERT(rows() == cols(), "The Jacobian matrix is not square");
    }

    MatExprType eval(const index_t k) const { return _G.data().jacInvTr.reshapeCol(k,cols(),rows()).transpose(); }

    index_t rows() const { return _G.source().domainDim(); }
    index_t cols() const { return _G.source().targetDim(); }

    void parse(gsExprHelper<Scalar> & evList) const
    {
        evList.add(_G);
        _G.data().flags |= NEED_GRAD_TRANSFORM;
    }

    const gsFeSpace<Scalar> & rowVar() const {return gsNullExpr<T>::get();}
    const gsFeSpace<Scalar> & colVar() const {return gsNullExpr<T>::get();}

    // todo mat_expr ?
    // tr() const --> _G.data().fundForm(k)

    void print(std::ostream &os) const { os << "jacInv("; _G.print(os); os <<")"; }
};

/*
  Expression for the Jacobian matrix of a FE variable
*/
template<class E>
class jac_expr : public _expr<jac_expr<E> >
{
    typename E::Nested_t _u;
public:
    enum {ColBlocks = (1==E::Space?1:0) };
    enum {Space = E::Space, ScalarValued= 0 };

    typedef typename E::Scalar Scalar;

    mutable gsMatrix<Scalar> res;

    jac_expr(const E & _u_)
    : _u(_u_) { }

    MatExprType eval(const index_t k) const
    {
        if (0!=Space)
        {
            // Dim x (numActive*Dim)
            res = _u.data().values[1].col(k).transpose().blockDiag(_u.dim());
        }
        else
        {
            res = _u.data().values[1]
                .reshapeCol(k, _u.parDim(), _u.targetDim()).transpose()
                .blockDiag(_u.dim());
        }
        return res;
    }

    const gsFeSpace<Scalar> & rowVar() const { return _u.rowVar(); }
    //const gsFeSpace<Scalar> & rowVar() const { return rowVar_impl<E>(); }
    const gsFeSpace<Scalar> & colVar() const { return gsNullExpr<Scalar>::get(); }

    index_t rows() const { return rows_impl(_u); }
    index_t cols() const { return cols_impl(_u); }

    // index_t rows() const { return _u.dim(); }
    // index_t cols() const { return _u.source().domainDim(); }

    index_t cardinality_impl() const
    {
        return _u.dim() * _u.data().actives.rows();
    }

    void parse(gsExprHelper<Scalar> & evList) const
    {
        evList.add(_u);
        _u.data().flags |= NEED_DERIV|NEED_ACTIVE;
        //note: cardinality() depends on actives
    }

    void print(std::ostream &os) const { os << "\u2207("; _u.print(os);os <<")"; }

private:

    // The jacobian is different for gsFeVariable, gsFeSolution and gsFeSpace
    // gsFeSolution: Does not work
    // gsFeVariable: dim()=1 and source().targetDim()=d
    // gsFeSpace: dim()=d and source().targetDim()=1
    template<class U> inline
    typename util::enable_if<!(util::is_same<U,gsFeSpace<Scalar> >::value), index_t >::type  // What about solution??
    rows_impl(const U & u)  const
    {
        return u.source().targetDim();
    }

    template<class U> inline
    typename util::enable_if< (util::is_same<U,gsFeSpace<Scalar> >::value), index_t >::type
    rows_impl(const U & u) const
    {
        return u.dim();
    }

    template<class U> inline
    typename util::enable_if<!(util::is_same<U,gsFeSpace<Scalar> >::value), index_t >::type
    cols_impl(const U & u)  const
    {
        return u.source().domainDim();
    }

    template<class U> inline
    typename util::enable_if< (util::is_same<U,gsFeSpace<Scalar> >::value), index_t >::type
    cols_impl(const U & u) const
    {
        return u.source().domainDim();
    }

    // The jacobian is different for gsFeVariable, gsFeSolution and gsFeSpace
    // gsFeSolution: Does not work
    // gsFeVariable: rowVar = NULL
    // gsFeSpace: rowVar = u
    template<class U> inline
    typename util::enable_if<!(util::is_same<U,gsFeSpace<Scalar> >::value), const gsFeSpace<Scalar> &  >::type
    rowVar_impl() const
    {
        return gsNullExpr<Scalar>::get();
    }

    template<class U> inline
    typename util::enable_if<(util::is_same<U,gsFeSpace<Scalar> >::value), const gsFeSpace<Scalar> &  >::type
    rowVar_impl() const
    {
        return _u;
    }
};

/*
  Expression for the Jacobian matrix of a geometry map
*/
template<class T>
class jac_expr<gsGeometryMap<T> > : public _expr<jac_expr<gsGeometryMap<T> > >
{
    typename gsGeometryMap<T>::Nested_t _G;

public:
    typedef T Scalar;
    enum {Space = 0, ScalarValued= 0, ColBlocks= 0};

    jac_expr(const gsGeometryMap<T> & G) : _G(G) { }
    MatExprType eval(const index_t k) const
    {
        // TarDim x ParDim
        return _G.data().values[1]
            .reshapeCol(k, _G.data().dim.first, _G.data().dim.second).transpose();
    }

    index_t rows() const { return _G.source().targetDim(); }

    index_t cols() const { return _G.source().domainDim(); }

    static const gsFeSpace<Scalar> & rowVar() { return gsNullExpr<Scalar>::get(); }
    static const gsFeSpace<Scalar> & colVar() { return gsNullExpr<Scalar>::get(); }

    void parse(gsExprHelper<Scalar> & evList) const
    {
        evList.add(_G);
        _G.data().flags |= NEED_DERIV;
    }

    meas_expr<T> absDet() const
    {
        GISMO_ASSERT(rows() == cols(), "The Jacobian matrix is not square");
        return meas_expr<T>(_G);
    }

    jacInv_expr<T> inv() const
    {
        GISMO_ASSERT(rows() == cols(), "The Jacobian matrix is not square");
        return jacInv_expr<T>(_G);
    }

    /// The generalized Jacobian matrix inverse, i.e.: (J^t J)^{-t} J^t
    jacInv_expr<T> ginv() const { return jacInv_expr<T>(_G); }

    void print(std::ostream &os) const { os << "\u2207("; _G.print(os); os <<")"; }
};

template<class E>
class hess_expr : public _expr<hess_expr<E> >
{
public:
    typedef typename E::Scalar Scalar;
private:
    typename E::Nested_t _u;
    mutable gsMatrix<Scalar> res;
public:
    enum {ScalarValued = 0, ColBlocks = (1==E::Space?1:0) };
    enum {Space = E::Space };

public:
    hess_expr(const E & u) : _u(u)
    {
        //gsInfo << "\n-expression is space ? "<<E::Space <<"\n"; _u.print(gsInfo);
        //GISMO_ASSERT(1==_u.dim(),"hess(.) requires 1D variable");
    }

    const gsMatrix<Scalar> & eval(const index_t k) const
    {
        const gsFuncData<Scalar> & dd = _u.data();
        const index_t sz = cardinality_impl();
        res.resize(dd.dim.first, sz*dd.dim.first);
        secDerToHessian(dd.values[2].col(k), dd.dim.first, res);
        res.resize(dd.dim.first, res.cols()*dd.dim.first);
        // Note: auto returns by value here,
        // in C++11 we may add in -> decltype(res) &
        return res;
    }

    index_t rows() const { return _u.data().dim.first; }
    index_t cols() const
    {   return rows();
        //return 2*_u.data().values[2].rows() / (1+_u.data().dim.first);
    }

    index_t cardinality_impl() const
    {
        return 2*_u.data().values[2].rows()/
            (_u.data().dim.first*(1+_u.data().dim.first));
        //gsDebugVar(_u.data().values.front().rows());//empty!
    }

    void parse(gsExprHelper<Scalar> & evList) const
    {
        evList.add(_u);
        _u.data().flags |= NEED_2ND_DER;
    }

    const gsFeSpace<Scalar> & rowVar() const { return _u.rowVar(); }
    const gsFeSpace<Scalar> & colVar() const { return gsNullExpr<Scalar>::get(); }

    void print(std::ostream &os) const
    //    { os << "hess("; _u.print(os);os <<")"; }
    { os << "\u210D(U)"; }
};

template<class T>
class hess_expr<gsFeSolution<T> > : public _expr<hess_expr<gsFeSolution<T> > >
{
protected:
    const gsFeSolution<T> _u;

public:
    typedef T Scalar;
    enum{Space = 0, ScalarValued = 0, ColBlocks = 0 };

    hess_expr(const gsFeSolution<T> & u) : _u(u) { }

    mutable gsMatrix<T> res;
    const gsMatrix<T> & eval(const index_t k) const
    {
        GISMO_ASSERT(1==_u.data().actives.cols(), "Single actives expected. Actives: \n"<<_u.data().actives);

        const gsDofMapper & map = _u.mapper();

        const index_t numActs = _u.data().values[0].rows();
        const index_t pdim = _u.parDim();
        index_t numDers = pdim*(pdim+1)/2;
        gsMatrix<T> deriv2;

        // In the scalar case, the hessian is returned as a pdim x pdim matrix
        if (1==_u.dim())
        {
            res.setZero(numDers,1);
            for (index_t i = 0; i!=numActs; ++i)
            {
                const index_t ii = map.index(_u.data().actives.at(i), _u.data().patchId,0);
                deriv2 = _u.data().values[2].block(i*numDers,k,numDers,1);
                if ( map.is_free_index(ii) ) // DoF value is in the solVector
                    res += _u.coefs().at(ii) * deriv2;
                else
                    res +=_u.fixedPart().at( map.global_to_bindex(ii) ) * deriv2;
            }
            secDerToHessian(res, pdim, deriv2);
            res.swap(deriv2);
            res.resize(pdim,pdim);
        }
        // In the vector case, the hessian is returned as a matrix where each row corresponds to the component of the solution and contains the derivatives in the columns
        else
        {
            res.setZero(rows(), numDers);
            for (index_t c = 0; c != _u.dim(); c++)
                for (index_t i = 0; i != numActs; ++i) {
                    const index_t ii = map.index(_u.data().actives.at(i), _u.data().patchId, c);
                    deriv2 = _u.space().data().values[2].block(i * numDers, k, numDers,
                                                               1).transpose(); // start row, start col, rows, cols
                    if (map.is_free_index(ii)) // DoF value is in the solVector
                        res.row(c) += _u.coefs().at(ii) * deriv2;
                    else
                        res.row(c) += _u.fixedPart().at(map.global_to_bindex(ii)) * deriv2;
                }
        }
        return res;
    }

    index_t rows() const
    {
        if (1==_u.dim())
            return _u.parDim();
        else
            return _u.dim(); //  number of components
    }
    index_t cols() const
    {
        if (1==_u.dim())
            return _u.parDim();
        // second derivatives in the columns; i.e. [d11, d22, d33, d12, d13, d23]
        else
            return _u.parDim() * (_u.parDim() + 1) / 2;
    }

    const gsFeSpace<Scalar> & rowVar() const { return gsNullExpr<Scalar>::get(); }
    const gsFeSpace<Scalar> & colVar() const { return gsNullExpr<Scalar>::get(); }

    void parse(gsExprHelper<Scalar> & evList) const
    {
        _u.parse(evList);                         // add symbol
        evList.add(_u.space());
        _u.data().flags |= NEED_ACTIVE | NEED_VALUE | NEED_DERIV2;
    }

    void print(std::ostream &os) const { os << "\u210D(s)"; }
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

    void parse(gsExprHelper<Scalar> & evList) const
    {
        evList.add(_G);
        _G.data().flags |= NEED_2ND_DER;
    }
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
    enum{ Space = 1, ScalarValued= 0, ColBlocks= 0};

    curl_expr(const gsFeVariable<T> & u) : _u(u)
    { GISMO_ASSERT(3==u.dim(),"curl(.) requires 3D variable."); }

    const gsMatrix<T> & eval(const index_t k) const
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

    void parse(gsExprHelper<Scalar> & evList) const
    {
        evList.add(_u);
        _u.data().flags |= NEED_GRAD;
    }

    const gsFeSpace<T> & rowVar() const { return _u.rowVar(); }
    const gsFeSpace<T> & colVar() const {return gsNullExpr<T>::get();}

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
    enum {Space = (int)E1::Space + (int)E2::Space };

    typedef typename E1::Scalar Scalar;

    mult_expr(_expr<E1> const& u,
              _expr<E2> const& v)
    : _u(u), _v(v) { }

    mutable Temporary_t tmp;
    const Temporary_t & eval(const index_t k) const
    {
        GISMO_ASSERT(0==_u.cols()*_v.rows() || _u.cols() == _v.rows(),
                     "Wrong dimensions "<<_u.cols()<<"!="<<_v.rows()<<" in * operation:\n"
                     << _u <<" times \n" << _v );

        // Note: a * b * c --> (a*b).eval()*c
        tmp = _u.eval(k) * _v.eval(k);
        return tmp; // assumes result is not scalarvalued
    }

    index_t rows() const { return E1::ScalarValued ? _v.rows()  : _u.rows(); }
    index_t cols() const { return E2::ScalarValued ? _u.cols()  : _v.cols(); }
    void parse(gsExprHelper<Scalar> & evList) const
    { _u.parse(evList); _v.parse(evList); }


    index_t cardinality_impl() const
    { return 0==E1::Space ? _v.cardinality(): _u.cardinality(); }

    const gsFeSpace<Scalar> & rowVar() const
    { return 0==E1::Space ? _v.rowVar() : _u.rowVar(); }
    const gsFeSpace<Scalar> & colVar() const
    { return 0==E2::Space ? _u.colVar() : _v.colVar(); }

    void print(std::ostream &os) const { _u.print(os); os<<"*"; _v.print(os); }
};

/*
  Expression for multiplication operation (second version)

  First argument E1 has ColBlocks = true

  Partial specialization for (right) blockwise multiplication
  [A1 A2 A3] * B = [A1*B  A2*B  A3*B]

  as well as

  both are ColBlocks: [A1 A2 A3] * [B1 B2 B3] = [A1*B1  A2*B2  A3*B3]
                                                [A2*B1 ..           ]
                                                [                   ]
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
    enum {ScalarValued = 0, ColBlocks = E1::ColBlocks}; //(!)
    enum {Space = (int)E1::Space + (int)E2::Space };

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
        const auto tmpA = _u.eval(k);
        const auto tmpB = _v.eval(k);

        const index_t vc = _v.cols();

        // either _v.cardinality()==1 or _v.cardinality()==_u.cardinality()
        if (  1 == _v.cardinality() ) //second is not ColBlocks
        {
            res.resize(ur, vc*nb);
            GISMO_ASSERT(tmpA.cols()==uc*nb, "Dimension error.. "<< tmpA.cols()<<"!="<<uc*nb );
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

    void parse(gsExprHelper<Scalar> & evList) const
    { _u.parse(evList); _v.parse(evList); }


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

    void parse(gsExprHelper<Scalar> & evList) const
    { _v.parse(evList); }

    index_t cardinality_impl() const
    { return _v.cardinality(); }

    const gsFeSpace<Scalar> & rowVar() const { return _v.rowVar(); }
    const gsFeSpace<Scalar> & colVar() const { return _v.colVar(); }

    void print(std::ostream &os) const { os << _c <<"*";_v.print(os); }
};

/*
  Expression for multiplication operation (fourth version)

  Multiplication by gsMatrix
*/
template <typename E2, bool rmult>
class mmult_expr
    : public _expr<mmult_expr<E2,rmult> >
{
public:
    typedef typename E2::Scalar Scalar;
private:
    const gsMatrix<Scalar> & _c;
    typename E2::Nested_t _v;

public:
    enum {ScalarValued = 0, ColBlocks = E2::ColBlocks};
    enum {Space = E2::Space};

    mmult_expr(gsMatrix<Scalar> const & c, _expr<E2> const& v)
    : _c(c), _v(v) { }

    EIGEN_STRONG_INLINE AutoReturn_t eval(const index_t k) const
    {
        return rmult ? ( _v.eval(k) * _c ) : ( _c * _v.eval(k) );
    }

    index_t rows() const { return _c.rows(); }
    index_t cols() const { return _v.cols(); }

    void parse(gsExprHelper<Scalar> & evList) const
    { _v.parse(evList); }

    index_t cardinality_impl() const
    { return _v.cardinality(); }

    const gsFeSpace<Scalar> & rowVar() const { return _v.rowVar(); }
    const gsFeSpace<Scalar> & colVar() const { return _v.colVar(); }

    void print(std::ostream &os) const { os << "M*";_v.print(os); }
};

template <typename E1, typename E2>
class collapse_expr : public _expr<collapse_expr<E1, E2> >
{
    typename E1::Nested_t _u;
    typename E2::Nested_t _v;

public:
    enum {ScalarValued = 0, ColBlocks = 0};
    enum { Space = (int)E1::Space + (int)E2::Space };

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
        const auto tmpA = _u.eval(k);
        const auto tmpB = _v.eval(k);

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

    void parse(gsExprHelper<Scalar> & evList) const
    { _u.parse(evList); _v.parse(evList); }

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

  [A1 A2 A3] . [B1 B2 B3]
  =
  [ A1.B1  A1.B2  A1.B3 ]
  [ A2.B1  A2.B2  A2.B3 ]
  [ A3.B1  A3.B2  A3.B3 ]
*/
template <typename E1, typename E2, bool = E2::ColBlocks>
class frprod_expr : public _expr<frprod_expr<E1, E2> >
{
public:
    typedef typename E1::Scalar Scalar;
    enum {ScalarValued = 0, ColBlocks=E2::ColBlocks};
    enum { Space = (int)E1::Space + (int)E2::Space };
    // E1 E2 this (16 cases..)
    // 0  0  0
    // 1  1
    // 2  2
    // 3  3

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
        const index_t rb = _u.rows();
        const index_t nb = _u.cardinality();
        auto A = _u.eval(k);
        auto B = _v.eval(k);
        res.resize(nb, nb);
        for (index_t i = 0; i!=nb; ++i) // all with all
            for (index_t j = 0; j!=nb; ++j)
                res(i,j) =
                    (A.middleCols(i*rb,rb).array() * B.middleCols(j*rb,rb).array()).sum();
        return res;
    }

    index_t rows() const { return _u.cols() / _u.rows(); }
    index_t cols() const { return _u.cols() / _u.rows(); }

    void parse(gsExprHelper<Scalar> & evList) const
    { _u.parse(evList); _v.parse(evList); }

    const gsFeSpace<Scalar> & rowVar() const { return _u.rowVar(); }
    const gsFeSpace<Scalar> & colVar() const { return _v.colVar(); }

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
    enum {ScalarValued = 0, Space = E1::Space, ColBlocks= E1::ColBlocks};

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
        auto A = _u.eval(k);
        auto B = _v.eval(k);
        const index_t rb = A.rows(); //==cb
        const index_t nb = _u.cardinality();
        res.resize(nb, 1);
        for (index_t i = 0; i!=nb; ++i) // all with all
            res(i,0) =
                (A.middleCols(i*rb,rb).array() * B.array()).sum();
        return res;
    }

    index_t rows() const { return _u.cols() / _u.rows(); }
    index_t cols() const { return 1; }

    void parse(gsExprHelper<Scalar> & evList) const
    { _u.parse(evList); _v.parse(evList); }


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
    enum {ScalarValued = E1::ScalarValued, ColBlocks= E2::ColBlocks};
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

    void parse(gsExprHelper<Scalar> & evList) const
    { _u.parse(evList); _v.parse(evList); }


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
    enum {Space= E1::Space, ScalarValued = E1::ScalarValued, ColBlocks= E1::ColBlocks};

    divide_expr(_expr<E1> const& u, Scalar const  c)
    : _u(u), _c(c) { }

    AutoReturn_t eval(const index_t k) const
    { return ( _u.eval(k) / _c ); }

    index_t rows() const { return _u.rows(); }
    index_t cols() const { return _u.cols(); }

    void parse(gsExprHelper<Scalar> & evList) const
    { _u.parse(evList); }


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
    enum {Space= 0, ScalarValued = 1, ColBlocks= 0};

    divide_expr(Scalar const c, _expr<E2> const& u)
    : _c(c), _u(u)
    { GISMO_STATIC_ASSERT(E2::ScalarValued, "The denominator needs to be scalar valued."); }

    Scalar eval(const index_t k) const
    { return ( _c / _u.val().eval(k) ); }

    index_t rows() const { return 0; }
    index_t cols() const { return 0; }

    void parse(gsExprHelper<Scalar> & evList) const
    { _u.parse(evList); }


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
        ColBlocks = E1::ColBlocks || E2::ColBlocks };
    enum {Space = E1::Space}; // == E2::Space

    typedef typename E1::Scalar Scalar;

    add_expr(_expr<E1> const& u, _expr<E2> const& v)
    : _u(u), _v(v)
    {
        GISMO_ENSURE((int)E1::Space == (int)E2::Space &&
                     _u.rowVar()==_v.rowVar() && _u.colVar()==_v.colVar(),
                     "Error: adding apples and oranges (use comma instead),"
                     " namely:\n" << _u <<"\n"<<_v<<
                     " \nvars:\n" << _u.rowVar().id()<<"!="<<_v.rowVar().id() <<", "<< _u.colVar().id()<<"!="<<_v.colVar().id()<<
                     " \nspaces:\n" << (int)E1::Space<< "!="<< (int)E2::Space
            );
    }

    mutable Temporary_t res;
    const Temporary_t & eval(const index_t k) const
    {
        GISMO_ASSERT(_u.rows() == _v.rows(),
                     "Wrong dimensions "<<_u.rows()<<"!="<<_v.rows()<<" in + operation:\n"
                     << _u <<" plus \n" << _v );
        GISMO_ASSERT(_u.cols() == _v.cols(),
                     "Wrong dimensions "<<_u.cols()<<"!="<<_v.cols()<<" in + operation:\n"
                     << _u <<" plus \n" << _v );
        res = _u.eval(k) + _v.eval(k);
        return res;
    }

    index_t rows() const { return _u.rows(); }
    index_t cols() const { return _u.cols(); }

    void parse(gsExprHelper<Scalar> & evList) const
    { _u.parse(evList); _v.parse(evList); }


    index_t cardinality_impl() const { return _u.cardinality_impl(); }

    const gsFeSpace<Scalar> & rowVar() const { return _u.rowVar(); }
    const gsFeSpace<Scalar> & colVar() const { return _v.colVar(); }

    void print(std::ostream &os) const
    { os << "("; _u.print(os);os <<" + ";_v.print(os);os << ")"; }
};

/*
  Expression for addition operation between an expression and a scalar
*/
template <typename E1>
class add_expr<typename E1::Scalar,E1>
    : public _expr<add_expr<typename E1::Scalar,E1> >
{
    typename E1::Nested_t _u;
    typename E1::Scalar   _v; 

public:
    enum {ScalarValued = 1, ColBlocks = 0};
    enum {Space = E1::Space}; // == E2::Space

    typedef typename E1::Scalar Scalar;

    add_expr(_expr<E1> const& u, Scalar v)
    : _u(u), _v(v)
    {
        GISMO_ENSURE( E1::ScalarValued && !E1::ColBlocks, 
        "Error, addition between expressions and numbers works only if the expression is scalar valued.[ try using: ("<<_u<<" ).val() ]"
        );
    }

    Scalar eval(const index_t k) const
    {
        return _u.eval(k) + _v;
    }

    index_t rows() const { return 1; }
    index_t cols() const { return 1; }

    void parse(gsExprHelper<Scalar> & evList) const
    { _u.parse(evList);}


    index_t cardinality_impl() const { return _u.cardinality_impl(); }

    const gsFeSpace<Scalar> & rowVar() const { return _u.rowVar(); }
    const gsFeSpace<Scalar> & colVar() const { return gsNullExpr<Scalar>::get(); }

    void print(std::ostream &os) const
    { os << "("; _u.print(os);os <<" + "<< _v << ")"; }
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

    enum {Space = E1::Space, ScalarValued= 0, ColBlocks= E2::ColBlocks};

    summ_expr(E1 const& u, E2 const& M) : _u(u), _M(M) { }

    const gsMatrix<Scalar> & eval(const index_t k) const
    {
        auto sl   = _u.eval(k);
        const index_t sr = sl.rows();
        auto ml   = _M.eval(k);
        const index_t mr = ml.rows();
        const index_t mb = _M.cardinality();

        GISMO_ASSERT(_M.cols()==_M.rows(),"Matrix must be square: "<< _M.rows()<<" x "<< _M.cols() << " expr: "<< _M );
        GISMO_ASSERT(mb==_u.cols(),"cardinality must match vector, but card(M)="<<_M.cardinality()<<" and cols(u)="<<_u.cols());

        res.setZero(mr, sr * mr);
        for (index_t i = 0; i!=sr; ++i)
            for (index_t t = 0; t!=mb; ++t) // lc
                res.middleCols(i*mr,mr) += sl(i,t) * ml.middleCols(t*mr,mr);
        return res;
    }

    index_t rows() const { return _M.rows(); }
    index_t cols() const { return _M.rows(); } //_u.rows()

    void parse(gsExprHelper<Scalar> & evList) const
    { _u.parse(evList); _M.parse(evList); }

    const gsFeSpace<Scalar> & rowVar() const { return _u.rowVar(); }
    const gsFeSpace<Scalar> & colVar() const { return gsNullExpr<Scalar>::get(); }

    index_t cardinality_impl() const
    { GISMO_ERROR("Something went terribly wrong"); }

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
        ColBlocks = E1::ColBlocks || E2::ColBlocks };
    enum {Space = E1::Space}; // == E2::Space

    typedef typename E1::Scalar Scalar;

    sub_expr(_expr<E1> const& u, _expr<E2> const& v)
    : _u(u), _v(v)
    {
        GISMO_ENSURE((int)E1::Space == (int)E2::Space &&
                     _u.rowVar()==_v.rowVar() && _u.colVar()==_v.colVar(),
                     "Error: substracting apples from oranges (use comma instead),"
                     " namely:\n" << _u <<"\n"<<_v);
    }

    mutable Temporary_t res;
    const Temporary_t & eval(const index_t k) const
    {
        // GISMO_ASSERT(_u.rowVar().id()==_v.rowVar().id() && _u.rowVar().isAcross()==_v.rowVar().isAcross(),
        //     "The row spaces are not split compatibly.");
        // GISMO_ASSERT(_u.colVar().id()==_v.colVar().id() && _u.colVar().isAcross()==_v.colVar().isAcross(),
        //     "The col spaces are not split compatibly.");
        GISMO_ASSERT(_u.rows() == _v.rows(),
                     "Wrong dimensions "<<_u.rows()<<"!="<<_v.rows()<<" in - operation:\n" << _u <<" minus \n" << _v );
        GISMO_ASSERT(_u.cols() == _v.cols(),
                     "Wrong dimensions "<<_u.cols()<<"!="<<_v.cols()<<" in - operation:\n" << _u <<" minus \n" << _v );
        GISMO_ASSERT(_u.cardinality() == _u.cardinality(),
                     "Cardinality "<< _u.cardinality()<<" != "<< _v.cardinality());
        //return (_u.eval(k) - _v.eval(k) ).eval();
        //return (_u.eval(k) - _v.eval(k) ); // any temporary matrices eval(.) will leak mem.
        res = _u.eval(k) - _v.eval(k);
        return res;
    }

    index_t rows() const { return _u.rows(); }
    index_t cols() const { return _u.cols(); }

    void parse(gsExprHelper<Scalar> & evList) const
    { _u.parse(evList); _v.parse(evList); }

    const gsFeSpace<Scalar> & rowVar() const { return _u.rowVar(); }
    const gsFeSpace<Scalar> & colVar() const { return _v.colVar(); }

    index_t cardinality_impl() const
    {
        GISMO_ASSERT(_u.cardinality() == _u.cardinality(),
                     "Cardinality "<< _u.cardinality()<<" != "<< _v.cardinality());
        return _u.cardinality();
    }

    void print(std::ostream &os) const
    { os << "("; _u.print(os); os<<" - ";_v.print(os); os << ")";}
};

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

    enum { Space = (0==E::Space ? 0 : E::Space), ScalarValued= E::ScalarValued, ColBlocks= E::ColBlocks };

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

    void parse(gsExprHelper<Scalar> & evList) const
    { _u.parse(evList); }

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
    enum { Space = (0==E::Space ? 0 : 3), ScalarValued=E::ScalarValued, ColBlocks= E::ColBlocks };
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

    void parse(gsExprHelper<Scalar> & evList) const
    { _u.parse(evList); }

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

/// The identity matrix of dimension \a dim
EIGEN_STRONG_INLINE idMat_expr id(const index_t dim) { return idMat_expr(dim); }

// Returns the unit as an expression
//EIGEN_STRONG_INLINE _expr<real_t> one() { return _expr<real_t,true>(1); }


/// Absolute value
template<class E> EIGEN_STRONG_INLINE
abs_expr<E> abs(const E & u) { return abs_expr<E>(u); }

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

template<class E> EIGEN_STRONG_INLINE
lapl_expr<E> lapl(const symbol_expr<E> & u) { return lapl_expr<E>(u); }

template<class T> EIGEN_STRONG_INLINE
lapl_expr<gsFeSolution<T> > lapl(const gsFeSolution<T> & u)
{ return lapl_expr<gsFeSolution<T> >(u); }

/// The second fundamental form of \a G
template<class T> EIGEN_STRONG_INLINE fform2nd_expr<T> fform2nd(const gsGeometryMap<T> & G)
{ return fform2nd_expr<T>(G); }

/// The Jacobian matrix of a FE variable
template<class E> EIGEN_STRONG_INLINE
jac_expr<E> jac(const symbol_expr<E> & u) { return jac_expr<E>(u); }

/// The Jacobian matrix of a geometry map
template<class T> EIGEN_STRONG_INLINE
jac_expr<gsGeometryMap<T> > jac(const gsGeometryMap<T> & G) {return jac_expr<gsGeometryMap<T> >(G);}

/// Jacobian matrix for a solution expression
template<class T> EIGEN_STRONG_INLINE
grad_expr<gsFeSolution<T> > jac(const gsFeSolution<T> & s) {return grad_expr<gsFeSolution<T> >(s);}

template<class E> EIGEN_STRONG_INLINE
hess_expr<E> hess(const symbol_expr<E> & u) { return hess_expr<E>(u); }

/// The hessian of a geometry map
template<class T> EIGEN_STRONG_INLINE
hess_expr<gsGeometryMap<T> > hess(const gsGeometryMap<T> & u) { return hess_expr<gsGeometryMap<T> >(u); }

/// The hessian of a solution variable
template<class T> EIGEN_STRONG_INLINE
hess_expr<gsFeSolution<T> > hess(const gsFeSolution<T> & u) { return hess_expr<gsFeSolution<T> >(u); }

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

template <typename E1> EIGEN_STRONG_INLINE
mult_expr<typename E1::Scalar,E1,false> const
operator*(_expr<E1> const& v, typename E1::Scalar const& u)
{ return mult_expr<typename E1::Scalar,E1, false>(u, v); }

template <typename E1> EIGEN_STRONG_INLINE
mult_expr<typename E1::Scalar,E1,false> const
operator-(_expr<E1> const& u)
{ return mult_expr<typename E1::Scalar,E1, false>(-1, u); }

template <typename E2> mmult_expr<E2,false> const
operator*(gsMatrix<typename E2::Scalar> const& u, _expr<E2> const& v)
{ return mmult_expr<E2, false>(u, v); }

template <typename E1> mmult_expr<E1,true> const
operator*(_expr<E1> const& v, gsMatrix<typename E1::Scalar> const& u)
{ return mmult_expr<E1, true>(u, v); }

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

/// Addition operator for expressions and numbers
template <typename E> EIGEN_STRONG_INLINE
add_expr<typename E::Scalar, E> const
operator+(_expr<E> const& u, const typename E::Scalar v)
{ return add_expr<typename E::Scalar, E>(u, v); }

/// Addition operator for expressions and numbers
template <typename E> EIGEN_STRONG_INLINE
add_expr<typename E::Scalar, E> const
operator+(const typename E::Scalar v, _expr<E> const& u)
{ return add_expr<typename E::Scalar, E>(u, v); }

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


// Shortcuts for common quantities, for instance function
// transformations by the geometry map \a G
#define GISMO_SHORTCUT_VAR_EXPRESSION(name,impl) template<class E> EIGEN_STRONG_INLINE \
    auto name(const E & u) -> decltype(impl) { return impl; }
#define GISMO_SHORTCUT_MAP_EXPRESSION(name,impl) template<class T> EIGEN_STRONG_INLINE \
    auto name(const gsGeometryMap<T> & G)  -> decltype(impl) { return impl; }
#define GISMO_SHORTCUT_PHY_EXPRESSION(name,impl) template<class E> EIGEN_STRONG_INLINE \
    auto name(const E & u, const gsGeometryMap<typename E::Scalar> & G)  -> decltype(impl) { return impl; }

// Divergence
GISMO_SHORTCUT_VAR_EXPRESSION(  div, jac(u).trace() )
GISMO_SHORTCUT_PHY_EXPRESSION( idiv, ijac(u,G).trace()    )

// The unit (normalized) boundary (outer pointing) normal
GISMO_SHORTCUT_MAP_EXPRESSION(unv, nv(G).normalized()   )
// The unit (normalized) boundary (surface) normal
GISMO_SHORTCUT_MAP_EXPRESSION(usn, sn(G).normalized()   )

GISMO_SHORTCUT_PHY_EXPRESSION(igrad, grad(u)*jac(G).ginv() ) // transpose() problem ??
GISMO_SHORTCUT_VAR_EXPRESSION(igrad, grad(u) ) // u is presumed to be defined over G

GISMO_SHORTCUT_PHY_EXPRESSION( ijac, jac(u) * jac(G).ginv())

// note and todo: does this work for non-scalar solutions?
GISMO_SHORTCUT_PHY_EXPRESSION(ihess,
                              jac(G).ginv().tr()*( hess(u) - summ(igrad(u,G),hess(G)) ) * jac(G).ginv() )
GISMO_SHORTCUT_VAR_EXPRESSION(ihess, hess(u) )

GISMO_SHORTCUT_PHY_EXPRESSION(ilapl, ihess(u,G).trace()   )
GISMO_SHORTCUT_VAR_EXPRESSION(ilapl, hess(u).trace() )

GISMO_SHORTCUT_VAR_EXPRESSION(fform, jac(u).tr()*jac(u) )

#undef GISMO_SHORTCUT_PHY_EXPRESSION
#undef GISMO_SHORTCUT_VAR_EXPRESSION
#undef GISMO_SHORTCUT_MAP_EXPRESSION

} // namespace expr

} //namespace gismo
