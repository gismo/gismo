/** @file embeddedCurve_example.cpp

    @brief Example for embedded curve

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): C. Chianese, H.M. Verhelst
*/

#include <iostream>

#include <gismo.h>

namespace gismo{
namespace expr{

template<class T>
class curve_tangent_expr : public _expr<curve_tangent_expr<T> >
{
public:
    typedef T Scalar;

    typename gsGeometryMap<Scalar>::Nested_t _S;
    typename gsGeometryMap<Scalar>::Nested_t _C;

public:
    enum{ Space = 0, ScalarValued= 0, ColBlocks= 0};

    curve_tangent_expr(const gsGeometryMap<Scalar> & S,
                       const gsGeometryMap<Scalar> & C)
    :
    _S(S),
    _C(C),
    Spatch(_S.source().piece(0)) //////// <<<<<<<<---------- REMOVE THIS! MAKE PATCH-INDEPENDENT PRE-COMPUTATION
    {
        GISMO_ASSERT(_S.source().domainDim() == 2, "The geometry must be a surface with domain dimension 2");
        GISMO_ASSERT(_C.source().domainDim() == 1, "The geometry must be a curve with domain dimension 1");
        GISMO_ASSERT(_C.source().targetDim() == _S.source().domainDim(), "Curve target dimension must match surface domain dimension");
    }

    mutable gsMatrix<Scalar> theta, dtheta, sJac, res;
    const gsFunctionSet<Scalar> & Spatch;
#   define Eigen gsEigen
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
#   undef Eigen


    const gsMatrix<Scalar> & eval(const index_t k) const
    {
        // NEEDED FOR ARTIFICIAL ADDITION OF EXTRA ZERO FOR 2D
        res.resize(3,1);

        // The parametric coordinates of the surface are the evaluation of the curve at point k
        theta = _C.data().values[0].col(k);
        dtheta = _C.data().values[1].reshapeCol(k,1,_C.data().dim.second).transpose();
        sJac  = Spatch.deriv(theta);
        sJac.resize(_S.source().domainDim(),_S.source().targetDim());
        sJac.transposeInPlace();

        if (_S.source().targetDim() == 2) //if the surface is a 2-D physical entity
        {
            res.col(0).head(2) = sJac*dtheta;
            res(2,0)           = 0;
        }
        else if (_S.source().targetDim()==3) //if the surface is a 3-D physical entity
        {
            res = sJac * dtheta;
        }
        else
        {
            GISMO_ERROR("The target dimension of the surface must be 2 or 3, but is "<<_S.source().targetDim());
        }

        return res;
    }

    index_t rows() const { return 3; }

    index_t cols() const { return 1; }

    void parse(gsExprHelper<Scalar> & evList) const
    {
        evList.add(_C);
        _C.data().flags |= NEED_VALUE | NEED_DERIV;
    }

    const gsFeSpace<Scalar> & rowVar() const {return gsNullExpr<Scalar>::get();}
    const gsFeSpace<Scalar> & colVar() const {return gsNullExpr<Scalar>::get();}

    void print(std::ostream &os) const { os << "ctv("; _S.print(os); os<<" , "; _C.print(os); os <<")"; }
};

template<class T>
class curve_normal_expr : public _expr<curve_normal_expr<T> >
{
public:
    typedef T Scalar;

    typename gsGeometryMap<Scalar>::Nested_t _S;
    typename gsGeometryMap<Scalar>::Nested_t _C;

public:
    enum{ Space = 0, ScalarValued= 0, ColBlocks= 0};

    curve_normal_expr(const gsGeometryMap<Scalar> & S,
                      const gsGeometryMap<Scalar> & C)
    :
    _S(S),
    _C(C),
    Spatch(S.source().piece(0)) //////// <<<<<<<<---------- REMOVE THIS! MAKE PATCH-INDEPENDENT PRE-COMPUTATION
    {
        GISMO_ASSERT(_S.source().domainDim() == 2, "The geometry must be a surface with domain dimension 2, but is "<<_S.source().domainDim());
        GISMO_ASSERT(_C.source().domainDim() == 1, "The geometry must be a curve with domain dimension 1, but is "<<_C.source().domainDim());
        GISMO_ASSERT(_C.source().targetDim() == _S.source().domainDim(), "Curve target dimension must match surface domain dimension");
    }

    mutable gsMatrix<Scalar> theta, sJac, res;
    const gsFunctionSet<Scalar> & Spatch;
#   define Eigen gsEigen
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
#   undef Eigen


    // Make this work for a 2D planar surface AND for a 3D surface:
    // - 2D: normal vector is always (0,0,1)
    // - 3D: normal vector is the cross product of the tangent vectors
    // Make the selection between 2D and 3D based on enable_if (and get the expression templated over d)
    // to avoid if-statement
    const gsMatrix<Scalar> & eval(const index_t k) const
    {
        if (_S.source().targetDim() == 2) //if the surface is a 2-D physical entity
        {
            res.resize(3,1);
            res << 0, 0, 1;
            return res;
        }
        else if (_S.source().targetDim()==3) //if the surface is a 3-D physical entity
        {
            // The parametric coordinates of the surface are the evaluation of the curve at point k
            theta = _C.data().values[0].col(k);
            sJac  = Spatch.deriv(theta);
            sJac.resize(_S.source().domainDim(),_S.source().targetDim());
            sJac.transposeInPlace();
            res = sJac.col3d(0).cross(sJac.col3d(1));
            //res.normalize();
            return res;
        }
        else
            GISMO_ERROR("The target dimension of the surface must be 2 or 3, but is "<<_S.source().targetDim());
    }

    index_t rows() const { return 3; }

    index_t cols() const { return 1; }

    void parse(gsExprHelper<Scalar> & evList) const
    {
        evList.add(_C);
        _C.data().flags |= NEED_VALUE;
    }

    const gsFeSpace<Scalar> & rowVar() const {return gsNullExpr<Scalar>::get();}
    const gsFeSpace<Scalar> & colVar() const {return gsNullExpr<Scalar>::get();}

    void print(std::ostream &os) const { os << "cnv("; _S.print(os); os<<" , "; _C.print(os); os <<")"; }
};

template<class T>
class curve_binormal_expr : public _expr<curve_binormal_expr<T> >
{
public:
    typedef T Scalar;

    typename gsGeometryMap<Scalar>::Nested_t _S;
    typename gsGeometryMap<Scalar>::Nested_t _C;

public:
    enum{ Space = 0, ScalarValued= 0, ColBlocks= 0};

    curve_binormal_expr(const gsGeometryMap<Scalar> & S,
                      const gsGeometryMap<Scalar> & C): _S(S),_C(C) {}

    mutable gsMatrix<Scalar> cnvec, ctvec, res;
    //const gsFunctionSet<Scalar> & Spatch;
#   define Eigen gsEigen
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
#   undef Eigen


    // Make the selection between 2D and 3D based on enable_if (and get the expression templated over d)
    // to avoid if-statement
    const gsMatrix<Scalar> & eval(const index_t k) const
    {
        curve_normal_expr<Scalar>  cnv(_S,_C);
        cnvec = cnv.eval(k);
        curve_tangent_expr<Scalar> ctv(_S,_C);
        ctvec = ctv.eval(k);
        res = ctvec.col3d(0).cross(cnvec.col3d(0));
        return res;
    }

    index_t rows() const { return 3; }

    index_t cols() const { return 1; }

    void parse(gsExprHelper<Scalar> & evList) const
    {
        // WHY IS THIS NEEDED?
        evList.add(_C);
        _C.data().flags |= NEED_VALUE;

        // For pre-computation of expressions for curve_normal_expr and curve_tangent_expr
        cnv(_S,_C).parse(evList);
        ctv(_S,_C).parse(evList);
    }

    const gsFeSpace<Scalar> & rowVar() const {return gsNullExpr<Scalar>::get();}
    const gsFeSpace<Scalar> & colVar() const {return gsNullExpr<Scalar>::get();}

    void print(std::ostream &os) const { os << "cnv("; _S.print(os); os<<" , "; _C.print(os); os <<")"; }
};

/**
 * @brief First variation of basis vectors wrt DOFs.
 */

template<class E>
class curve_tangent_var1_expr : public _expr<curve_tangent_var1_expr<E> >
{
public:
    typedef typename E::Scalar Scalar;

    typename E::Nested_t _u;
    typename gsGeometryMap<Scalar>::Nested_t _S;
    typename gsGeometryMap<Scalar>::Nested_t _C;

public:
    enum{ Space = E::Space, ScalarValued= 0, ColBlocks= 0};

    curve_tangent_var1_expr(const E & u,
                            const gsGeometryMap<Scalar> & S,
                            const gsGeometryMap<Scalar> & C): _u(u), _S(S), _C(C)
    {
        GISMO_ASSERT(_S.source().domainDim() == 2, "The geometry must be a surface with domain dimension 2");
        GISMO_ASSERT(_C.source().domainDim() == 1, "The geometry must be a curve with domain dimension 1");
        GISMO_ASSERT(_C.source().targetDim() == _S.source().domainDim(), "Curve target dimension must match surface domain dimension");
    }

    mutable gsMatrix<Scalar> theta, dtheta, bGrads, res;

#   define Eigen gsEigen
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
#   undef Eigen

    // helper function
    static inline gsVector<Scalar,3> vecFun(index_t pos, Scalar val)
    {
        gsVector<Scalar,3> result = gsVector<Scalar,3>::Zero();
        result[pos] = val;
        return result;
    }

    const gsMatrix<Scalar> & eval(const index_t k) const
    {
        theta = _C.data().values[0].col(k);
        dtheta = _C.data().values[1].col(k);

        const gsMultiBasis<Scalar> & mb = static_cast<const gsMultiBasis<Scalar>&>(_u.source());
        GISMO_ASSERT( dynamic_cast<const gsMultiBasis<Scalar>*>(&_u.source()), "error");

        const index_t A = mb.basis(0).active(theta).rows(); // _u.data().actives.rows()

        res.resize(A*_u.dim(),cols()); // cols==3
        res.setZero();

        bGrads = mb.basis(0).deriv(theta); // bGrads = _u.data().values[1].col(k);

        for (index_t d = 0; d!= _u.dim(); ++d) // for all basis function components
        {
            const short_t s = d*A;
            for (index_t j = 0; j!= A; ++j) // for all actives
            {
                res.row(s+j) = vecFun(d, bGrads.at(2*j  ))* dtheta(0,0)
                             + vecFun(d, bGrads.at(2*j+1))* dtheta(1,0);
            }
        }
        return res;
    }

    index_t rows() const { return 1; }

    index_t cols() const { return 3; }

    index_t cardinality_impl() const { return _u.cardinality_impl(); }

    void parse(gsExprHelper<Scalar> & evList) const
    {
        evList.add(_C);
        _C.data().flags |= NEED_VALUE | NEED_DERIV;
    }

    const gsFeSpace<Scalar> & rowVar() const {return _u.rowVar();}
    const gsFeSpace<Scalar> & colVar() const {return gsNullExpr<Scalar>::get();}

    void print(std::ostream &os) const { os << "ctv("; _S.print(os); os<<" , "; _C.print(os); os <<")"; }
};

template<class E>
class curve_normal_var1_expr : public _expr<curve_normal_var1_expr<E> >
{
public:
    typedef typename E::Scalar Scalar;

private:
    typename E::Nested_t _u;
    typename gsGeometryMap<Scalar>::Nested_t _S;
    typename gsGeometryMap<Scalar>::Nested_t _C;

public:
    enum{ Space = E::Space, ScalarValued= 0, ColBlocks= 0};

    curve_normal_var1_expr(const E & u, const gsGeometryMap<Scalar> & S, const gsGeometryMap<Scalar> & C): _u(u), _S(S), _C(C),
    Spatch(_S.source().piece(0)) {} //////// <<<<<<<<---------- REMOVE THIS! MAKE PATCH-INDEPENDENT PRE-COMPUTATION

    mutable gsMatrix<Scalar> theta, bGrads, sJac, cnvec, res, var;
    const gsFunctionSet<Scalar> & Spatch;

#   define Eigen gsEigen
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
#   undef Eigen

    // helper function
    static inline gsVector<Scalar,3> vecFun(index_t pos, Scalar val)
    {
        gsVector<Scalar,3> result = gsVector<Scalar,3>::Zero();
        result[pos] = val;
        return result;
    }

    const gsMatrix<Scalar> & eval(const index_t k) const
    {
        theta  = _C.data().values[0].col(k);

        const gsMultiBasis<Scalar> & mb = static_cast<const gsMultiBasis<Scalar>&>(_u.source());
        GISMO_ASSERT( dynamic_cast<const gsMultiBasis<Scalar>*>(&_u.source()), "error");
        const index_t A = mb.basis(0).active(theta).rows(); // _u.data().actives.rows()

        res.resize(A*_u.dim(), cols());

        if (_S.source().targetDim() == 2) // if the surface is a 2-D physical entity
        {
            res.setZero();
            return res;
        }
        else if (_S.source().targetDim()==3) //if the surface is a 3-D physical entity
        {
            bGrads = mb.basis(0).deriv(theta);
            sJac  = Spatch.deriv(theta);
            sJac.resize(_S.source().domainDim(),_S.source().targetDim());
            sJac.transposeInPlace();
            curve_normal_expr<Scalar> cnv(_S,_C);
            cnvec = cnv.eval(k);
            const Scalar measure = cnvec.norm();
            cnvec.normalize();

            for (index_t d = 0; d!= cols(); ++d) // for all basis function components
            {
                const short_t s = d*A;
                for (index_t j = 0; j!= A; ++j) // for all active basis functions
                {
                    //first variation of non-unit normal vector
                    res.row(s+j) =  (vecFun(d, bGrads.at(2*j  )).cross(sJac.col3d(1))
                                   - vecFun(d, bGrads.at(2*j+1)).cross(sJac.col3d(0))).transpose();

                    //first variation of non-unit normal vector (divided by measure)
                    // var =  (vecFun(d, bGrads.at(2*j  )).cross(sJac.col3d(1))
                    //       - vecFun(d, bGrads.at(2*j+1)).cross(sJac.col3d(0)))/measure;

                    //first variation of unit normal vector
                    // res.row(s+j) = (var - ((cnvec.col3d(0)*var.transpose()) * cnvec.col3d(0))).transpose();
                }
            }
            return res;
        }
        else
            GISMO_ERROR("The target dimension of the surface must be 2 or 3, but is "<<_S.source().targetDim());
    }

    index_t rows() const { return 1; }

    index_t cols() const { return 3; }

    void parse(gsExprHelper<Scalar> & evList) const
    {
        evList.add(_C);
        _C.data().flags |= NEED_VALUE;
        cnv(_S,_C).parse(evList);
    }

    const gsFeSpace<Scalar> & rowVar() const { return _u.rowVar(); }
    const gsFeSpace<Scalar> & colVar() const {return gsNullExpr<Scalar>::get();}
    index_t cardinality_impl() const { return _u.cardinality_impl(); }

    void print(std::ostream &os) const { os << "var1("; _u.print(os); os <<")"; }

};

template<class E>
class curve_binormal_var1_expr : public _expr<curve_binormal_var1_expr<E> >
{
public:
    typedef typename E::Scalar Scalar;

private:
    typename E::Nested_t _u;
    typename gsGeometryMap<Scalar>::Nested_t _S;
    typename gsGeometryMap<Scalar>::Nested_t _C;

public:
    enum{ Space = E::Space, ScalarValued= 0, ColBlocks= 0};

    curve_binormal_var1_expr(const E & u, const gsGeometryMap<Scalar> & S, const gsGeometryMap<Scalar> & C): _u(u), _S(S), _C(C)
    {}

    mutable gsMatrix<Scalar> ctvec, cnvec, cbvec, ctvec_var1, cnvec_var1, theta, var, res;

#   define Eigen gsEigen
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
#   undef Eigen

    // helper function
    static inline gsVector<Scalar,3> vecFun(index_t pos, Scalar val)
    {
        gsVector<Scalar,3> result = gsVector<Scalar,3>::Zero();
        result[pos] = val;
        return result;
    }

    // Make the selection between 2D and 3D based on enable_if (and get the expression templated over d) to avoid if-statement
    const gsMatrix<Scalar> & eval(const index_t k) const
    {
        theta = _C.data().values[0].col(k);

        const gsMultiBasis<Scalar> & mb = static_cast<const gsMultiBasis<Scalar>&>(_u.source());
        GISMO_ASSERT( dynamic_cast<const gsMultiBasis<Scalar>*>(&_u.source()), "error");
        const index_t A = mb.basis(0).active(theta).rows(); // _u.data().actives.rows()

        var.resize(A*_u.dim(), cols());
        res.resize(A*_u.dim(), cols());

        curve_tangent_expr<Scalar> ctv(_S,_C);
        ctvec = ctv.eval(k);
        curve_normal_expr<Scalar> cnv(_S,_C);
        cnvec = cnv.eval(k);
        curve_tangent_var1_expr<E> ctv_var1(_u,_S,_C);
        ctvec_var1 = ctv_var1.eval(k).transpose();
        curve_normal_var1_expr<E> cnv_var1(_u,_S,_C);
        cnvec_var1 = cnv_var1.eval(k).transpose();
        curve_binormal_expr<Scalar> cbv(_S,_C);
        cbvec = cbv.eval(k);
        const Scalar measure = cbvec.norm();
        cbvec.normalize();

        // Loop over rows of (ctv1Mat or cnv1Mat) == cardinality of space
        for (index_t i = 0; i!= A*_u.dim(); ++i)
        {
            //first variation of non-unit binormal vector
            res.row(i) = ctvec_var1.col3d(i).cross(cnvec.col3d(0)) +
                         ctvec.col3d(0).cross(cnvec_var1.col3d(i));

            //first variation of non-unit binormal vector (divided by measure)
            // var.row(i) = (ctvec_var1.col3d(i).cross(cnvec.col3d(0)) +
                          // ctvec.col3d(0).cross(cnvec_var1.col3d(i)))/measure;

            //first variation of unit binormal vector
            // res.row(i) = (var.row(i).transpose() - (cbvec.col3d(0)*var.row(i)) *cbvec.col3d(0)).transpose();
        }
        return res;
    }

    index_t rows() const { return 1; }

    index_t cols() const { return 3; }

    void parse(gsExprHelper<Scalar> & evList) const
    {
        evList.add(_C);
        _C.data().flags |= NEED_VALUE;
        ctv(_S,_C).parse(evList);
        cnv(_S,_C).parse(evList);
        ctv_var1(_u,_S,_C).parse(evList);
        cnv_var1(_u,_S,_C).parse(evList);
        cbv(_S,_C).parse(evList);
    }

    const gsFeSpace<Scalar> & rowVar() const { return _u.rowVar(); }
    const gsFeSpace<Scalar> & colVar() const {return gsNullExpr<Scalar>::get();}
    index_t cardinality_impl() const { return _u.cardinality_impl(); }

    void print(std::ostream &os) const { os << "var1("; _u.print(os); os <<")"; }

};

/**
 * @brief Derivative of basis vectors along the curve.
 */

template<class T>
class curve_tangent_varalong_expr : public _expr<curve_tangent_varalong_expr<T> >
{
public:
    typedef T Scalar;

    typename gsGeometryMap<Scalar>::Nested_t _S;
    typename gsGeometryMap<Scalar>::Nested_t _C;

public:
    enum{ Space = 0, ScalarValued= 0, ColBlocks= 0};

    curve_tangent_varalong_expr(const gsGeometryMap<Scalar> & S, const gsGeometryMap<Scalar> & C): _S(S), _C(C),
    Spatch(S.source().piece(0)) //////// <<<<<<<<---------- REMOVE THIS! MAKE PATCH-INDEPENDENT PRE-COMPUTATION
    {}

    mutable gsMatrix<Scalar> theta, ddtheta, sJac, sHess, res;
    mutable gsVector<Scalar,3> vec;
    mutable gsEigen::Array<Scalar,2,1> dtheta;
    const gsFunctionSet<Scalar> & Spatch;
#   define Eigen gsEigen
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
#   undef Eigen

    const gsMatrix<Scalar> & eval(const index_t k) const
    {
        index_t numRows = _S.source().domainDim()*(_S.source().domainDim()+1)/2;
        theta = _C.data().values[0].col(k);
        dtheta = _C.data().values[1].col(k);
        ddtheta = _C.data().values[2].col(k);
        sJac  = Spatch.deriv(theta);
        sJac.resize(_S.source().domainDim(),_S.source().targetDim());
        sJac.transposeInPlace();
        sHess = Spatch.deriv2(theta);
        sHess.resize(numRows,_S.source().targetDim());
        sHess.transposeInPlace();
        vec.head(2) = dtheta.pow(2);
        vec(2) = 2*dtheta(0,0)*dtheta(1,0);
        
        res = sJac*ddtheta + sHess*vec;
        //res.normalize();
        return res;
    }

    index_t rows() const { return 3; }

    index_t cols() const { return 1; }

    void parse(gsExprHelper<Scalar> & evList) const
    {
        evList.add(_C);
        _C.data().flags |= NEED_VALUE | NEED_DERIV | NEED_2ND_DER;
    }

    const gsFeSpace<Scalar> & rowVar() const {return gsNullExpr<Scalar>::get();}
    const gsFeSpace<Scalar> & colVar() const {return gsNullExpr<Scalar>::get();}

    void print(std::ostream &os) const { os << "ctv_varalong("; _S.print(os); os<<" , "; _C.print(os); os <<")"; }
};

template<class T>
class curve_normal_varalong_expr : public _expr<curve_normal_varalong_expr<T> >
{
public:
    typedef T Scalar;

    typename gsGeometryMap<Scalar>::Nested_t _S;
    typename gsGeometryMap<Scalar>::Nested_t _C;

public:
    enum{ Space = 0, ScalarValued= 0, ColBlocks= 0};

    curve_normal_varalong_expr(const gsGeometryMap<Scalar> & S,
                                const gsGeometryMap<Scalar> & C): _S(S), _C(C),
    Spatch(S.source().piece(0)) //////// <<<<<<<<---------- REMOVE THIS! MAKE PATCH-INDEPENDENT PRE-COMPUTATION
    {}

    mutable gsMatrix<Scalar> theta, dtheta, sJac, sHess, res, vecA, vecB;
    const gsFunctionSet<Scalar> & Spatch;
#   define Eigen gsEigen
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
#   undef Eigen

    const gsMatrix<Scalar> & eval(const index_t k) const
    {
        index_t numRows = _S.source().domainDim()*(_S.source().domainDim()+1)/2;
        theta = _C.data().values[0].col(k);
        dtheta = _C.data().values[1].col(k);
        sJac  = Spatch.deriv(theta);
        sJac.resize(_S.source().domainDim(),_S.source().targetDim());
        sJac.transposeInPlace();
        sHess = Spatch.deriv2(theta);
        sHess.resize(numRows,_S.source().targetDim());
        sHess.transposeInPlace();
        sHess.col(2).swap(sHess.col(1));
        vecA = sHess.leftCols(2)*dtheta;
        vecB = sHess.rightCols(2)*dtheta;

        res  = vecA.col3d(0).cross(sJac.col3d(1)) + sJac.col3d(0).cross(vecB.col3d(0));
        //res.normalize();
        return res;
    }

    index_t rows() const { return 3; }

    index_t cols() const { return 1; }

    void parse(gsExprHelper<Scalar> & evList) const
    {
        evList.add(_C);
        _C.data().flags |= NEED_VALUE | NEED_DERIV;
    }

    const gsFeSpace<Scalar> & rowVar() const {return gsNullExpr<Scalar>::get();}
    const gsFeSpace<Scalar> & colVar() const {return gsNullExpr<Scalar>::get();}

    void print(std::ostream &os) const { os << "cnv_varalong("; _S.print(os); os<<" , "; _C.print(os); os <<")"; }
};

template<class T>
class curve_binormal_varalong_expr : public _expr<curve_binormal_varalong_expr<T> >
{
public:
    typedef T Scalar;

    typename gsGeometryMap<Scalar>::Nested_t _S;
    typename gsGeometryMap<Scalar>::Nested_t _C;

public:
    enum{ Space = 0, ScalarValued= 0, ColBlocks= 0};

    curve_binormal_varalong_expr(const gsGeometryMap<Scalar> & S, const gsGeometryMap<Scalar> & C): _S(S), _C(C) {}

    mutable gsMatrix<Scalar> ctvec, ctvec_vara, cnvec, cnvec_vara, res;
#   define Eigen gsEigen
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
#   undef Eigen

    const gsMatrix<Scalar> & eval(const index_t k) const
    {
        curve_tangent_varalong_expr<T> ctv_vara(_S,_C);
        ctvec_vara = ctv_vara.eval(k);
        curve_normal_varalong_expr<T> cnv_vara(_S,_C);
        cnvec_vara = cnv_vara.eval(k);
        curve_tangent_expr<T> ctv(_S,_C);
        ctvec = ctv.eval(k);
        curve_normal_expr<T> cnv(_S,_C);
        cnvec = cnv.eval(k);

        res = ctvec_vara.col3d(0).cross(cnvec.col3d(0)) +
              ctvec.col3d(0).cross(cnvec_vara.col3d(0));
        //res.normalize();
        return res;
    }

    index_t rows() const { return 3; }

    index_t cols() const { return 1; }

    void parse(gsExprHelper<Scalar> & evList) const
    {
        // WHY IS THIS NEEDED?
        evList.add(_C);
        _C.data().flags |= NEED_VALUE;

        // For pre-computation of expressions for curve_normal_expr and curve_tangent_expr
        ctv(_S,_C).parse(evList);
        ctv_vara(_S,_C).parse(evList);
        cnv(_S,_C).parse(evList);
        cnv_vara(_S,_C).parse(evList);
    }

    const gsFeSpace<Scalar> & rowVar() const {return gsNullExpr<Scalar>::get();}
    const gsFeSpace<Scalar> & colVar() const {return gsNullExpr<Scalar>::get();}

    void print(std::ostream &os) const { os << "cbv_varalong("; _S.print(os); os<<" , "; _C.print(os); os <<")"; }
};

/**
 * @brief First variation of the derivative of basis vectors along the curve wrt DOFs.
 */

template<class E>
class curve_tangent_varalong_var1_expr : public _expr<curve_tangent_varalong_var1_expr<E> >
{
public:
    typedef typename E::Scalar Scalar;

private:
    typename E::Nested_t _u;
    typename gsGeometryMap<Scalar>::Nested_t _S;
    typename gsGeometryMap<Scalar>::Nested_t _C;

public:
    enum{ Space = E::Space, ScalarValued= 0, ColBlocks= 0};

    curve_tangent_varalong_var1_expr(const E & u, const gsGeometryMap<Scalar> & S, const gsGeometryMap<Scalar> & C): _u(u), _S(S), _C(C) {}

    mutable gsMatrix<Scalar> theta, res;
    mutable gsEigen::Array<Scalar,2,1> dtheta;
    mutable gsVector<Scalar> bGrads, bHess, ddtheta;
    mutable gsVector<Scalar,3> vec;
    mutable Scalar dotA, dotB;
#   define Eigen gsEigen
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
#   undef Eigen

    // helper function
    static inline gsVector<Scalar,3> vecFun(index_t pos, Scalar val)
    {
        gsVector<Scalar,3> result = gsVector<Scalar,3>::Zero();
        result[pos] = val;
        return result;
    }

    const gsMatrix<Scalar> & eval(const index_t k) const
    {
        theta  = _C.data().values[0].col(k);
        dtheta = _C.data().values[1].col(k);
        ddtheta = _C.data().values[2].col(k);

        const gsMultiBasis<Scalar> & mb = static_cast<const gsMultiBasis<Scalar>&>(_u.source());
        GISMO_ASSERT( dynamic_cast<const gsMultiBasis<Scalar>*>(&_u.source()), "error");
        const index_t A = mb.basis(0).active(theta).rows(); 
        bGrads = mb.basis(0).deriv(theta);
        bHess =  mb.basis(0).deriv2(theta);
        vec.head(2) = dtheta.pow(2);
        vec(2) = 2*dtheta(0,0)*dtheta(1,0);

        res.resize(A*_u.dim(), cols());

        for (index_t d = 0; d!= cols(); ++d) // for all basis function components
        {
            const short_t s = d*A;
            for (index_t j = 0; j!= A; ++j) // for all active basis functions
            {
                dotA = bHess.middleRows(3*j,3).dot(vec);
                dotB = bGrads.middleRows(2*j,2).dot(ddtheta);
                res.row(s+j) = vecFun(d, dotA) + vecFun(d, dotB);
            }
        }
        return res;
    }

    index_t rows() const { return 1; }
    index_t cols() const { return 3; }

    void parse(gsExprHelper<Scalar> & evList) const
    {
        evList.add(_C);
        _C.data().flags |= NEED_VALUE | NEED_DERIV | NEED_2ND_DER;
    }

    const gsFeSpace<Scalar> & rowVar() const { return _u.rowVar(); }
    const gsFeSpace<Scalar> & colVar() const {return gsNullExpr<Scalar>::get();}
    index_t cardinality_impl() const { return _u.cardinality_impl(); }

    void print(std::ostream &os) const { os << "curve_tangent_varalong_var1("; _u.print(os); os <<")"; }

};

template<class E>
class curve_normal_varalong_var1_expr : public _expr<curve_normal_varalong_var1_expr<E> >
{
public:
    typedef typename E::Scalar Scalar;

private:
    typename E::Nested_t _u;
    typename gsGeometryMap<Scalar>::Nested_t _S;
    typename gsGeometryMap<Scalar>::Nested_t _defS;
    typename gsGeometryMap<Scalar>::Nested_t _C;

public:
    enum{ Space = E::Space, ScalarValued= 0, ColBlocks= 0};

    curve_normal_varalong_var1_expr(const E & u, const gsGeometryMap<Scalar> & S,
                                     const gsGeometryMap<Scalar> & defS,
                                     const gsGeometryMap<Scalar> & C): _u(u), _S(S), _defS(defS), _C(C),
                                     Spatch(S.source().piece(0)),defSpatch(defS.source().piece(0)) //////// <<<<<<<<---------- REMOVE THIS! MAKE PATCH-INDEPENDENT PRE-COMPUTATION
                                     {}

    mutable gsMatrix<Scalar> theta, dtheta, bGrads, bHess, localbHess, uJac, uHess, res;
    mutable gsMatrix<Scalar> SJac, defSJac, SHess, defSHess;
    mutable gsVector<Scalar,2> dot;
    mutable gsVector<Scalar,3> dotA, dotB;
    const gsFunctionSet<Scalar> & Spatch;
    const gsFunctionSet<Scalar> & defSpatch;

#   define Eigen gsEigen
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
#   undef Eigen

    // helper function
    static inline gsVector<Scalar,3> vecFun(index_t pos, Scalar val)
    {
        gsVector<Scalar,3> result = gsVector<Scalar,3>::Zero();
        result[pos] = val;
        return result;
    }

    const gsMatrix<Scalar> & eval(const index_t k) const
    {
        theta  = _C.data().values[0].col(k);
        dtheta = _C.data().values[1].col(k);

        const gsMultiBasis<Scalar> & mb = static_cast<const gsMultiBasis<Scalar>&>(_u.source());
        GISMO_ASSERT( dynamic_cast<const gsMultiBasis<Scalar>*>(&_u.source()), "error");
        const index_t A = mb.basis(0).active(theta).rows(); // _u.data().actives.rows()
        bGrads = mb.basis(0).deriv(theta);
        bHess =  mb.basis(0).deriv2(theta);

        index_t numRows = _S.source().domainDim()*(_S.source().domainDim()+1)/2;
        SJac  = Spatch.deriv(theta);
        SJac.resize(_S.source().domainDim(),_S.source().targetDim());
        SJac.transposeInPlace();
        defSJac  = defSpatch.deriv(theta);
        defSJac.resize(_defS.source().domainDim(),_defS.source().targetDim());
        defSJac.transposeInPlace();
        uJac = defSJac - SJac;
        SHess = Spatch.deriv2(theta);
        SHess.resize(numRows,_S.source().targetDim());
        SHess.transposeInPlace();
        defSHess = defSpatch.deriv2(theta);
        defSHess.resize(numRows,_defS.source().targetDim());
        defSHess.transposeInPlace();
        uHess = defSHess - SHess;

        uHess.col(2).swap(uHess.col(1));
        dotA = uHess.leftCols(2) * dtheta;
        dotB = uHess.rightCols(2) * dtheta;

        res.resize(A*_u.dim(), cols());

        for (index_t d = 0; d!= cols(); ++d) // for all basis function components
        {
            const short_t s = d*A;
            for (index_t j = 0; j!= A; ++j) // for all active basis functions
            {
                secDerToHessian(bHess.col(0).segment(3*j,3),2,localbHess);
                localbHess.resize(2,2);
                dot = localbHess * dtheta;

                //not normalized result
                res.row(s+j) =  vecFun(d,dot(0)).cross(uJac.col3d(1)) +
                                dotA.cross(vecFun(d,bGrads.at(2*j+1))) +
                                vecFun(d,bGrads.at(2*j  )).cross(dotB) +
                                uJac.col3d(0).cross(vecFun(d,dot(1)));
            }
        }
        return res;
    }

    index_t rows() const { return 1; }

    index_t cols() const { return 3; }

    void parse(gsExprHelper<Scalar> & evList) const
    {
        evList.add(_C);
        _C.data().flags |= NEED_VALUE | NEED_DERIV;
        cnv_vara(_S,_C).parse(evList);
    }

    const gsFeSpace<Scalar> & rowVar() const { return _u.rowVar(); }
    const gsFeSpace<Scalar> & colVar() const {return gsNullExpr<Scalar>::get();}
    index_t cardinality_impl() const { return _u.cardinality_impl(); }

    void print(std::ostream &os) const { os << "curve_normal_varalong_var1("; _u.print(os); os <<")"; }

};

template<class E>
class curve_normal_varalong_var1_normalized_expr : public _expr<curve_normal_varalong_var1_normalized_expr<E> >
{
public:
    typedef typename E::Scalar Scalar;

private:
    typename E::Nested_t _u;
    typename gsGeometryMap<Scalar>::Nested_t _S;
    typename gsGeometryMap<Scalar>::Nested_t _defS;
    typename gsGeometryMap<Scalar>::Nested_t _C;

public:
    enum{ Space = E::Space, ScalarValued= 0, ColBlocks= 0};

    curve_normal_varalong_var1_normalized_expr(const E & u, const gsGeometryMap<Scalar> & S,
                                                const gsGeometryMap<Scalar> & defS,
                                                const gsGeometryMap<Scalar> & C): _u(u), _S(S), _defS(defS), _C(C),
                                                Spatch(S.source().piece(0)),defSpatch(defS.source().piece(0)) //////// <<<<<<<<---------- REMOVE THIS! MAKE PATCH-INDEPENDENT PRE-COMPUTATION
                                                {}

    mutable gsMatrix<Scalar> theta, dtheta, bGrads, bHess, localbHess, uJac, uHess, res;
    mutable gsMatrix<Scalar> SJac, defSJac, SHess, defSHess;
    mutable gsMatrix<Scalar> cnvec_vara, var;
    mutable gsVector<Scalar,2> dot;
    mutable gsVector<Scalar,3> dotA, dotB;
    const gsFunctionSet<Scalar> & Spatch;
    const gsFunctionSet<Scalar> & defSpatch;

#   define Eigen gsEigen
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
#   undef Eigen

    // helper function
    static inline gsVector<Scalar,3> vecFun(index_t pos, Scalar val)
    {
        gsVector<Scalar,3> result = gsVector<Scalar,3>::Zero();
        result[pos] = val;
        return result;
    }

    const gsMatrix<Scalar> & eval(const index_t k) const
    {
        theta  = _C.data().values[0].col(k);
        dtheta = _C.data().values[1].col(k);

        const gsMultiBasis<Scalar> & mb = static_cast<const gsMultiBasis<Scalar>&>(_u.source());
        GISMO_ASSERT( dynamic_cast<const gsMultiBasis<Scalar>*>(&_u.source()), "error");
        const index_t A = mb.basis(0).active(theta).rows(); // _u.data().actives.rows()
        bGrads = mb.basis(0).deriv(theta);
        bHess =  mb.basis(0).deriv2(theta);

        index_t numRows = _S.source().domainDim()*(_S.source().domainDim()+1)/2;
        SJac  = Spatch.deriv(theta);
        SJac.resize(_S.source().domainDim(),_S.source().targetDim());
        SJac.transposeInPlace();
        defSJac  = defSpatch.deriv(theta);
        defSJac.resize(_defS.source().domainDim(),_defS.source().targetDim());
        defSJac.transposeInPlace();
        uJac = defSJac - SJac;
        SHess = Spatch.deriv2(theta);
        SHess.resize(numRows,_S.source().targetDim());
        SHess.transposeInPlace();
        defSHess = defSpatch.deriv2(theta);
        defSHess.resize(numRows,_defS.source().targetDim());
        defSHess.transposeInPlace();
        uHess = defSHess - SHess;

        uHess.col(2).swap(uHess.col(1));
        dotA = uHess.leftCols(2) * dtheta;
        dotB = uHess.rightCols(2) * dtheta;

        curve_normal_varalong_expr<Scalar> cnv_vara(_S,_C);
        cnvec_vara = cnv_vara.eval(k);
        const Scalar measure = cnvec_vara.norm();
        cnvec_vara.normalize();

        res.resize(A*_u.dim(), cols());
        var.resize(A*_u.dim(), cols());

        for (index_t d = 0; d!= cols(); ++d) // for all basis function components
        {
            const short_t s = d*A;
            for (index_t j = 0; j!= A; ++j) // for all active basis functions
            {
                secDerToHessian(bHess.col(0).segment(3*j,3),2,localbHess);
                localbHess.resize(2,2);
                dot = localbHess * dtheta;
           
                //not normalized result (divided by measure)
                var.row(s+j) =  (vecFun(d,dot(0)).cross(uJac.col3d(1)) +
                                dotA.cross(vecFun(d,bGrads.at(2*j+1))) +
                                vecFun(d,bGrads.at(2*j  )).cross(dotB) +
                                uJac.col3d(0).cross(vecFun(d,dot(1))))/measure;
                
                //normalized result
                res.row(s+j) = (var.row(s+j).transpose() - (cnvec_vara.col3d(0)*var.row(s+j)) * cnvec_vara.col3d(0)).transpose();
            }
        }
        return res;
    }

    index_t rows() const { return 1; }

    index_t cols() const { return 3; }

    void parse(gsExprHelper<Scalar> & evList) const
    {
        evList.add(_C);
        _C.data().flags |= NEED_VALUE | NEED_DERIV;
        cnv_vara(_S,_C).parse(evList);
    }

    const gsFeSpace<Scalar> & rowVar() const { return _u.rowVar(); }
    const gsFeSpace<Scalar> & colVar() const {return gsNullExpr<Scalar>::get();}
    index_t cardinality_impl() const { return _u.cardinality_impl(); }

    void print(std::ostream &os) const { os << "curve_normal_varalong_var1_normalized("; _u.print(os); os <<")"; }

};

template<class E>
class curve_binormal_varalong_var1_expr : public _expr<curve_binormal_varalong_var1_expr<E> >
{
public:
    typedef typename E::Scalar Scalar;

private:
    typename E::Nested_t _u;
    typename gsGeometryMap<Scalar>::Nested_t _S;
    typename gsGeometryMap<Scalar>::Nested_t _defS;
    typename gsGeometryMap<Scalar>::Nested_t _C;

public:
    enum{ Space = E::Space, ScalarValued= 0, ColBlocks= 0};

    curve_binormal_varalong_var1_expr(const E & u, const gsGeometryMap<Scalar> & S,
                                       const gsGeometryMap<Scalar> & defS, const gsGeometryMap<Scalar> & C):
                                       _u(u), _S(S), _defS(defS), _C(C) {}

    mutable gsMatrix<Scalar> ctvec, ctvec_vara, ctvec_var1, ctvec_vara_var1, theta;
    mutable gsMatrix<Scalar> cnvec, cnvec_vara, cnvec_var1, cnvec_vara_var1, res; 

#   define Eigen gsEigen
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
#   undef Eigen

    const gsMatrix<Scalar> & eval(const index_t k) const
    {   
        theta  = _C.data().values[0].col(k);
        const gsMultiBasis<Scalar> & mb = static_cast<const gsMultiBasis<Scalar>&>(_u.source());
        GISMO_ASSERT( dynamic_cast<const gsMultiBasis<Scalar>*>(&_u.source()), "error");
        const index_t A = mb.basis(0).active(theta).rows(); // _u.data().actives.rows()
        
        curve_tangent_expr<Scalar> ctv(_S,_C);
        ctvec = ctv.eval(k);
        curve_tangent_varalong_expr<Scalar> ctv_vara(_S,_C);
        ctvec_vara = ctv_vara.eval(k);
        curve_tangent_var1_expr<E> ctv_var1(_u,_S,_C);
        ctvec_var1 = ctv_var1.eval(k).transpose();
        curve_tangent_varalong_var1_expr<E> ctv_vara_var1(_u,_S,_C);
        ctvec_vara_var1 = ctv_vara_var1.eval(k).transpose();

        curve_normal_expr<Scalar> cnv(_S,_C);
        cnvec = cnv.eval(k);
        curve_normal_varalong_expr<Scalar> cnv_vara(_S,_C);
        cnvec_vara = cnv_vara.eval(k);
        curve_normal_var1_expr<E> cnv_var1(_u,_S,_C);
        cnvec_var1 = cnv_var1.eval(k).transpose();
        curve_normal_varalong_var1_expr<E> cnv_vara_var1(_u,_S,_defS,_C);
        cnvec_vara_var1 = cnv_vara_var1.eval(k).transpose();

        res.resize(A*_u.dim(),cols());

        for (index_t i = 0; i!= A*_u.dim(); ++i)
        {
            //not normalized result
            res.row(i) =  ctvec_vara_var1.col3d(i).cross(cnvec.col3d(0)) +
                          ctvec_vara.col3d(0).cross(cnvec_var1.col3d(i)) +
                          ctvec_var1.col3d(i).cross(cnvec_vara.col3d(0)) +
                          ctvec.col3d(0).cross(cnvec_vara_var1.col3d(i));
        }
        return res;
    }

    index_t rows() const { return 1; }

    index_t cols() const { return 3; }

    void parse(gsExprHelper<Scalar> & evList) const
    {
        evList.add(_C);
        _C.data().flags |= NEED_VALUE;
        ctv(_S,_C).parse(evList);
        ctv_vara(_S,_C).parse(evList);
        ctv_var1(_u,_S,_C).parse(evList);
        ctv_vara_var1(_u,_S,_C).parse(evList);
        cnv(_S,_C).parse(evList);
        cnv_vara(_S,_C).parse(evList);
        cnv_var1(_u,_S,_C).parse(evList);
        cnv_vara_var1(_u,_S,_defS,_C).parse(evList);
    }

    const gsFeSpace<Scalar> & rowVar() const { return _u.rowVar(); }
    const gsFeSpace<Scalar> & colVar() const {return gsNullExpr<Scalar>::get();}
    index_t cardinality_impl() const { return _u.cardinality_impl(); }

    void print(std::ostream &os) const { os << "curve_binormal_varalong_var1("; _u.print(os); os <<")"; }

};

template<class E>
class curve_binormal_varalong_var1_normalized_expr : public _expr<curve_binormal_varalong_var1_normalized_expr<E> >
{
public:
    typedef typename E::Scalar Scalar;

private:
    typename E::Nested_t _u;
    typename gsGeometryMap<Scalar>::Nested_t _S;
    typename gsGeometryMap<Scalar>::Nested_t _defS;
    typename gsGeometryMap<Scalar>::Nested_t _C;

public:
    enum{ Space = E::Space, ScalarValued= 0, ColBlocks= 0};

    curve_binormal_varalong_var1_normalized_expr(const E & u, const gsGeometryMap<Scalar> & S,
                                                  const gsGeometryMap<Scalar> & defS, const gsGeometryMap<Scalar> & C):
                                                  _u(u), _S(S), _defS(defS), _C(C) {}

    mutable gsMatrix<Scalar> ctvec, ctvec_vara, ctvec_var1, ctvec_vara_var1, theta;
    mutable gsMatrix<Scalar> cnvec, cnvec_vara, cnvec_var1, cnvec_vara_var1, res; 
    mutable gsMatrix<Scalar> cbvec_vara, var;
#   define Eigen gsEigen
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
#   undef Eigen

    const gsMatrix<Scalar> & eval(const index_t k) const
    {   
        theta  = _C.data().values[0].col(k);
        const gsMultiBasis<Scalar> & mb = static_cast<const gsMultiBasis<Scalar>&>(_u.source());
        GISMO_ASSERT( dynamic_cast<const gsMultiBasis<Scalar>*>(&_u.source()), "error");
        const index_t A = mb.basis(0).active(theta).rows();
        
        curve_tangent_expr<Scalar> ctv(_S,_C);
        ctvec = ctv.eval(k);
        curve_tangent_varalong_expr<Scalar> ctv_vara(_S,_C);
        ctvec_vara = ctv_vara.eval(k);
        curve_tangent_var1_expr<E> ctv_var1(_u,_S,_C);
        ctvec_var1 = ctv_var1.eval(k).transpose();
        curve_tangent_varalong_var1_expr<E> ctv_vara_var1(_u,_S,_C);
        ctvec_vara_var1 = ctv_vara_var1.eval(k).transpose();

        curve_normal_expr<Scalar> cnv(_S,_C);
        cnvec = cnv.eval(k);
        curve_normal_varalong_expr<Scalar> cnv_vara(_S,_C);
        cnvec_vara = cnv_vara.eval(k);
        curve_normal_var1_expr<E> cnv_var1(_u,_S,_C);
        cnvec_var1 = cnv_var1.eval(k).transpose();
        curve_normal_varalong_var1_expr<E> cnv_vara_var1(_u,_S,_defS,_C);
        cnvec_vara_var1 = cnv_vara_var1.eval(k).transpose();
        curve_binormal_varalong_expr<Scalar> cbv_vara(_S,_C);
        cbvec_vara = cbv_vara.eval(k);
        const Scalar measure = cbvec_vara.norm();
        cbvec_vara.normalize();

        var.resize(A*_u.dim(),cols());
        res.resize(A*_u.dim(),cols());

        for (index_t i = 0; i!= A*_u.dim(); ++i)
        {

            //not normalized result (divided by measure)
            var.row(i) = (ctvec_vara_var1.col3d(i).cross(cnvec.col3d(0)) +
                          ctvec_vara.col3d(0).cross(cnvec_var1.col3d(i)) +
                          ctvec_var1.col3d(i).cross(cnvec_vara.col3d(0)) +
                          ctvec.col3d(0).cross(cnvec_vara_var1.col3d(i)))/measure;

            //normalized result
            res.row(i) = (var.row(i).transpose() - (cbvec_vara.col3d(0)*var.row(i)) * cbvec_vara.col3d(0)).transpose();
        }
        return res;
    }

    index_t rows() const { return 1; }

    index_t cols() const { return 3; }

    void parse(gsExprHelper<Scalar> & evList) const
    {
        evList.add(_C);
        _C.data().flags |= NEED_VALUE;
        ctv(_S,_C).parse(evList);
        ctv_vara(_S,_C).parse(evList);
        ctv_var1(_u,_S,_C).parse(evList);
        ctv_vara_var1(_u,_S,_C).parse(evList);
        cnv(_S,_C).parse(evList);
        cnv_vara(_S,_C).parse(evList);
        cnv_var1(_u,_S,_C).parse(evList);
        cnv_vara_var1(_u,_S,_defS,_C).parse(evList);
        cbv_vara(_S,_C).parse(evList);
    }

    const gsFeSpace<Scalar> & rowVar() const { return _u.rowVar(); }
    const gsFeSpace<Scalar> & colVar() const {return gsNullExpr<Scalar>::get();}
    index_t cardinality_impl() const { return _u.cardinality_impl(); }

    void print(std::ostream &os) const { os << "curve_binormal_varalong_var1_normalized("; _u.print(os); os <<")"; }

};

/**
 * @brief Second variation of derivative of basis vector wrt DOFs times a vector.
 */

template<class E1, class E2, class E3>
class curve_normal_var2dot_expr : public _expr<curve_normal_var2dot_expr<E1,E2,E3> >
{
public:
    typedef typename E1::Scalar Scalar;

private:
    typename E1::Nested_t _u;
    typename E2::Nested_t _v;
    typename gsGeometryMap<Scalar>::Nested_t _C;
    typename E3::Nested_t _vec;

public:
    enum{ Space = 3, ScalarValued= 0, ColBlocks= 0};

    curve_normal_var2dot_expr(const E1 & u, const E2 & v, const gsGeometryMap<Scalar> & C, _expr<E3> const& vec):
                              _u(u), _v(v), _C(C), _vec(vec) {}

    mutable gsMatrix<Scalar> theta, bGradsu, bGradsv, evec, res;

#   define Eigen gsEigen
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
#   undef Eigen

    // helper function
    static inline gsVector<Scalar,3> vecFun(index_t pos, Scalar val)
    {
        gsVector<Scalar,3> result = gsVector<Scalar,3>::Zero();
        result[pos] = val;
        return result;
    }

    const gsMatrix<Scalar> & eval(const index_t k) const
    {
        theta  = _C.data().values[0].col(k);
        evec = _vec.eval(k);

        const gsMultiBasis<Scalar> & mbu = static_cast<const gsMultiBasis<Scalar>&>(_u.source());
        const gsMultiBasis<Scalar> & mbv = static_cast<const gsMultiBasis<Scalar>&>(_v.source());
        const index_t Au = mbu.basis(0).active(theta).rows(); // _u.data().actives.rows()
        const index_t Av = mbv.basis(0).active(theta).rows(); // _v.data().actives.rows()

        bGradsu = mbu.basis(0).deriv(theta);
        bGradsv = mbv.basis(0).deriv(theta);

        res.resize(Au*_u.dim(), Av*_v.dim());
        
        for (index_t d = 0; d!= _u.dim(); ++d) // for all components
        {
            const short_t s = d*Au;
            for (index_t j = 0; j!= Au; ++j) // for all actives
            {
                for (index_t e = 0; e!= _v.dim(); ++e) // for all components
                {
                    const short_t r = e*Av;
                    for (index_t k = 0; k!= Av; ++k) // for all actives
                    {   
                        res(s+j,r+k) = (vecFun(d,bGradsu.at(2*j )).cross(vecFun(e,bGradsv.at(2*k+1))) +
                                        vecFun(e,bGradsv.at(2*k )).cross(vecFun(d,bGradsu.at(2*j+1)))).dot(evec.col(0));
                    }
                }
            }
        }
        return res;
    }

    index_t rows() const { return 1; }

    index_t cols() const { return 1; }

    void parse(gsExprHelper<Scalar> & evList) const
    {
        evList.add(_C);
        _C.data().flags |= NEED_VALUE;
        _vec.parse(evList);
    }

    const gsFeSpace<Scalar> & rowVar() const { return _u.rowVar(); }
    const gsFeSpace<Scalar> & colVar() const {return _v.rowVar();}
    index_t cardinality_impl() const { return _u.cardinality_impl(); }

    void print(std::ostream &os) const { os << "curve_normal_var2dot("; _u.print(os); os <<")"; }

};

template<class E1, class E2, class E3>
class curve_binormal_var2dot_expr : public _expr<curve_binormal_var2dot_expr<E1,E2,E3> >
{
public:
    typedef typename E1::Scalar Scalar;

private:
    typename E1::Nested_t _u;
    typename E2::Nested_t _v;
    typename gsGeometryMap<Scalar>::Nested_t _S;
    typename gsGeometryMap<Scalar>::Nested_t _C;
    typename E3::Nested_t _vec;

public:
    enum{ Space = 3, ScalarValued= 0, ColBlocks= 0};

    curve_binormal_var2dot_expr(const E1 & u, const E2 & v, const gsGeometryMap<Scalar> & S, const gsGeometryMap<Scalar> & C,
                                _expr<E3> const& vec): _u(u), _v(v), _S(S), _C(C), _vec(vec) {}

    mutable gsMatrix<Scalar> theta, bGradsu, bGradsv, evec, res;
    mutable gsMatrix<Scalar> ctvec, ctvec_var1u, ctvec_var1v, cnvec_var1u, cnvec_var1v;
    mutable gsVector<Scalar,3> cnvec_var2;

#   define Eigen gsEigen
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
#   undef Eigen

    // helper function
    static inline gsVector<Scalar,3> vecFun(index_t pos, Scalar val)
    {
        gsVector<Scalar,3> result = gsVector<Scalar,3>::Zero();
        result[pos] = val;
        return result;
    }

    const gsMatrix<Scalar> & eval(const index_t k) const
    {
        theta  = _C.data().values[0].col(k);
        evec = _vec.eval(k);

        const gsMultiBasis<Scalar> & mbu = static_cast<const gsMultiBasis<Scalar>&>(_u.source());
        const gsMultiBasis<Scalar> & mbv = static_cast<const gsMultiBasis<Scalar>&>(_v.source());
        const index_t Au = mbu.basis(0).active(theta).rows();
        const index_t Av = mbv.basis(0).active(theta).rows();

        bGradsu = mbu.basis(0).deriv(theta);
        bGradsv = mbv.basis(0).deriv(theta);

        curve_tangent_expr<Scalar> ctv(_S,_C);
        ctvec = ctv.eval(k);
        curve_tangent_var1_expr<E1> ctv_var1u(_u,_S,_C);
        ctvec_var1u = ctv_var1u.eval(k).transpose();
        curve_tangent_var1_expr<E2> ctv_var1v(_v,_S,_C);
        ctvec_var1v = ctv_var1v.eval(k).transpose();
        curve_normal_var1_expr<E1> cnv_var1u(_u,_S,_C);
        cnvec_var1u = cnv_var1u.eval(k).transpose();
        curve_normal_var1_expr<E2> cnv_var1v(_v,_S,_C);
        cnvec_var1v = cnv_var1v.eval(k).transpose();
               
        res.resize(Au*_u.dim(), Av*_v.dim());
        
        for (index_t d = 0; d!= _u.dim(); ++d) // for all components
        {
            const short_t s = d*Au;
            for (index_t j = 0; j!= Au; ++j) // for all actives
            {
                for (index_t e = 0; e!= _v.dim(); ++e) // for all components
                {
                    const short_t r = e*Av;
                    for (index_t k = 0; k!= Av; ++k) // for all actives
                    {   
                        cnvec_var2 = vecFun(d,bGradsu.at(2*j )).cross(vecFun(e,bGradsv.at(2*k+1))) +
                                     vecFun(e,bGradsv.at(2*k )).cross(vecFun(d,bGradsu.at(2*j+1)));

                        res(s+j,r+k) = (ctvec_var1u.col3d(s+j).cross(cnvec_var1v.col3d(r+k)) +
                                       ctvec_var1v.col3d(r+k).cross(cnvec_var1u.col3d(s+j)) + 
                                       ctvec.col3d(0).cross(cnvec_var2)).dot(evec.col3d(0));
                    }
                }
            }
        }
        return res;
    }

    index_t rows() const { return 1; }
    index_t cols() const { return 1; }

    void parse(gsExprHelper<Scalar> & evList) const
    {
        evList.add(_C);
        _C.data().flags |= NEED_VALUE;
        ctv(_S,_C).parse(evList);
        ctv_var1(_u,_S,_C).parse(evList);
        ctv_var1(_v,_S,_C).parse(evList);
        cnv_var1(_u,_S,_C).parse(evList);
        cnv_var1(_v,_S,_C).parse(evList);
        _vec.parse(evList);
    }

    const gsFeSpace<Scalar> & rowVar() const { return _u.rowVar(); }
    const gsFeSpace<Scalar> & colVar() const {return _v.rowVar();}
    index_t cardinality_impl() const { return _u.cardinality_impl(); }

    void print(std::ostream &os) const { os << "curve_binormal_var2dot("; _u.print(os); os <<")"; }
};

/**
 * @brief Second variation of derivative of basis vector along the curve wrt DOFs times a vector.
 */

template<class E1, class E2, class E3>
class curve_normal_varalong_var2dot_expr : public _expr<curve_normal_varalong_var2dot_expr<E1,E2,E3> >
{
public:
    typedef typename E1::Scalar Scalar;

private:
    typename E1::Nested_t _u;
    typename E2::Nested_t _v;
    typename gsGeometryMap<Scalar>::Nested_t _S;
    typename gsGeometryMap<Scalar>::Nested_t _defS;
    typename gsGeometryMap<Scalar>::Nested_t _C;
    typename E3::Nested_t _vec;

public:
    enum{ Space = 3, ScalarValued= 0, ColBlocks= 0};

    curve_normal_varalong_var2dot_expr(const E1 &u, const E2 &v, const gsGeometryMap<Scalar> &S, const gsGeometryMap<Scalar> &defS,
                                       const gsGeometryMap<Scalar> &C, _expr<E3> const &vec): _u(u), _v(v), _S(S), _defS(defS), _C(C), _vec(vec)
                                       {}

    mutable gsMatrix<Scalar> theta, dtheta, bGradsu, bGradsv, bHessu, bHessv, localbHess, res;
    mutable gsMatrix<Scalar> dotu, dotv, cnvec_vara, cnvec_vara_var1u, cnvec_vara_var1v, evec;
    mutable gsMatrix<gsVector3d<Scalar>> matrix;
    mutable Scalar measure1u, measure1v, measure2;

#   define Eigen gsEigen
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
#   undef Eigen

    // helper function
    static inline gsVector<Scalar,3> vecFun(index_t pos, Scalar val)
    {
        gsVector<Scalar,3> result = gsVector<Scalar,3>::Zero();
        result[pos] = val;
        return result;
    }

    const gsMatrix<Scalar> & eval(const index_t k) const
    {
        theta  = _C.data().values[0].col(k);
        dtheta = _C.data().values[1].col(k);
        evec = _vec.eval(k);

        const gsMultiBasis<Scalar> & mbu = static_cast<const gsMultiBasis<Scalar>&>(_u.source());
        const gsMultiBasis<Scalar> & mbv = static_cast<const gsMultiBasis<Scalar>&>(_v.source());
        const index_t Au = mbu.basis(0).active(theta).rows();
        const index_t Av = mbv.basis(0).active(theta).rows();

        bGradsu = mbu.basis(0).deriv(theta);
        bGradsv = mbv.basis(0).deriv(theta);
        bHessu =  mbu.basis(0).deriv2(theta);
        bHessv =  mbv.basis(0).deriv2(theta);

        curve_normal_varalong_expr<Scalar> cnv_vara(_S,_C);
        cnvec_vara = cnv_vara.eval(k);
        const Scalar measure = cnvec_vara.norm();
        curve_normal_varalong_var1_expr<E1> cnv_vara_var1u(_u,_S,_defS,_C);  
        cnvec_vara_var1u = cnv_vara_var1u.eval(k).transpose();
        curve_normal_varalong_var1_expr<E2> cnv_vara_var1v(_v,_S,_defS,_C);
        cnvec_vara_var1v = cnv_vara_var1v.eval(k).transpose();

        matrix.resize(Au*_u.dim(), Av*_v.dim());
        res.resize(Au*_u.dim(), Av*_v.dim());

        for (index_t d = 0; d!= _u.dim(); ++d) // for all components
        {
            const short_t s = d*Au;
            for (index_t j = 0; j!= Au; ++j) // for all actives
            {
                secDerToHessian(bHessu.col(0).segment(3*j,3),2,localbHess);
                localbHess.resize(2,2);
                dotu = localbHess * dtheta;
                measure1u = (cnvec_vara_var1u.col3d(s+j).dot(cnvec_vara.col3d(0)))/measure;

                for (index_t e = 0; e!= _v.dim(); ++e) // for all components
                {
                    const short_t r = e*Av;
                    for (index_t k = 0; k!= Av; ++k) // for all actives
                    {
                        secDerToHessian(bHessv.col(0).segment(3*j,3),2,localbHess);
                        localbHess.resize(2,2);
                        dotv = localbHess * dtheta;
                        measure1v = (cnvec_vara_var1v.col3d(r+k).dot(cnvec_vara.col3d(0)))/measure;

                        //second-variation of derivative of normal basis vector along curve (not normalized)
                        matrix(s+j,r+k) = vecFun(d,dotu(0)).cross(vecFun(e,bGradsv.at(2*k+1))) +
                                          vecFun(e,dotv(0)).cross(vecFun(d,bGradsu.at(2*j+1))) +
                                          vecFun(d,bGradsu.at(2*j )).cross(vecFun(e,dotv(1))) +
                                          vecFun(e,bGradsv.at(2*k )).cross(vecFun(d,dotu(1)));

                        measure2 = ((matrix(s+j,r+k).dot(cnvec_vara.col3d(0)) + 
                                     cnvec_vara_var1u.col3d(s+j).dot(cnvec_vara_var1v.col3d(r+k)))*measure -
                                    (cnvec_vara_var1u.col3d(s+j).dot(cnvec_vara.col3d(0)))*measure1v)/std::pow(measure,2);

                        //second-variation of derivative of normal basis vector along curve (normalized) times a vector
                        res(s+j,r+k) = ((matrix(s+j,r+k)*measure + cnvec_vara_var1u.col3d(s+j)*measure1v -
                                         cnvec_vara_var1v.col3d(r+k)*measure1u - cnvec_vara*measure2)/std::pow(measure,2) -
                                         2*measure1v * (cnvec_vara_var1u.col3d(s+j)*measure - cnvec_vara*measure1u)/
                                         std::pow(measure,3)).dot(evec.col3d(0));
                    }
                }
            }
        }
        return res;
    }

    index_t rows() const { return 1; }

    index_t cols() const { return 1; }

    void parse(gsExprHelper<Scalar> & evList) const
    {
        evList.add(_C);
        _C.data().flags |= NEED_VALUE | NEED_DERIV;
        cnv_vara(_S,_C).parse(evList);
        cnv_vara_var1(_u,_S,_defS,_C).parse(evList); 
        cnv_vara_var1(_v,_S,_defS,_C).parse(evList);
        _vec.parse(evList);
    }

    const gsFeSpace<Scalar> & rowVar() const { return _u.rowVar(); }
    const gsFeSpace<Scalar> & colVar() const {return _u.rowVar();}
    index_t cardinality_impl() const { return _u.cardinality_impl(); }

    void print(std::ostream &os) const { os << "curve_normal_varalong_var2dot("; _u.print(os); os <<")"; }
};

template<class E1, class E2, class E3>
class curve_binormal_varalong_var2dot_expr : public _expr<curve_binormal_varalong_var2dot_expr<E1,E2,E3> >
{
public:
    typedef typename E1::Scalar Scalar;

private:
    typename E1::Nested_t _u;
    typename E2::Nested_t _v;
    typename gsGeometryMap<Scalar>::Nested_t _S;
    typename gsGeometryMap<Scalar>::Nested_t _defS;
    typename gsGeometryMap<Scalar>::Nested_t _C;
    typename E3::Nested_t _vec;

public:
    enum{ Space = 3, ScalarValued= 0, ColBlocks= 0};

    curve_binormal_varalong_var2dot_expr(const E1 &u, const E2 &v, const gsGeometryMap<Scalar> &S, const gsGeometryMap<Scalar> &defS,
                                         const gsGeometryMap<Scalar> &C, _expr<E3> const &vec): _u(u), _v(v), _S(S), _defS(defS), _C(C), _vec(vec)
                                         {}

    mutable gsMatrix<Scalar> theta, dtheta, bGradsu, bGradsv, bHessu, bHessv, localbHess, dotu, dotv, evec, res;
    mutable gsMatrix<Scalar> ctvec, ctvec_vara, ctvec_var1u, ctvec_var1v, ctvec_vara_var1u, ctvec_vara_var1v;
    mutable gsMatrix<Scalar> cnvec_var1u, cnvec_var1v, cnvec_vara_var1u, cnvec_vara_var1v, cbvec_vara_var1u;
    mutable gsVector<Scalar,3> cnvec_var2, cnvec_vara_var2;
    mutable gsMatrix<gsVector<Scalar,3>> matrix;
    mutable Scalar measure;

#   define Eigen gsEigen
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
#   undef Eigen

    // helper function
    static inline gsVector<Scalar,3> vecFun(index_t pos, Scalar val)
    {
        gsVector<Scalar,3> result = gsVector<Scalar,3>::Zero();
        result[pos] = val;
        return result;
    }

    const gsMatrix<Scalar> & eval(const index_t k) const
    {
        theta  = _C.data().values[0].col(k);
        dtheta = _C.data().values[1].col(k);
        evec = _vec.eval(k);

        const gsMultiBasis<Scalar> & mbu = static_cast<const gsMultiBasis<Scalar>&>(_u.source());
        const gsMultiBasis<Scalar> & mbv = static_cast<const gsMultiBasis<Scalar>&>(_v.source());
        const index_t Au = mbu.basis(0).active(theta).rows(); 
        const index_t Av = mbv.basis(0).active(theta).rows();

        bGradsu = mbu.basis(0).deriv(theta);
        bGradsv = mbv.basis(0).deriv(theta);
        bHessu =  mbu.basis(0).deriv2(theta);
        bHessv =  mbv.basis(0).deriv2(theta);

        curve_tangent_expr<Scalar> ctv(_S,_C);
        ctvec = ctv.eval(k);
        curve_tangent_varalong_expr<Scalar> ctv_vara(_S,_C);
        ctvec_vara = ctv_vara.eval(k);
        curve_tangent_var1_expr<E1> ctv_var1u(_u,_S,_C);    
        ctvec_var1u = ctv_var1u.eval(k).transpose();
        curve_tangent_var1_expr<E2> ctv_var1v(_v,_S,_C);    
        ctvec_var1v = ctv_var1v.eval(k).transpose();
        curve_tangent_varalong_var1_expr<E1> ctv_vara_var1u(_u,_S,_C);
        ctvec_vara_var1u = ctv_vara_var1u.eval(k).transpose();
        curve_tangent_varalong_var1_expr<E2> ctv_vara_var1v(_v,_S,_C);
        ctvec_vara_var1v = ctv_vara_var1v.eval(k).transpose();
        curve_normal_var1_expr<E1> cnv_var1u(_u,_S,_C);
        cnvec_var1u = cnv_var1u.eval(k).transpose();
        curve_normal_var1_expr<E2> cnv_var1v(_v,_S,_C);
        cnvec_var1v = cnv_var1v.eval(k).transpose();
        curve_normal_varalong_var1_expr<E1> cnv_vara_var1u(_u,_S,_defS,_C);
        cnvec_vara_var1u = cnv_vara_var1u.eval(k).transpose();
        curve_normal_varalong_var1_expr<E2> cnv_vara_var1v(_v,_S,_defS,_C);    
        cnvec_vara_var1v = cnv_vara_var1v.eval(k).transpose();
        curve_binormal_varalong_var1_expr<E1> cbv_vara_var1u(_u,_S,_defS,_C);
        cbvec_vara_var1u = cbv_vara_var1u.eval(k).transpose();

        matrix.resize(Au*_u.dim(), Av*_v.dim());
        res.resize(Au*_u.dim(), Av*_v.dim());

        for (index_t d = 0; d!= _u.dim(); ++d) // for all components
        {
            const short_t s = d*Au;
            for (index_t j = 0; j!= Au; ++j) // for all actives
            {
                secDerToHessian(bHessu.col(0).segment(3*j,3),2,localbHess);
                localbHess.resize(2,2);
                dotu = localbHess * dtheta;
                measure = cbvec_vara_var1u.col3d(s+j).norm();

                for (index_t e = 0; e!= _v.dim(); ++e) // for all components
                {
                    const short_t r = e*Av;
                    for (index_t k = 0; k!= Av; ++k) // for all actives
                    {
                        secDerToHessian(bHessv.col(0).segment(3*j,3),2,localbHess);
                        localbHess.resize(2,2);
                        dotv = localbHess * dtheta;

                        cnvec_var2 = vecFun(d,bGradsu.at(2*j )).cross(vecFun(e,bGradsv.at(2*k+1))) +
                                     vecFun(e,bGradsv.at(2*k )).cross(vecFun(d,bGradsu.at(2*j+1)));
                        
                        cnvec_vara_var2 = vecFun(d,dotu(0)).cross(vecFun(e,bGradsv.at(2*k+1))) +
                                          vecFun(e,dotv(0)).cross(vecFun(d,bGradsu.at(2*j+1))) +
                                          vecFun(d,bGradsu.at(2*j )).cross(vecFun(e,dotv(1))) +
                                          vecFun(e,bGradsv.at(2*k )).cross(vecFun(d,dotu(1)));
                        
                        //second-variation of derivative of binormal basis vector along curve (not normalized)
                        matrix(s+j,r+k) = ctvec_vara_var1u.col3d(s+j).cross(cnvec_var1v.col3d(r+k)) +
                                          ctvec_vara_var1v.col3d(r+k).cross(cnvec_var1u.col3d(s+j)) +
                                          ctvec_vara.col3d(0).cross(cnvec_var2) + ctvec.col3d(0).cross(cnvec_vara_var2) +
                                          ctvec_var1u.col3d(s+j).cross(cnvec_vara_var1v.col3d(r+k)) +
                                          ctvec_var1v.col3d(r+k).cross(cnvec_vara_var1u.col3d(s+j));
                    
                        //second-variation of derivative of binormal basis vector along curve (normalized) times a vector
                        res(s+j,r+k) = ((matrix(s+j,r+k)/measure) - (cbvec_vara_var1u.col3d(s+j)/std::pow(measure,3)) *
                                        (matrix(s+j,r+k).dot(cbvec_vara_var1u.col3d(s+j)))).dot(evec.col(0));
                    }
                }
            }
        }
        return res;
    }

    index_t rows() const { return 1; }

    index_t cols() const { return 1; }

    void parse(gsExprHelper<Scalar> & evList) const
    {
        evList.add(_C);
        _C.data().flags |= NEED_VALUE | NEED_DERIV;
        ctv(_S,_C).parse(evList);
        ctv_vara(_S,_C).parse(evList);
        ctv_var1(_u,_S,_C).parse(evList); 
        ctv_var1(_v,_S,_C).parse(evList);  
        ctv_vara_var1(_u,_S,_C).parse(evList);  
        ctv_vara_var1(_v,_S,_C).parse(evList);     
        cnv_var1(_u,_S,_C).parse(evList);   
        cnv_var1(_v,_S,_C).parse(evList);   
        cnv_vara_var1(_u,_S,_defS,_C).parse(evList); 
        cnv_vara_var1(_v,_S,_defS,_C).parse(evList);  
        cbv_vara_var1(_u,_S,_defS,_C).parse(evList);
        _vec.parse(evList);
    }

    const gsFeSpace<Scalar> & rowVar() const { return _u.rowVar(); }
    const gsFeSpace<Scalar> & colVar() const {return _v.rowVar();}
    index_t cardinality_impl() const { return _u.cardinality_impl(); }

    void print(std::ostream &os) const { os << "curve_binormal_varalong_var2dot("; _u.print(os); os <<")"; }
};

template<class E1, class E2, class E3, class E4>
class var1_dot_othervar1_expr : public _expr<var1_dot_othervar1_expr<E1,E2,E3,E4>>
{
public:
    typedef typename E1::Scalar Scalar;

private:
    typename E1::Nested_t _u;
    typename E2::Nested_t _v;
    typename gsGeometryMap<Scalar>::Nested_t _C;
    typename E3::Nested_t _that;
    typename E4::Nested_t _other;
    
public:
    enum{ Space = 3, ScalarValued= 0, ColBlocks= 0};

    var1_dot_othervar1_expr(const E1 &u, const E2 &v, const gsGeometryMap<Scalar> &C, _expr<E3> const &that, _expr<E4> const &other): 
                            _u(u), _v(v), _C(C), _that(that), _other(other){}

    mutable gsMatrix<Scalar> ethat, eother, theta, res;

#   define Eigen gsEigen
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
#   undef Eigen

    const gsMatrix<Scalar> & eval(const index_t k) const
    {
        theta  = _C.data().values[0].col(k);
        ethat = _that.eval(k).transpose();
        eother = _other.eval(k).transpose();

        const gsMultiBasis<Scalar> & mbu = static_cast<const gsMultiBasis<Scalar>&>(_u.source());
        const gsMultiBasis<Scalar> & mbv = static_cast<const gsMultiBasis<Scalar>&>(_v.source());
        const index_t Au = mbu.basis(0).active(theta).rows();
        const index_t Av = mbv.basis(0).active(theta).rows();

        res.resize(Au*_u.dim(), Av*_v.dim());

        for (index_t i = 0; i!= Au*_u.dim(); ++i)
        {
            for (index_t j = 0; j!= Av*_v.dim(); ++j)
            {
                res(i,j) = ethat.col3d(i).dot(eother.col3d(j));
            }
        }
        return res;
    }

    index_t rows() const { return 1; }
    index_t cols() const { return 1; }

    void parse(gsExprHelper<Scalar> & evList) const
    {
         evList.add(_C);
         _C.data().flags |= NEED_VALUE;
         _that.parse(evList);
         _other.parse(evList);
    }

    const gsFeSpace<Scalar> & rowVar() const { return _u.rowVar(); }
    const gsFeSpace<Scalar> & colVar() const {return _u.rowVar();}
    index_t cardinality_impl() const { return _u.cardinality_impl(); }

    void print(std::ostream &os) const { os << "var1_dot_othervar1("; _u.print(os); os <<")"; }
};

template<class T> EIGEN_STRONG_INLINE
curve_tangent_expr<T> ctv(const gsGeometryMap<T> &S,
                          const gsGeometryMap<T> &C) { return curve_tangent_expr<T>(S,C); }

template<class E> EIGEN_STRONG_INLINE
curve_tangent_var1_expr<E> ctv_var1(const E &u, const gsGeometryMap<typename E::Scalar> &S,
                           const gsGeometryMap<typename E::Scalar> &C) { return curve_tangent_var1_expr<E>(u,S,C); }

template<class T> EIGEN_STRONG_INLINE
curve_tangent_varalong_expr<T> ctv_vara(const gsGeometryMap<T> &S,
                               const gsGeometryMap<T> &C) { return curve_tangent_varalong_expr<T>(S,C); }

template<class E> EIGEN_STRONG_INLINE
curve_tangent_varalong_var1_expr<E> ctv_vara_var1(const E &u, const gsGeometryMap<typename E::Scalar> &S,
                                    const gsGeometryMap<typename E::Scalar> &C) { return curve_tangent_varalong_var1_expr<E>(u,S,C); }

template<class T> EIGEN_STRONG_INLINE
curve_normal_expr<T>  cnv(const gsGeometryMap<T> &S,
                          const gsGeometryMap<T> &C) { return curve_normal_expr<T>(S,C); }

template<class E> EIGEN_STRONG_INLINE
curve_normal_var1_expr<E> cnv_var1(const E &u, const gsGeometryMap<typename E::Scalar> &S,
                          const gsGeometryMap<typename E::Scalar> &C) { return curve_normal_var1_expr<E>(u,S,C); }

template<class T> EIGEN_STRONG_INLINE
curve_normal_varalong_expr<T>  cnv_vara(const gsGeometryMap<T> &S,
                               const gsGeometryMap<T> &C) { return curve_normal_varalong_expr<T>(S,C); }

template<class E> EIGEN_STRONG_INLINE
curve_normal_varalong_var1_expr<E> cnv_vara_var1(const E &u, const gsGeometryMap<typename E::Scalar> &S,
                                   const gsGeometryMap<typename E::Scalar> &defS,
                                   const gsGeometryMap<typename E::Scalar> &C)
                                   { return curve_normal_varalong_var1_expr<E>(u,S,defS,C); }

template<class E> EIGEN_STRONG_INLINE
curve_normal_varalong_var1_normalized_expr<E> cnv_vara_var1_normalized(const E &u, const gsGeometryMap<typename E::Scalar> &S,
                                              const gsGeometryMap<typename E::Scalar> &defS,
                                              const gsGeometryMap<typename E::Scalar> &C)
                                              { return curve_normal_varalong_var1_normalized_expr<E>(u,S,defS,C); }

template<class E1, class E2, class E3> EIGEN_STRONG_INLINE
curve_normal_var2dot_expr<E1,E2,E3> cnv_var2dot(const E1 &u, const E2 &v, const gsGeometryMap<typename E1::Scalar> &C, const E3 &vec)
                                               { return curve_normal_var2dot_expr<E1,E2,E3>(u,v,C,vec); }

template<class E1, class E2, class E3> EIGEN_STRONG_INLINE
curve_normal_varalong_var2dot_expr<E1,E2,E3> cnv_vara_var2dot(const E1 &u, const E2 &v, const gsGeometryMap<typename E1::Scalar> &S,
                                             const gsGeometryMap<typename E1::Scalar> &defS, const gsGeometryMap<typename E1::Scalar> &C,
                                             const E3 &vec) { return curve_normal_varalong_var2dot_expr<E1,E2,E3>(u,v,S,defS,C,vec); }

template<class T> EIGEN_STRONG_INLINE
curve_binormal_expr<T>  cbv(const gsGeometryMap<T> &S,
                            const gsGeometryMap<T> &C) { return curve_binormal_expr<T>(S,C); }

template<class E> EIGEN_STRONG_INLINE
curve_binormal_var1_expr<E> cbv_var1(const E &u, const gsGeometryMap<typename E::Scalar> &S,
                            const gsGeometryMap<typename E::Scalar> &C) { return curve_binormal_var1_expr<E>(u,S,C); }

template<class T> EIGEN_STRONG_INLINE
curve_binormal_varalong_expr<T> cbv_vara(const gsGeometryMap<T> &S,
                                const gsGeometryMap<T> &C) { return curve_binormal_varalong_expr<T>(S,C); }

template<class E> EIGEN_STRONG_INLINE
curve_binormal_varalong_var1_expr<E> cbv_vara_var1(const E &u, const gsGeometryMap<typename E::Scalar> &S,
                                     const gsGeometryMap<typename E::Scalar> &defS, const gsGeometryMap<typename E::Scalar> &C)
                                     { return curve_binormal_varalong_var1_expr<E>(u,S,defS,C); }

template<class E> EIGEN_STRONG_INLINE
curve_binormal_varalong_var1_normalized_expr<E> cbv_vara_var1_normalized(const E &u, const gsGeometryMap<typename E::Scalar> &S,
                                                const gsGeometryMap<typename E::Scalar> &defS, const gsGeometryMap<typename E::Scalar> &C)
                                                { return curve_binormal_varalong_var1_normalized_expr<E>(u,S,defS,C); }

template<class E1, class E2, class E3> EIGEN_STRONG_INLINE
curve_binormal_var2dot_expr<E1,E2,E3> cbv_var2dot(const E1 &u, const E2 &v, const gsGeometryMap<typename E1::Scalar> &S,
                                                  const gsGeometryMap<typename E1::Scalar> &C, const E3 &vec) 
                                                 { return curve_binormal_var2dot_expr<E1,E2,E3>(u,v,S,C,vec); }

template<class E1, class E2, class E3> EIGEN_STRONG_INLINE
curve_binormal_varalong_var2dot_expr<E1,E2,E3> cbv_vara_var2dot(const E1 &u, const E2 &v, const gsGeometryMap<typename E1::Scalar> &S,
                                               const gsGeometryMap<typename E1::Scalar> &defS, const gsGeometryMap<typename E1::Scalar> &C,
                                               const E3 &vec) { return curve_binormal_varalong_var2dot_expr<E1,E2,E3>(u,v,S,defS,C,vec); }

template<class E1, class E2, class E3, class E4> EIGEN_STRONG_INLINE
var1_dot_othervar1_expr<E1,E2,E3,E4> var1_dot_var1(const E1 &u, const E2 &v, const gsGeometryMap<typename E1::Scalar> &C,
                                                   const E3 &that, const E4 &other)
                                                  {return var1_dot_othervar1_expr<E1,E2,E3,E4>(u,v,C,that,other);}
}
}

using namespace gismo;

template <class T>
T curveCoordinate(const gsGeometry<T> & geometry, const T & value, const short_t & dir);

template <class T = real_t>
class G_fun : public gsFunction<T>
{
public:
    // Constructor: curve geometry map
    G_fun(const gsMultiPatch<T> & mp_curve,
          const gsVector<T> & gamma)
    :
    m_curve(mp_curve),
    m_gamma(gamma)
    {
        GISMO_ASSERT(m_curve.nPatches()==1, "The curve must be a single patch");
        GISMO_ASSERT(m_gamma.rows()==1, "The gamma vector must be scalar coordinates");
    }

    short_t domainDim() const override { return m_curve.patch(0).coefs().rows()*m_curve.patch(0).coefs().cols(); }

    short_t targetDim() const override { return 1; } //used when calling FDM-based deriv_into

    void eval_into(const gsMatrix<T>  & u, gsMatrix<T> &result) const override
    {
        result.resize(m_curve.patch(0).geoDim(),u.cols());

        for (index_t k = 0; k!= u.cols(); ++k)
        {
            // Perturb the curve with control point displacements
            gsMultiPatch<T> curve_def = m_curve;
            curve_def.patch(0).coefs() += u.reshapeCol(k,m_curve.patch(0).coefs().rows(),m_curve.patch(0).coefs().cols());

            // Construct expression evaluator
            gsExprEvaluator<T> exprEvaluator;
            auto G_curve = exprEvaluator.getMap(curve_def);

            // Evaluate curve geometry map expression at point gamma
            result.col(k) = exprEvaluator.eval(G_curve, m_gamma);
        }
    }

protected:

    const gsMultiPatch<T> & m_curve;
    const gsVector<T>     & m_gamma;

};

template <class T = real_t>
class cnv_fun : public gsFunction<T>
{
public:
    // Constructor: normal basis vector
    cnv_fun(const gsMultiPatch<T> & mp_surf,
            const gsMultiPatch<T> & mp_curve,
            const gsVector<T> & gamma)
    :
    m_surf(mp_surf),
    m_curve(mp_curve),
    m_gamma(gamma)
    {
        GISMO_ASSERT(m_surf.nPatches()==1, "The surface must be a single patch");
        GISMO_ASSERT(m_curve.nPatches()==1, "The curve must be a single patch");
        GISMO_ASSERT(m_gamma.rows()==1, "The gamma vector must be  scalar coordinates");
    }

    short_t domainDim() const override { return m_surf.patch(0).coefs().rows()*m_surf.patch(0).coefs().cols(); }

    short_t targetDim() const override { return 3; } //used when calling FDM-based deriv_into

    void eval_into(const gsMatrix<T>  & u, gsMatrix<T> &result) const override
    {
        result.resize(m_surf.patch(0).geoDim(),u.cols());

        for (index_t k = 0; k!= u.cols(); ++k)
        {
            // Perturb the surface with control point displacements
            gsMultiPatch<T> surf_def = m_surf;
            surf_def.patch(0).coefs() += u.reshapeCol(k,m_surf.patch(0).coefs().rows(),m_surf.patch(0).coefs().cols());

            // Construct expression evaluator
            gsExprEvaluator<T> exprEvaluator;
            auto G_curve = exprEvaluator.getMap(m_curve);
            auto G_surf  = exprEvaluator.getMap(surf_def);

            // Evaluate normal basis vector expression at point gamma
            result.col(k) = exprEvaluator.eval( gismo::expr::cnv(G_surf, G_curve), m_gamma);
        }
    }

protected:

    const gsMultiPatch<T> & m_surf;
    const gsMultiPatch<T> & m_curve;
    const gsVector<T>     & m_gamma;

};

template <class T = real_t>
class ctv_fun : public gsFunction<T>
{
public:
    // Constructor: normal basis vector
    ctv_fun(const gsMultiPatch<T> & mp_surf,
            const gsMultiPatch<T> & mp_curve,
            const gsVector<T> & gamma)
    :
    m_surf(mp_surf),
    m_curve(mp_curve),
    m_gamma(gamma)
    {
        GISMO_ASSERT(m_surf.nPatches()==1, "The surface must be a single patch");
        GISMO_ASSERT(m_curve.nPatches()==1, "The curve must be a single patch");
        GISMO_ASSERT(m_gamma.rows()==1, "The gamma vector must be scalar coordinates");
    }

    short_t domainDim() const override { return m_surf.patch(0).coefs().rows()*m_surf.patch(0).coefs().cols(); }

    short_t targetDim() const override { return 3; } //used when calling FDM-based deriv_into

    void eval_into(const gsMatrix<T>  & u, gsMatrix<T> &result) const override
    {
        result.resize(m_surf.patch(0).geoDim(),u.cols());

        for (index_t k = 0; k!= u.cols(); ++k)
        {
            // Perturb the surface with control point displacements
            gsMultiPatch<T> surf_def = m_surf;
            surf_def.patch(0).coefs() += u.reshapeCol(k,m_surf.patch(0).coefs().rows(),m_surf.patch(0).coefs().cols());

            // Construct expression evaluator
            gsExprEvaluator<T> exprEvaluator;
            auto G_curve = exprEvaluator.getMap(m_curve);
            auto G_surf  = exprEvaluator.getMap(surf_def);

            // Evaluate tangent basis vector expression at point gamma
            result.col(k) = exprEvaluator.eval( gismo::expr::ctv(G_surf, G_curve), m_gamma);
        }
    }

protected:

    const gsMultiPatch<T> & m_surf;
    const gsMultiPatch<T> & m_curve;
    const gsVector<T>     & m_gamma;

};

template <class T = real_t>
class cbv_fun : public gsFunction<T>
{
public:
    // Constructor: binormal basis vector
    cbv_fun(const gsMultiPatch<T> & mp_surf,
            const gsMultiPatch<T> & mp_curve,
            const gsVector<T> & gamma)
    :
    m_surf(mp_surf),
    m_curve(mp_curve),
    m_gamma(gamma)
    {
        GISMO_ASSERT(m_surf.nPatches()==1, "The surface must be a single patch");
        GISMO_ASSERT(m_curve.nPatches()==1, "The curve must be a single patch");
        GISMO_ASSERT(m_gamma.rows()==1, "The gamma vector must be  scalar coordinates");
    }

    short_t domainDim() const override { return m_surf.patch(0).coefs().rows()*m_surf.patch(0).coefs().cols(); }

    short_t targetDim() const override { return 3; } //used when calling FDM-based deriv_into

    void eval_into(const gsMatrix<T>  & u, gsMatrix<T> &result) const override
    {
        result.resize(m_surf.patch(0).geoDim(),u.cols());

        for (index_t k = 0; k!= u.cols(); ++k)
        {
            // Perturb the surface with control point displacements
            gsMultiPatch<T> surf_def = m_surf;
            surf_def.patch(0).coefs() += u.reshapeCol(k,m_surf.patch(0).coefs().rows(),m_surf.patch(0).coefs().cols());

            // Construct expression evaluator
            gsExprEvaluator<T> exprEvaluator;
            auto G_curve = exprEvaluator.getMap(m_curve);
            auto G_surf  = exprEvaluator.getMap(surf_def);
            
            // Evaluate binormal basis vector expression at point gamma
            result.col(k) = exprEvaluator.eval( gismo::expr::cbv(G_surf, G_curve), m_gamma);
        }
    }

protected:

    const gsMultiPatch<T> & m_surf;
    const gsMultiPatch<T> & m_curve;
    const gsVector<T>     & m_gamma;

};

template <class T = real_t>
class ctv_along_fun : public gsFunction<T>
{
public:
    ctv_along_fun(const gsMultiPatch<T> & mp_surf, // deformed surface
                  const gsMultiPatch<T> & mp_curve)
    :
    m_surf(mp_surf),
    m_curve(mp_curve)
    {
        GISMO_ASSERT(m_surf.nPatches()==1, "The surface must be a single patch");
        GISMO_ASSERT(m_curve.nPatches()==1, "The curve must be a single patch");
    }

    short_t domainDim() const override { return 1; }

    short_t targetDim() const override { return 3; } //used when calling FDM-based deriv_into

    // u is a matrix of gamma values per point (in each column), 1xu.cols()
    void eval_into(const gsMatrix<T>  & u, gsMatrix<T> &result) const override
    {
        result.resize(m_surf.patch(0).geoDim(),u.cols());

        for (index_t k = 0; k!= u.cols(); ++k)
        {
            // Construct expression evaluator
            gsExprEvaluator<T> exprEvaluator;
            auto G_curve = exprEvaluator.getMap(m_curve);
            auto G_surf  = exprEvaluator.getMap(m_surf);

            // Evaluate tangent basis vector expression at point gamma
            result.col(k) = exprEvaluator.eval( gismo::expr::ctv(G_surf, G_curve), u.col(k));
        }
    }

protected:

    const gsMultiPatch<T> & m_surf;
    const gsMultiPatch<T> & m_curve;

};

template <class T = real_t>
class cnv_along_fun : public gsFunction<T>
{
public:
    cnv_along_fun(const gsMultiPatch<T> & mp_surf,
                  const gsMultiPatch<T> & mp_curve)
    :
    m_surf(mp_surf),
    m_curve(mp_curve)
    {
        GISMO_ASSERT(m_surf.nPatches()==1, "The surface must be a single patch");
        GISMO_ASSERT(m_curve.nPatches()==1, "The curve must be a single patch");
    }

    short_t domainDim() const override { return 1; }

    short_t targetDim() const override { return 3; }

    void eval_into(const gsMatrix<T>  & u, gsMatrix<T> &result) const override
    {
        result.resize(m_surf.patch(0).geoDim(),u.cols());

        for (index_t k = 0; k!= u.cols(); ++k)
        {
            // Construct expression evaluator
            gsExprEvaluator<T> exprEvaluator;
            auto G_curve = exprEvaluator.getMap(m_curve);
            auto G_surf  = exprEvaluator.getMap(m_surf);

            // Evaluate tangent basis vector expression at point gamma
            result.col(k) = exprEvaluator.eval( gismo::expr::cnv(G_surf, G_curve), u.col(k));
        }
    }

protected:

    const gsMultiPatch<T> & m_surf;
    const gsMultiPatch<T> & m_curve;

};

template <class T = real_t>
class cbv_along_fun : public gsFunction<T>
{
public:
    cbv_along_fun(const gsMultiPatch<T> & mp_surf,
                  const gsMultiPatch<T> & mp_curve)
    :
    m_surf(mp_surf),
    m_curve(mp_curve)
    {
        GISMO_ASSERT(m_surf.nPatches()==1, "The surface must be a single patch");
        GISMO_ASSERT(m_curve.nPatches()==1, "The curve must be a single patch");
    }

    short_t domainDim() const override { return 1; }

    short_t targetDim() const override { return 3; }

    void eval_into(const gsMatrix<T>  & u, gsMatrix<T> &result) const override
    {
        result.resize(m_surf.patch(0).geoDim(),u.cols());

        for (index_t k = 0; k!= u.cols(); ++k)
        {
            // Construct expression evaluator
            gsExprEvaluator<T> exprEvaluator;
            auto G_curve = exprEvaluator.getMap(m_curve);
            auto G_surf  = exprEvaluator.getMap(m_surf);

            // Evaluate tangent basis vector expression at point gamma
            result.col(k) = exprEvaluator.eval( gismo::expr::cbv(G_surf, G_curve), u.col(k));
        }
    }

protected:

    const gsMultiPatch<T> & m_surf;
    const gsMultiPatch<T> & m_curve;

};

template <class T = real_t>
class ctv_along_var_fun : public gsFunction<T>
{
public:
    // Constructor: normal basis vector
    ctv_along_var_fun(const gsMultiPatch<T> & mp_surf,
                   const gsMultiPatch<T> & mp_curve,
                   const gsVector<T> & gamma)
    :
    m_surf(mp_surf),
    m_curve(mp_curve),
    m_gamma(gamma)
    {
        GISMO_ASSERT(m_surf.nPatches()==1, "The surface must be a single patch");
        GISMO_ASSERT(m_curve.nPatches()==1, "The curve must be a single patch");
        GISMO_ASSERT(m_gamma.rows()==1, "The gamma vector must be scalar coordinates");
    }

    short_t domainDim() const override { return m_surf.patch(0).coefs().rows()*m_surf.patch(0).coefs().cols(); }

    short_t targetDim() const override { return 3; } //used when calling FDM-based deriv_into

    void eval_into(const gsMatrix<T>  & u, gsMatrix<T> &result) const override
    {
        result.resize(m_surf.patch(0).geoDim(),u.cols());

        for (index_t k = 0; k!= u.cols(); ++k)
        {
            // Perturb the surface with control point displacements
            gsMultiPatch<T> surf_def = m_surf;
            surf_def.patch(0).coefs() += u.reshapeCol(k,m_surf.patch(0).coefs().rows(),m_surf.patch(0).coefs().cols());

            // Construct expression evaluator
            gsExprEvaluator<T> exprEvaluator;
            auto G_curve = exprEvaluator.getMap(m_curve);
            auto G_surf  = exprEvaluator.getMap(surf_def);

            // Evaluate tangent basis vector expression along curve at point gamma
            result.col(k) = exprEvaluator.eval( gismo::expr::ctv_vara(G_surf, G_curve), m_gamma);
        }
    }

protected:

    const gsMultiPatch<T> & m_surf;
    const gsMultiPatch<T> & m_curve;
    const gsVector<T>     & m_gamma;

};

template <class T = real_t>
class cnv_along_var_fun : public gsFunction<T>
{
public:
    // Constructor: normal basis vector
    cnv_along_var_fun(const gsMultiPatch<T> & mp_surf,
                      const gsMultiPatch<T> & mp_curve,
                      const gsVector<T> & gamma)
    :
    m_surf(mp_surf),
    m_curve(mp_curve),
    m_gamma(gamma)
    {
        GISMO_ASSERT(m_surf.nPatches()==1, "The surface must be a single patch");
        GISMO_ASSERT(m_curve.nPatches()==1, "The curve must be a single patch");
        GISMO_ASSERT(m_gamma.rows()==1, "The gamma vector must be scalar coordinates");
    }

    short_t domainDim() const override { return m_surf.patch(0).coefs().rows()*m_surf.patch(0).coefs().cols(); }

    short_t targetDim() const override { return 3; } //used when calling FDM-based deriv_into

    void eval_into(const gsMatrix<T>  & u, gsMatrix<T> &result) const override
    {
        result.resize(m_surf.patch(0).geoDim(),u.cols());

        for (index_t k = 0; k!= u.cols(); ++k)
        {
            // Perturb the surface with control point displacements
            gsMultiPatch<T> surf_def = m_surf;
            surf_def.patch(0).coefs() += u.reshapeCol(k,m_surf.patch(0).coefs().rows(),m_surf.patch(0).coefs().cols());

            // Construct expression evaluator
            gsExprEvaluator<T> exprEvaluator;
            auto G_curve = exprEvaluator.getMap(m_curve);
            auto G_surf  = exprEvaluator.getMap(surf_def);

            // Evaluate normal basis vector expression along curve at point gamma
            result.col(k) = exprEvaluator.eval( gismo::expr::cnv_vara(G_surf, G_curve), m_gamma);
        }
    }

protected:

    const gsMultiPatch<T> & m_surf;
    const gsMultiPatch<T> & m_curve;
    const gsVector<T>     & m_gamma;

};

template <class T = real_t>
class cbv_along_var_fun : public gsFunction<T>
{
public:
    // Constructor: normal basis vector
    cbv_along_var_fun(const gsMultiPatch<T> & mp_surf,
                      const gsMultiPatch<T> & mp_curve,
                      const gsVector<T> & gamma)
    :
    m_surf(mp_surf),
    m_curve(mp_curve),
    m_gamma(gamma)
    {
        GISMO_ASSERT(m_surf.nPatches()==1, "The surface must be a single patch");
        GISMO_ASSERT(m_curve.nPatches()==1, "The curve must be a single patch");
        GISMO_ASSERT(m_gamma.rows()==1, "The gamma vector must be scalar coordinates");
    }

    short_t domainDim() const override { return m_surf.patch(0).coefs().rows()*m_surf.patch(0).coefs().cols(); }

    short_t targetDim() const override { return 3; } //used when calling FDM-based deriv_into

    void eval_into(const gsMatrix<T>  & u, gsMatrix<T> &result) const override
    {
        result.resize(m_surf.patch(0).geoDim(),u.cols());

        for (index_t k = 0; k!= u.cols(); ++k)
        {
            // Perturb the surface with control point displacements
            gsMultiPatch<T> surf_def = m_surf;
            surf_def.patch(0).coefs() += u.reshapeCol(k,m_surf.patch(0).coefs().rows(),m_surf.patch(0).coefs().cols());

            // Construct expression evaluator
            gsExprEvaluator<T> exprEvaluator;
            auto G_curve = exprEvaluator.getMap(m_curve);
            auto G_surf  = exprEvaluator.getMap(surf_def);

            // Evaluate normal basis vector expression along curve at point gamma
            result.col(k) = exprEvaluator.eval( gismo::expr::cbv_vara(G_surf, G_curve), m_gamma);
        }
    }

protected:

    const gsMultiPatch<T> & m_surf;
    const gsMultiPatch<T> & m_curve;
    const gsVector<T>     & m_gamma;

};

int main(int argc, char *argv[])
{

    gsInfo<<"To do list:\n";
    gsInfo<<"- [ ] Parse expressions as compositions, i.e. pre-compute the coordinates of the curve and pre-compute the surface at this point.\n";
    gsInfo<<"- [ ] Parse dependent expressions (e.g., the use of curve_normal_expr in curve_tangent_expr in the curve_binormal_expr)\n";
    gsInfo<<"- [ ] Make the expressions work for 2D and 3D surfaces using enable_if\n";

    // bool plot = false; // If set to true, paraview file is generated and launched on exit
    // bool trim = false; // If set to true, trim/merge operations are displayed
    // bool intersect = false; // If set to true, intersection example is displayed

    bool surface = false;
    // std::string filename;

    gsCmdLine cmd("TODO");
    // cmd.addSwitch("plot", "Plot result in ParaView format", plot);
    // cmd.addSwitch("trim", "Basic trim/merge operations", trim);
    cmd.addSwitch("surface", "Model the geometry as a surface", surface);
    // cmd.addString("file", "Path to the file to be read", filename);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    // Make a BSpline surface
    gsKnotVector<> kv_u(0,1,2,3);
    gsKnotVector<> kv_v(0,1,2,3);
    gsTensorBSplineBasis<2, real_t> basis_s(kv_u, kv_v);
    gsMatrix<> coefs_s (basis_s.size(), 3);
    coefs_s << 0, 0, 0,
               0.166, 0, 0.130,
               0.5, 0, 0.215,
               0.833, 0, 0.130,
               1, 0, 0,
               0, 0.166, 0.130,
               0.166, 0.166, 0.260,
               0.5, 0.166, 0.347,
               0.833, 0.166, 0.260,
               1, 0.166, 0.130,
               0, 0.5, 0.215,
               0.166, 0.5, 0.347,
               0.5, 0.5, 0.433,
               0.833, 0.5, 0.347,
               1, 0.5, 0.215,
               0, 0.833, 0.130,
               0.166, 0.833, 0.260,
               0.5, 0.833, 0.347,
               0.833, 0.833, 0.260,
               1, 0.833, 0.130,
               0, 1, 0,
               0.166, 1, 0.130,
               0.5, 1, 0.215,
               0.833, 1, 0.130,
               1, 1, 0;

    gsTensorBSpline<2, real_t>  surf(basis_s, coefs_s); //original surface

    // Make a BSpline curve within the surface source domain
    gsKnotVector<> kv_c(0, 1, 2, 3); //start,end,interior knots, start/end multiplicity
    gsMatrix<> coefs_c(5, 2); //u,v;..
    coefs_c << 0, 0,
             0.2618, 0.053,
             0.5, 0.5,
             0.738, 0.9465,
             1,1;

    coefs_c.col(0) = coefs_c.col(0) * kv_u.last();
    coefs_c.col(1) = coefs_c.col(1) * kv_v.last();

    gsBSpline<> curve( kv_c, give(coefs_c));

    // 2D->2D
    gsMultiPatch<> mp_surf;
    mp_surf.addPatch(surf);
    if (surface)
        mp_surf.embed(3);

    gsInfo<<"The surface (R^"<<mp_surf.patch(0).domainDim()<<" -> R^"<<mp_surf.patch(0).targetDim()<<") is:\n"<<mp_surf.patch(0)<<"\n";
    // 1D->2D
    gsMultiPatch<> mp_curve;
    mp_curve.addPatch(curve);
    gsInfo<<"The curve (R^"<<mp_curve.patch(0).domainDim()<<" -> R^"<<mp_curve.patch(0).targetDim()<<") is:\n"<<mp_curve.patch(0)<<"\n";

    gsMultiBasis<> basis_curve(mp_curve);
    gsMultiBasis<> basis_surf(mp_surf);

    // Declare the expression assembler
    gsExprAssembler<> exprAssembler(1,1);
    // Set integration elements
    exprAssembler.setIntegrationElements(basis_curve); // maybe surf?
    // Register expressions inside assembler
    auto G_curve = exprAssembler.getMap(mp_curve);
    auto G_surf  = exprAssembler.getMap(mp_surf);
    auto defG_surf = exprAssembler.getMap(mp_surf); //NOTE: deformed surface must be set as solution field
    auto u_curve = exprAssembler.getSpace(basis_curve); // NOTE: needed only for validation
    auto u_surf  = exprAssembler.getSpace(basis_surf,mp_surf.geoDim());

    // Declare expression evaluator from assembler
    gsExprEvaluator<> exprEvaluator(exprAssembler);

    gsInfo<<">>>>>>>>>>>>>>>>>>>>>>>>>>> VALIDATION >>>>>>>>>>>>>>>>>>>>>>>>>>>>\n";

    gsVector<> gamma(1);
    gamma << 0.5;

    //control displacement vector for numerical derivation via central finite difference
    gsMatrix<> u(mp_curve.patch(0).coefs().rows()*mp_curve.patch(0).coefs().cols(),1);
    u.setZero();

    gsInfo<<"Curve geometry map at gamma\n";
    G_fun<real_t> GFun(mp_curve,gamma);
    gsDebug<<GFun.eval(u).transpose()<<"\n should match \n";
    gsDebug<<exprEvaluator.eval( G_curve, gamma ).transpose() << "\n";
    gsInfo<<"First variation of curve geometry map at gamma\n";
    gsDebug<<GFun.deriv(u).transpose()<<"\n should match \n";
    gsDebug<<exprEvaluator.eval( u_curve, gamma ).transpose() << "\n";

    //control displacement vector for numerical derivation via central finite difference
    u.clear();
    u.resize(mp_surf.patch(0).coefs().rows()*mp_surf.patch(0).coefs().cols(),1);
    u.setZero();

    gsInfo <<"Normal basis vector at gamma\n";
    cnv_fun<real_t> cnvFun(mp_surf,mp_curve,gamma);
    gsDebug<<cnvFun.eval(u).transpose()<<"\n should match \n";
    gsDebug<<exprEvaluator.eval( gismo::expr::cnv(G_surf, G_curve), gamma ).transpose() << "\n";
    gsInfo <<"First variation of normal basis vector at gamma\n";
    gsDebug<<cnvFun.deriv(u)<<"\n should match \n";
    gsDebug<<exprEvaluator.eval( gismo::expr::cnv_var1(u_surf,G_surf,G_curve), gamma ) << "\n";

    gsInfo <<"Tangent basis vector at gamma\n";
    ctv_fun<real_t> ctvFun(mp_surf,mp_curve,gamma);
    gsDebug<<ctvFun.eval(u).transpose()<<"\n should match \n";
    gsDebug<<exprEvaluator.eval( gismo::expr::ctv(G_surf, G_curve), gamma ).transpose() << "\n";
    gsInfo <<"First variation of tangent basis vector at gamma\n";
    gsDebug<<ctvFun.deriv(u)<<"\n should match \n";
    gsDebug<<exprEvaluator.eval( gismo::expr::ctv_var1(u_surf,G_surf,G_curve), gamma ) << "\n";

    gsInfo <<"Binormal basis vector at gamma\n";
    cbv_fun<real_t> cbvFun(mp_surf,mp_curve,gamma);
    gsDebug<<cbvFun.eval(u).transpose()<<"\n should match \n";
    gsDebug<<exprEvaluator.eval( gismo::expr::cbv(G_surf, G_curve), gamma ).transpose() << "\n";
    gsInfo <<"First variation of binormal basis vector at gamma\n";
    gsDebug<<cbvFun.deriv(u)<<"\n should match \n";
    gsDebug<<exprEvaluator.eval( gismo::expr::cbv_var1(u_surf,G_surf,G_curve), gamma) << "\n";

    gsInfo <<"Derivative of tangent vector along curve\n";
    ctv_along_fun<real_t> ctv_alongFun(mp_surf,mp_curve);
    gsDebug<<ctv_alongFun.deriv(gamma).transpose()<<"\n should match \n";
    gsDebug<<exprEvaluator.eval( gismo::expr::ctv_vara(G_surf, G_curve), gamma ).transpose() << "\n";

    gsInfo <<"Derivative of normal vector along curve\n";
    cnv_along_fun<real_t> cnv_alongFun(mp_surf,mp_curve);
    gsDebug<<cnv_alongFun.deriv(gamma).transpose()<<"\n should match \n";
    gsDebug<<exprEvaluator.eval( gismo::expr::cnv_vara(G_surf, G_curve), gamma ).transpose() << "\n";

    gsInfo <<"Derivative of binormal vector along curve\n";
    cbv_along_fun<real_t> cbv_alongFun(mp_surf,mp_curve);
    gsDebug<<cbv_alongFun.deriv(gamma).transpose()<<"\n should match \n";
    gsDebug<<exprEvaluator.eval( gismo::expr::cbv_vara(G_surf, G_curve), gamma ).transpose() << "\n";

    gsInfo <<"First variation of derivative of tangent basis vector along curve\n";
    ctv_along_var_fun<real_t> ctvalongvarFun(mp_surf,mp_curve,gamma);
    gsDebug<<ctvalongvarFun.deriv(u)<<"\n should match \n";
    gsDebug<<exprEvaluator.eval( gismo::expr::ctv_vara_var1(u_surf,G_surf,G_curve), gamma) << "\n";

    gsInfo <<"First variation of derivative of normal basis vector along curve\n";
    cnv_along_var_fun<real_t> cnvalongvarFun(mp_surf,mp_curve,gamma);
    gsDebug<<cnvalongvarFun.deriv(u)<<"\n should match \n";
    gsDebug<<exprEvaluator.eval( gismo::expr::cnv_vara_var1(u_surf,G_surf,defG_surf,G_curve), gamma) << "\n";

    gsInfo <<"First variation of derivative of binormal basis vector along curve\n";
    cbv_along_var_fun<real_t> cbvalongvarFun(mp_surf,mp_curve,gamma);
    gsDebug<<cbvalongvarFun.deriv(u)<<"\n should match \n";
    gsDebug<<exprEvaluator.eval( gismo::expr::cbv_vara_var1(u_surf,G_surf,defG_surf,G_curve), gamma) << "\n";

    //gsDebug<<exprEvaluator.eval( gismo::expr::cnv_vara_var2dot(u_surf,u_surf,G_surf,G_surf,G_curve,gismo::expr::ctv(G_surf, G_curve)), gamma) << "\n";
    
    /*
        The following functions are provided:
            eps         axial strain component.
            eps_der     first variation of eps. 
            eps_der2    second variation of eps.
            k21         flexural curvature.  
            k21_der     first variationo of k21.  
            k31         flexural curvature.
            k31_der     first variation of k31.
            k23         torsional curvature.
            k23_der     first variation of k23.
            k32         torsional curvature.
            k32_der     first variation of k32.
    **/

    real_t E_modulus_b = 2.10E5; // [MPa]
    real_t PoissonRatio = 0.3;
    real_t thickness_b = 100; // [mm]
    real_t height_b = 300; // [mm]
    auto G_modulus_b = 0.5 * E_modulus_b / (1+PoissonRatio);
    auto area_b = height_b * thickness_b;
    //auto rho = (0.00785*9.81)*area_b; // material density per unit beam length [N/mm]
    auto inertiamin_b = (height_b * pow(thickness_b,3))/12; // minimum principal inertia moment of beam cross-section
    auto inertiamax_b = (thickness_b * pow(height_b,3))/12; // maximum principal inertia moment of beam cross-section
    auto inertiap_b = inertiamin_b + inertiamax_b; // polar inertia moment of beam cross-section

    auto eps = 0.5 * (ctv(defG_surf, G_curve).adj().tr()*ctv(defG_surf, G_curve) - ctv(G_surf, G_curve).adj().tr()*ctv(G_surf, G_curve));
    auto eps_der = ctv_var1(u_surf,defG_surf,G_curve) * ctv(defG_surf, G_curve);
    auto eps_der2 = ctv_var1(u_surf,defG_surf,G_curve)*ctv_var1(u_surf,defG_surf,G_curve).tr();
    
    auto k21 = cnv_vara(defG_surf, G_curve).normalized().adj().tr()*ctv(defG_surf, G_curve) - cnv_vara(G_surf, G_curve).normalized().adj().tr()*ctv(G_surf, G_curve);
    auto k21_der = cnv_vara_var1_normalized(u_surf,G_surf,defG_surf,G_curve)*ctv(defG_surf, G_curve) + ctv_var1(u_surf,defG_surf,G_curve)*cnv_vara(defG_surf, G_curve).normalized();
    auto k21_der2 = cnv_vara_var2dot(u_surf,u_surf,G_surf,defG_surf,G_curve,ctv(defG_surf, G_curve)) +
                    2*var1_dot_var1(u_surf,u_surf,G_curve,cnv_vara_var1_normalized(u_surf,G_surf,defG_surf,G_curve),ctv_var1(u_surf,defG_surf,G_curve));      
    
    auto k31 = cbv_vara(defG_surf, G_curve).normalized().adj().tr()*ctv(defG_surf, G_curve) - cbv_vara(G_surf, G_curve).normalized().adj().tr()*ctv(G_surf, G_curve);
    auto k31_der = cbv_vara_var1(u_surf,G_surf,defG_surf,G_curve)*ctv(defG_surf, G_curve) + ctv_var1(u_surf,defG_surf,G_curve)*cbv_vara(defG_surf, G_curve).normalized();
    auto k31_der2 = cbv_vara_var2dot(u_surf,u_surf,G_surf,defG_surf,G_curve,ctv(defG_surf, G_curve)) +
                    2*var1_dot_var1(u_surf,u_surf,G_curve,cbv_vara_var1_normalized(u_surf,G_surf,defG_surf,G_curve),ctv_var1(u_surf,defG_surf,G_curve));
    
    auto k23 = cnv_vara(defG_surf, G_curve).normalized().adj().tr()*cbv(defG_surf, G_curve) - cnv_vara(G_surf, G_curve).normalized().adj().tr()*cbv(G_surf, G_curve);
    auto k23_der = cnv_vara_var1_normalized(u_surf,G_surf,defG_surf,G_curve)*cbv(defG_surf, G_curve) + cbv_var1(u_surf,defG_surf,G_curve)*cnv_vara(defG_surf, G_curve).normalized();
    auto k23_der2 = cnv_vara_var2dot(u_surf,u_surf,G_surf,defG_surf,G_curve,cbv(defG_surf, G_curve)) +
                    2*var1_dot_var1(u_surf,u_surf,G_curve,cnv_vara_var1_normalized(u_surf,G_surf,defG_surf,G_curve),cbv_var1(u_surf,defG_surf,G_curve)) +
                    cbv_var2dot(u_surf,u_surf,defG_surf,G_curve,cnv_vara(defG_surf,G_curve));
    
    auto k32 = cbv_vara(defG_surf, G_curve).normalized().adj().tr()*cnv(defG_surf, G_curve) - cbv_vara(G_surf, G_curve).normalized().adj().tr()*cnv(G_surf, G_curve);
    auto k32_der = cbv_vara_var1(u_surf,G_surf,defG_surf,G_curve)*cnv(defG_surf, G_curve) + cnv_var1(u_surf,defG_surf,G_curve)*cbv_vara(defG_surf, G_curve).normalized();
    auto k32_der2 = cbv_vara_var2dot(u_surf,u_surf,G_surf,defG_surf,G_curve,cnv(defG_surf,G_curve)) +
                    2*var1_dot_var1(u_surf,u_surf,G_curve,cbv_vara_var1_normalized(u_surf,G_surf,defG_surf,G_curve),cbv_var1(u_surf,defG_surf,G_curve)) +
                    cnv_var2dot(u_surf,u_surf,G_curve,cbv_vara(defG_surf,G_curve));

    // integrands in beam stiffness matrix
    auto KintBendingAxial_b = E_modulus_b/pow(ctv(G_surf, G_curve).norm(),4) * 
                          (area_b * eps_der * eps_der.tr() + inertiamax_b * k21_der * k21_der.tr() + inertiamin_b * k31_der * k31_der.tr() +
                           area_b * eps.val() * eps_der2 + inertiamax_b * k21.val() * k21_der2 + inertiamin_b * k31.val() * k31_der2);
    auto KintTorsion_b = G_modulus_b*inertiap_b/(4*pow(ctv(G_surf, G_curve).norm(),2)) *
                         (-k32_der * k32_der.tr() + k23_der * k23_der.tr() + k23.val() * k23_der2 + (-k32.val()) * k32_der2);

    // integrands in beam residual force vector
    auto RintBendingAxial_b = E_modulus_b/pow(ctv(G_surf, G_curve).norm(),4) * 
                         (area_b * eps.val() * eps_der + inertiamax_b * k21.val() * k21_der + inertiamin_b * k31.val() * k31_der);
    auto RintTorsion_b = G_modulus_b*inertiap_b/(4*pow(ctv(G_surf, G_curve).norm(),2)) *
                         (-k32.val() * k32_der + k23.val() * k23_der);

    gsDebugVar(exprEvaluator.eval(k21_der2,gamma));
    // gsDebugVar(KintTorsion_b.rows());
    // gsDebugVar(KintTorsion_b.cols());

    gsWriteParaview(mp_surf,  "surf",  1000, true,  false); // multi-patch, file name, number of points, mesh, control net
    gsWriteParaview(mp_curve, "curve", 1000, false, false); // multi-patch, file name, number of points, mesh, control net

    return EXIT_SUCCESS;
}

template <class T>
T curveCoordinate(const gsGeometry<T> & geometry, const T & value, const short_t & dir)
{
    GISMO_ASSERT(geometry.parDim() == 1, "The geometry must be a curve with parameter dimension 1");
    GISMO_ASSERT(dir < geometry.targetDim(), "The direction must be less than the target dimension of the geometry");
    gsMatrix<T> u(1,1);
    u(0,0) = value;
    return geometry.eval(u)(dir,0);
    // return geometry.deriv(u)(dir,0);
};