/** @file gsEmbeddingUtils.h

    @brief Utilities for gsEmbedding. Mainly expressions

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s):
        C. Chianese     (2025 UniNa)
        H.M. Verhelst   (2025 TU Delft)
*/

#pragma once

//! [Include namespace]
#include <gsCore/gsLinearAlgebra.h>
#include <gsCore/gsBasis.h>
#include <gsCore/gsFuncData.h>
#include <gsCore/gsDofMapper.h>

#include <gsPde/gsBoundaryConditions.h>

#include <gsUtils/gsPointGrid.h>

#include <gsAssembler/gsAssemblerOptions.h>
#include <gsAssembler/gsExpressions.h>


#  define MatExprType  auto

namespace gismo{
namespace expr{

/**
 * @brief      expressions of local basis vectors on the embedded curve
 */

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
    :_S(S),_C(C),Spatch(_S.source().piece(0)) //////// <<<<<<<<---------- REMOVE THIS! MAKE PATCH-INDEPENDENT PRE-COMPUTATION
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

        theta = _C.data().values[0].col(k);
        dtheta = _C.data().values[1].reshapeCol(k,1,_C.data().dim.second).transpose();
        sJac  = Spatch.deriv(theta);
        sJac.resize(_S.source().domainDim(),_S.source().targetDim());
        sJac.transposeInPlace();

        if (_S.source().targetDim() == 2)
        {
            res.col(0).head(2) = sJac*dtheta;
            res(2,0)           = 0;
        }
        else if (_S.source().targetDim()==3)
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
    :_S(S),_C(C),Spatch(S.source().piece(0)) //////// <<<<<<<<---------- REMOVE THIS! MAKE PATCH-INDEPENDENT PRE-COMPUTATION
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

    // TODO: Make this work for a 2D planar surface AND for a 3D surface:
    // - 2D: normal vector is always (0,0,1)
    // - 3D: normal vector is the cross product of the tangent vectors
    // Make the selection between 2D and 3D based on enable_if (and get the expression templated over d)
    // to avoid if-statement
    const gsMatrix<Scalar> & eval(const index_t k) const
    {
        if (_S.source().targetDim() == 2)
        {
            res.resize(3,1);
            res << 0, 0, 1;
            return res;
        }
        else if (_S.source().targetDim()==3)
        {
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


    // TODO: Make the selection between 2D and 3D based on enable_if (and get the expression templated over d) to avoid if-statement
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

        cnv(_S,_C).parse(evList);
        ctv(_S,_C).parse(evList);
    }

    const gsFeSpace<Scalar> & rowVar() const {return gsNullExpr<Scalar>::get();}
    const gsFeSpace<Scalar> & colVar() const {return gsNullExpr<Scalar>::get();}

    void print(std::ostream &os) const { os << "cnv("; _S.print(os); os<<" , "; _C.print(os); os <<")"; }
};

/**
 * @brief First variation of local basis vectors wrt DOFs.
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

        res.resize(A*_u.dim(),cols());
        res.setZero();

        bGrads = mb.basis(0).deriv(theta); // bGrads = _u.data().values[1].col(k);

        for (index_t d = 0; d!= _u.dim(); ++d)
        {
            const short_t s = d*A;
            for (index_t j = 0; j!= A; ++j)
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

        if (_S.source().targetDim() == 2)
        {
            res.setZero();
            return res;
        }
        else if (_S.source().targetDim()==3)
        {
            bGrads = mb.basis(0).deriv(theta);
            sJac  = Spatch.deriv(theta);
            sJac.resize(_S.source().domainDim(),_S.source().targetDim());
            sJac.transposeInPlace();
            curve_normal_expr<Scalar> cnv(_S,_C);
            cnvec = cnv.eval(k);
            const Scalar measure = cnvec.norm();
            cnvec.normalize();

            for (index_t d = 0; d!= cols(); ++d)
            {
                const short_t s = d*A;
                for (index_t j = 0; j!= A; ++j)
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

    // TODO: Make the selection between 2D and 3D based on enable_if (and get the expression templated over d) to avoid if-statement
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
 * @brief Derivative of local basis vectors along the embedded curve.
 */

template<class T>
class curve_tangent_var_along_expr : public _expr<curve_tangent_var_along_expr<T> >
{
public:
    typedef T Scalar;

    typename gsGeometryMap<Scalar>::Nested_t _S;
    typename gsGeometryMap<Scalar>::Nested_t _C;

public:
    enum{ Space = 0, ScalarValued= 0, ColBlocks= 0};

    curve_tangent_var_along_expr(const gsGeometryMap<Scalar> & S,
                                 const gsGeometryMap<Scalar> & C): _S(S), _C(C),
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

    void print(std::ostream &os) const { os << "ctv_var_along("; _S.print(os); os<<" , "; _C.print(os); os <<")"; }
};

template<class T>
class curve_normal_var_along_expr : public _expr<curve_normal_var_along_expr<T> >
{
public:
    typedef T Scalar;

    typename gsGeometryMap<Scalar>::Nested_t _S;
    typename gsGeometryMap<Scalar>::Nested_t _C;

public:
    enum{ Space = 0, ScalarValued= 0, ColBlocks= 0};

    curve_normal_var_along_expr(const gsGeometryMap<Scalar> & S,
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

    void print(std::ostream &os) const { os << "cnv_var_along("; _S.print(os); os<<" , "; _C.print(os); os <<")"; }
};

template<class T>
class curve_binormal_var_along_expr : public _expr<curve_binormal_var_along_expr<T> >
{
public:
    typedef T Scalar;

    typename gsGeometryMap<Scalar>::Nested_t _S;
    typename gsGeometryMap<Scalar>::Nested_t _C;

public:
    enum{ Space = 0, ScalarValued= 0, ColBlocks= 0};

    curve_binormal_var_along_expr(const gsGeometryMap<Scalar> & S,
                                  const gsGeometryMap<Scalar> & C): _S(S), _C(C) {}

    mutable gsMatrix<Scalar> ctvec, ctvec_vara, cnvec, cnvec_vara, res;
#   define Eigen gsEigen
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
#   undef Eigen

    const gsMatrix<Scalar> & eval(const index_t k) const
    {
        curve_tangent_var_along_expr<T> ctv_vara(_S,_C);
        ctvec_vara = ctv_vara.eval(k);
        curve_normal_var_along_expr<T> cnv_vara(_S,_C);
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

        ctv(_S,_C).parse(evList);
        ctv_vara(_S,_C).parse(evList);
        cnv(_S,_C).parse(evList);
        cnv_vara(_S,_C).parse(evList);
    }

    const gsFeSpace<Scalar> & rowVar() const {return gsNullExpr<Scalar>::get();}
    const gsFeSpace<Scalar> & colVar() const {return gsNullExpr<Scalar>::get();}

    void print(std::ostream &os) const { os << "cbv_var_along("; _S.print(os); os<<" , "; _C.print(os); os <<")"; }
};

/**
 * @brief First variation of the derivative of local basis vectors along the embedded curve wrt DOFs.
 */

template<class E>
class curve_tangent_var_along_var1_expr : public _expr<curve_tangent_var_along_var1_expr<E> >
{
public:
    typedef typename E::Scalar Scalar;

private:
    typename E::Nested_t _u;
    typename gsGeometryMap<Scalar>::Nested_t _S;
    typename gsGeometryMap<Scalar>::Nested_t _C;

public:
    enum{ Space = E::Space, ScalarValued= 0, ColBlocks= 0};

    curve_tangent_var_along_var1_expr(const E & u, const gsGeometryMap<Scalar> & S, const gsGeometryMap<Scalar> & C): _u(u), _S(S), _C(C) {}

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
        const index_t A = mb.basis(0).active(theta).rows(); // _u.data().actives.rows()
        bGrads = mb.basis(0).deriv(theta);
        bHess =  mb.basis(0).deriv2(theta);
        vec.head(2) = dtheta.pow(2);
        vec(2) = 2*dtheta(0,0)*dtheta(1,0);

        res.resize(A*_u.dim(), cols());

        for (index_t d = 0; d!= cols(); ++d)
        {
            const short_t s = d*A;
            for (index_t j = 0; j!= A; ++j)
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

    void print(std::ostream &os) const { os << "curve_tangent_var_along_var1("; _u.print(os); os <<")"; }

};

template<class E>
class curve_normal_var_along_var1_expr : public _expr<curve_normal_var_along_var1_expr<E> >
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

    curve_normal_var_along_var1_expr(const E & u, const gsGeometryMap<Scalar> & S,
                                     const gsGeometryMap<Scalar> & defS,
                                     const gsGeometryMap<Scalar> & C): _u(u), _S(S), _defS(defS), _C(C),
                                     Spatch(S.source().piece(0)),defSpatch(defS.source().piece(0)) //////// <<<<<<<<---------- REMOVE THIS! MAKE PATCH-INDEPENDENT PRE-COMPUTATION
                                     {}

    mutable gsMatrix<Scalar> theta, dtheta, bGrads, bHess, localbHess, uJac, uHess, sol, res;
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

        curve_normal_var_along_expr<Scalar> cnv_vara(_S,_C);
        cnvec_vara = cnv_vara.eval(k);
        const Scalar measure = cnvec_vara.norm();
        cnvec_vara.normalize();

        res.resize(A*_u.dim(), cols());
        var.resize(A*_u.dim(), cols());

        for (index_t d = 0; d!= cols(); ++d)
        {
            const short_t s = d*A;
            for (index_t j = 0; j!= A; ++j)
            {
                secDerToHessian(bHess.col(0).segment(3*j,3),2,localbHess);
                localbHess.resize(2,2);
                dot = localbHess * dtheta;

                //not normalized result
                res.row(s+j) =  vecFun(d,dot(0)).cross(uJac.col3d(1)) +
                                dotA.cross(vecFun(d,bGrads.at(2*j+1))) +
                                vecFun(d,bGrads.at(2*j  )).cross(dotB) +
                                uJac.col3d(0).cross(vecFun(d,dot(1)));
           
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

    void print(std::ostream &os) const { os << "curve_normal_var_along_var1("; _u.print(os); os <<")"; }

};

template<class E>
class curve_binormal_var_along_var1_expr : public _expr<curve_binormal_var_along_var1_expr<E> >
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

    curve_binormal_var_along_var1_expr(const E & u, const gsGeometryMap<Scalar> & S,
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
        const index_t A = mb.basis(0).active(theta).rows(); // _u.data().actives.rows()
        
        curve_tangent_expr<Scalar> ctv(_S,_C);
        ctvec = ctv.eval(k);
        curve_tangent_var_along_expr<Scalar> ctv_vara(_S,_C);
        ctvec_vara = ctv_vara.eval(k);
        curve_tangent_var1_expr<E> ctv_var1(_u,_S,_C);
        ctvec_var1 = ctv_var1.eval(k).transpose();
        curve_tangent_var_along_var1_expr<E> ctv_vara_var1(_u,_S,_C);
        ctvec_vara_var1 = ctv_vara_var1.eval(k).transpose();

        curve_normal_expr<Scalar> cnv(_S,_C);
        cnvec = cnv.eval(k);
        curve_normal_var_along_expr<Scalar> cnv_vara(_S,_C);
        cnvec_vara = cnv_vara.eval(k);
        curve_normal_var1_expr<E> cnv_var1(_u,_S,_C);
        cnvec_var1 = cnv_var1.eval(k).transpose();
        curve_normal_var_along_var1_expr<E> cnv_vara_var1(_u,_S,_defS,_C);
        cnvec_vara_var1 = cnv_vara_var1.eval(k).transpose();
        curve_binormal_var_along_expr<Scalar> cbv_vara(_S,_C);
        cbvec_vara = cbv_vara.eval(k);
        const Scalar measure = cbvec_vara.norm();
        cbvec_vara.normalize();

        var.resize(A*_u.dim(),cols());
        res.resize(A*_u.dim(),cols());

        for (index_t i = 0; i!= A*_u.dim(); ++i)
        {
            //not normalized result
            // res.row(i) =  ctvec_vara_var1.col3d(i).cross(cnvec.col3d(0)) +
            //               ctvec_vara.col3d(0).cross(cnvec_var1.col3d(i)) +
            //               ctvec_var1.col3d(i).cross(cnvec_vara.col3d(0)) +
            //               ctvec.col3d(0).cross(cnvec_vara_var1.col3d(i));

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
        // WHY IS THIS NEEDED?
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

    void print(std::ostream &os) const { os << "curve_binormal_var_along_var1("; _u.print(os); os <<")"; }
};

template<class T> EIGEN_STRONG_INLINE
curve_tangent_expr<T> ctv(const gsGeometryMap<T> & S,
                          const gsGeometryMap<T> & C) { return curve_tangent_expr<T>(S,C); }

template<class E> EIGEN_STRONG_INLINE
curve_tangent_var1_expr<E> ctv_var1(const E &u, const gsGeometryMap<typename E::Scalar> & S,
                          const gsGeometryMap<typename E::Scalar> & C) { return curve_tangent_var1_expr<E>(u,S,C); }

template<class T> EIGEN_STRONG_INLINE
curve_tangent_var_along_expr<T> ctv_vara(const gsGeometryMap<T> & S,
                                const gsGeometryMap<T> & C) { return curve_tangent_var_along_expr<T>(S,C); }

template<class E> EIGEN_STRONG_INLINE
curve_tangent_var_along_var1_expr<E> ctv_vara_var1(const E &u, const gsGeometryMap<typename E::Scalar> & S,
                                     const gsGeometryMap<typename E::Scalar> & C) { return curve_tangent_var_along_var1_expr<E>(u,S,C); }

template<class T> EIGEN_STRONG_INLINE
curve_normal_expr<T>  cnv(const gsGeometryMap<T> & S,
                          const gsGeometryMap<T> & C) { return curve_normal_expr<T>(S,C); }

template<class E> EIGEN_STRONG_INLINE
curve_normal_var1_expr<E> cnv_var1(const E &u, const gsGeometryMap<typename E::Scalar> & S,
                          const gsGeometryMap<typename E::Scalar> & C) { return curve_normal_var1_expr<E>(u,S,C); }

template<class T> EIGEN_STRONG_INLINE
curve_normal_var_along_expr<T>  cnv_vara(const gsGeometryMap<T> & S,
                                const gsGeometryMap<T> & C) { return curve_normal_var_along_expr<T>(S,C); }

template<class E> EIGEN_STRONG_INLINE
curve_normal_var_along_var1_expr<E> cnv_vara_var1(const E &u, const gsGeometryMap<typename E::Scalar> & S,
                                    const gsGeometryMap<typename E::Scalar> & defS,
                                    const gsGeometryMap<typename E::Scalar> & C)
                                    { return curve_normal_var_along_var1_expr<E>(u,S,defS,C); }

template<class T> EIGEN_STRONG_INLINE
curve_binormal_expr<T>  cbv(const gsGeometryMap<T> & S,
                            const gsGeometryMap<T> & C) { return curve_binormal_expr<T>(S,C); }

template<class E> EIGEN_STRONG_INLINE
curve_binormal_var1_expr<E> cbv_var1(const E &u, const gsGeometryMap<typename E::Scalar> & S,
                          const gsGeometryMap<typename E::Scalar> & C) { return curve_binormal_var1_expr<E>(u,S,C); }

template<class T> EIGEN_STRONG_INLINE
curve_binormal_var_along_expr<T> cbv_vara(const gsGeometryMap<T> & S,
                                 const gsGeometryMap<T> & C) { return curve_binormal_var_along_expr<T>(S,C); }

template<class E> EIGEN_STRONG_INLINE
curve_binormal_var_along_var1_expr<E> cbv_vara_var1(const E &u, const gsGeometryMap<typename E::Scalar> & S,
                                      const gsGeometryMap<typename E::Scalar> & defS,
                                      const gsGeometryMap<typename E::Scalar> & C) { return curve_binormal_var_along_var1_expr<E>(u,S,defS,C); }
}
}