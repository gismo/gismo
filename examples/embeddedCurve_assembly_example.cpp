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
        GISMO_ASSERT(_C.source().targetDim() == _S.source().domainDim(), "The curve and surface must have the same target dimension");
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
        dtheta = _C.data().values[1].reshapeCol(k,1,_C.data().dim.second).transpose(); // VERIFY
        // We take the gradient of the surface at the parametric coordinates theta
        sJac  = Spatch.deriv(theta);
        // HUGO: Check conventions
        // dG = [[dx/dxi,  dy/dxi,  dz/dxi],
        //       [dx/deta, dy/deta, dz/deta]]
        sJac.transposeInPlace();
        sJac.resize(_S.source().domainDim(),_S.source().targetDim());

        if (_S.source().targetDim() == 2)
        {
            res.col(0).head(2) = sJac*dtheta;
            res(2,0)           = 0;
        }
        else if (_S.source().targetDim()==3)
            res = sJac*dtheta;
        else
            GISMO_ERROR("The target dimension of the surface must be 2 or 3, but is "<<_S.source().targetDim());

        return res;
    }

    index_t rows() const { return _C.source().targetDim(); }

    index_t cols() const { return 1; }

    void parse(gsExprHelper<Scalar> & evList) const
    {
        // evList.add(_S);
        // _S.data().flags |= NEED_DERIV;

        evList.add(_C);
        _C.data().flags |= NEED_VALUE | NEED_DERIV;
    }

    const gsFeSpace<Scalar> & rowVar() const {return gsNullExpr<Scalar>::get();}
    const gsFeSpace<Scalar> & colVar() const {return gsNullExpr<Scalar>::get();}

    void print(std::ostream &os) const { os << "ctv("; _S.print(os); os<<" , "; _C.print(os); os <<")"; }
};

template<class E>
class curve_tangent_var1_expr : public _expr<curve_tangent_var1_expr<E> >
{
    typedef typename E::Scalar Scalar;

    typename E::Nested_t _u;
    typename gsGeometryMap<Scalar>::Nested_t _S;
    typename gsGeometryMap<Scalar>::Nested_t _C;

public:
    enum{ Space = E::Space, ScalarValued= 0, ColBlocks= 0};

    curve_tangent_var1_expr(const E & u,
                            const gsGeometryMap<Scalar> & S,
                            const gsGeometryMap<Scalar> & C)
    :
    _u(u),
    _S(S),
    _C(C),
    Spatch(_S.source().piece(0)) //////// <<<<<<<<---------- REMOVE THIS! MAKE PATCH-INDEPENDENT PRE-COMPUTATION
    {
        GISMO_ASSERT(_S.source().domainDim() == 2, "The geometry must be a surface with domain dimension 2");
        GISMO_ASSERT(_C.source().domainDim() == 1, "The geometry must be a curve with domain dimension 1");
        GISMO_ASSERT(_C.source().targetDim() == _S.source().domainDim(), "The curve and surface must have the same target dimension");
    }

    mutable gsMatrix<Scalar> theta, dtheta, bGrads, sJac, res;
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
        // The parametric coordinates of the surface are the evaluation of the curve at point k
        theta = _C.data().values[0].col(k);
        dtheta = _C.data().values[1].reshapeCol(k,1,_C.data().dim.second).transpose(); // VERIFY

        const gsMultiBasis<Scalar> & mb = static_cast<const gsMultiBasis<Scalar>&>(_u.source());
        GISMO_ASSERT( dynamic_cast<const gsMultiBasis<Scalar>*>(&_u.source()), "error");

        const index_t A = mb.basis(0).active(theta).rows(); // _u.data().actives.rows()

        res.resize(A*_u.dim(),cols()); // cols==3
        res.setZero();

        // Take the gradient of the basis functions
        bGrads = mb.basis(0).deriv(theta);

        // bGrads = _u.data().values[1].col(k);
        for (index_t d = 0; d!= _u.dim(); ++d) // for all basis function components
        {
            /// NOT USED YET
            sJac  = Spatch.deriv(theta);
            sJac.transposeInPlace();
            sJac.resize(_S.source().domainDim(),_S.source().targetDim());

            const short_t s = d*A;
            for (index_t j = 0; j!= A; ++j) // for all actives
            {
                // a_{1,r} = vecFun(d, bGrads.at(2*j  ))
                // a_{2,r} = vecFun(d, bGrads.at(2*j+1))
                // dtheta_1= dtheta(0,0);
                // dtheta_2= dtheta(1,0);
                res.row(s+j) = vecFun(d, bGrads.at(2*j  ))* dtheta(0,0)
                             + vecFun(d, bGrads.at(2*j+1))* dtheta(1,0);

                             /// FORGOTTEN: THE DIVISION BY THE LENGTH OF THE VECTOR
            }
        }
        return res;
    }

    index_t rows() const { return 1; }

    index_t cols() const { return 3; }

    index_t cardinality_impl() const { return _u.cardinality_impl(); }

    void parse(gsExprHelper<Scalar> & evList) const
    {
        // CANNOT BE PRECOMPUTED SINCE IT HAS TO BE COMPUTED ON C
        // evList.add(_u);
        // _u.data().flags |= NEED_ACTIVE | NEED_VALUE | NEED_GRAD; // need actives for cardinality

        // evList.add(_S);
        // _S.data().flags |= NEED_DERIV;

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
        GISMO_ASSERT(_C.source().targetDim() == _S.source().domainDim(), "The curve and surface must have the same target dimension, but are "<<_C.source().targetDim()<<" and "<<_S.source().domainDim());
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
        if (_S.source().targetDim() == 2)
        {
            res.resize(3,1);
            res << 0, 0, 1;
            return res;
        }
        else if (_S.source().targetDim()==3)
        {
            // The parametric coordinates of the surface are the evaluation of the curve at point k
            theta = _C.data().values[0].col(k);
            // We take the gradient of the surface at the parametric coordinates theta
            sJac  = Spatch.deriv(theta);
            // HUGO: Check conventions
            // dG = [[dx/dxi,  dy/dxi,  dz/dxi],
            //       [dx/deta, dy/deta, dz/deta]]
            sJac.transposeInPlace();
            sJac.resize(_S.source().domainDim(),_S.source().targetDim());
            gsDebugVar(sJac);
            res = sJac.col3d(0).cross(sJac.col3d(1));
            return res;
        }
        else
            GISMO_ERROR("The target dimension of the surface must be 2 or 3, but is "<<_S.source().targetDim());
    }

    index_t rows() const { return 3; }

    index_t cols() const { return 1; }

    void parse(gsExprHelper<Scalar> & evList) const
    {
        // evList.add(_S);
        // _S.data().flags |= NEED_DERIV;

        evList.add(_C);
        _C.data().flags |= NEED_VALUE;
    }

    const gsFeSpace<Scalar> & rowVar() const {return gsNullExpr<Scalar>::get();}
    const gsFeSpace<Scalar> & colVar() const {return gsNullExpr<Scalar>::get();}

    void print(std::ostream &os) const { os << "cnv("; _S.print(os); os<<" , "; _C.print(os); os <<")"; }
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

    mutable gsMatrix<Scalar> theta, bGrads, sJac, cnvMat, res, var;
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
        res.resize(_u.cardinality(), cols());

        if (_S.source().targetDim() == 2)
        {
            res.setZero();
            return res;
        }
        else if (_S.source().targetDim()==3)
        {
            const index_t A = _u.cardinality()/_u.dim();
            bGrads = _u.data().values[1].col(k);
            theta  = _C.data().values[0].col(k);
            sJac  = Spatch.deriv(theta);
            sJac.transposeInPlace();
            sJac.resize(_S.source().domainDim(),_S.source().targetDim());
            curve_normal_expr<Scalar> cnv(_S,_C);
            cnvMat = cnv.eval(k);
            const Scalar measure = cnvMat.norm();
            cnvMat.normalize();

            for (index_t d = 0; d!= cols(); ++d) // for all basis function components
            {
                const short_t s = d*A;
                for (index_t j = 0; j!= A; ++j) // for all active basis functions
                {
                    //first variation of non-unit normal vector (divided by measure)
                    var =  (vecFun(d, bGrads.at(2*j  )).cross(sJac.col3d(1))
                        - vecFun(d, bGrads.at(2*j+1)).cross(sJac.col3d(0)))/measure;

                    //first variation of unit normal vector
                    res.row(s+j) = (var - ((cnvMat.col3d(0)*var.transpose()) * cnvMat.col3d(0))).transpose();
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
        evList.add(_u);
        _u.data().flags |= NEED_ACTIVE | NEED_GRAD;

        evList.add(_S);
        _S.data().flags |= NEED_DERIV;

        evList.add(_C);
        _C.data().flags |= NEED_VALUE;
    }

    const gsFeSpace<Scalar> & rowVar() const { return _u.rowVar(); }
    const gsFeSpace<Scalar> & colVar() const {return gsNullExpr<Scalar>::get();}
    index_t cardinality_impl() const { return _u.cardinality_impl(); }

    void print(std::ostream &os) const { os << "var1("; _u.print(os); os <<")"; }

};

template<class T>
class curve_binormal_expr : public _expr<curve_binormal_expr<T> >
{
    typedef T Scalar;

    typename gsGeometryMap<Scalar>::Nested_t _S;
    typename gsGeometryMap<Scalar>::Nested_t _C;

public:
    enum{ Space = 0, ScalarValued= 0, ColBlocks= 0};

    curve_binormal_expr(const gsGeometryMap<Scalar> & S,
                      const gsGeometryMap<Scalar> & C)
    :
    _S(S),
    _C(C),
    Spatch(S.source().piece(0)) //////// <<<<<<<<---------- REMOVE THIS! MAKE PATCH-INDEPENDENT PRE-COMPUTATION
    {
        GISMO_ASSERT(_S.source().domainDim() == 2, "The geometry must be a surface with domain dimension 2, but is "<<_S.source().domainDim());
        GISMO_ASSERT(_C.source().domainDim() == 1, "The geometry must be a curve with domain dimension 1, but is "<<_C.source().domainDim());
        GISMO_ASSERT(_C.source().targetDim() == _S.source().domainDim(), "The curve and surface must have the same target dimension, but are "<<_C.source().targetDim()<<" and "<<_S.source().domainDim());
    }

    mutable gsMatrix<Scalar> cnvMat, ctvMat, res;
    const gsFunctionSet<Scalar> & Spatch;
#   define Eigen gsEigen
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
#   undef Eigen


    // Make the selection between 2D and 3D based on enable_if (and get the expression templated over d)
    // to avoid if-statement
    const gsMatrix<Scalar> & eval(const index_t k) const
    {
        curve_normal_expr<T>  cnv(_S,_C);
        cnvMat = cnv.eval(k);
        curve_tangent_expr<T> ctv(_S,_C);
        ctvMat = ctv.eval(k);
        res = cnvMat.col3d(0).cross(ctvMat.col3d(0));
        return res;
    }

    index_t rows() const { return 3; }

    index_t cols() const { return 1; }

    void parse(gsExprHelper<Scalar> & evList) const
    {
        // evList.add(_S);
        // _S.data().flags |= NEED_DERIV;

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

    curve_binormal_var1_expr(const E & u, const gsGeometryMap<Scalar> & S, const gsGeometryMap<Scalar> & C): _u(u), _S(S), _C(C),
    Spatch(_S.source().piece(0)) //////// <<<<<<<<---------- REMOVE THIS! MAKE PATCH-INDEPENDENT PRE-COMPUTATION
    {}

    mutable gsMatrix<Scalar> ctvMat, cnvMat, cbvMat, ctv1Mat, cnv1Mat, var, res;
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

    // Make the selection between 2D and 3D based on enable_if (and get the expression templated over d)
    // to avoid if-statement
    const gsMatrix<Scalar> & eval(const index_t k) const
    {
        var.resize(_u.cardinality(), cols());
        res.resize(_u.cardinality(), cols());

        curve_tangent_expr<Scalar> ctv(_S,_C);
        ctvMat = ctv.eval(k);
        gsDebug<<"1\n";
        curve_normal_expr<Scalar> cnv(_S,_C);
        cnvMat = cnv.eval(k);
        gsDebug<<"2\n";
        curve_tangent_var1_expr<E> ctv_var1(_u,_S,_C);
        ctv1Mat = ctv_var1.eval(k).transpose();
        gsDebug<<"3\n";
        curve_normal_var1_expr<E> cnv_var1(_u,_S,_C);
        cnv1Mat = cnv_var1.eval(k).transpose();
        gsDebug<<"4\n";
        curve_binormal_expr<Scalar> cbv(_S,_C);
        cbvMat = cbv.eval(k);
        gsDebug<<"5\n";
        const Scalar measure = cbvMat.norm();
        cbvMat.normalize();

        // Loop over rows of (ctv1Mat or cnv1Mat) == cardinality of space
        for (index_t i = 0; i!= _u.cardinality(); ++i)
        {
            //first variation of non-unit binormal vector (divided by measure)
            var.row(i) = (cnv1Mat.col3d(i).cross(ctvMat.col3d(0)) +
                         cnvMat.col3d(0).cross(ctv1Mat.col3d(i)))/measure;

            //first variation of unit binormal vector
            res.row(i) = (var.row(i).transpose() - (cbvMat.col3d(0)*var.row(i)) *cbvMat.col3d(0)).transpose();
        }
        return res;
    }

    index_t rows() const { return 1; }

    index_t cols() const { return 3; }

    void parse(gsExprHelper<Scalar> & evList) const
    {
        // For pre-computation of expressions for curve_normal_expr, curve_tangent_expr, curve_tangent_var1_expr, curve_normal_var1_expr, curve_binormal_expr
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


template<class T> EIGEN_STRONG_INLINE
curve_tangent_expr<T> ctv(const gsGeometryMap<T> & S,
                          const gsGeometryMap<T> & C) { return curve_tangent_expr<T>(S,C); }

template<class E> EIGEN_STRONG_INLINE
curve_tangent_var1_expr<E> ctv_var1(const E &u, const gsGeometryMap<typename E::Scalar> & S,
                          const gsGeometryMap<typename E::Scalar> & C) { return curve_tangent_var1_expr<E>(u,S,C); }

template<class T> EIGEN_STRONG_INLINE
curve_normal_expr<T>  cnv(const gsGeometryMap<T> & S,
                          const gsGeometryMap<T> & C) { return curve_normal_expr<T>(S,C); }

template<class E> EIGEN_STRONG_INLINE
curve_normal_var1_expr<E> cnv_var1(const E &u, const gsGeometryMap<typename E::Scalar> & S,
                          const gsGeometryMap<typename E::Scalar> & C) { return curve_normal_var1_expr<E>(u,S,C); }

template<class T> EIGEN_STRONG_INLINE
curve_binormal_expr<T>  cbv(const gsGeometryMap<T> & S,
                            const gsGeometryMap<T> & C) { return curve_binormal_expr<T>(S,C); }

template<class E> EIGEN_STRONG_INLINE
curve_binormal_var1_expr<E> cbv_var1(const E &u, const gsGeometryMap<typename E::Scalar> & S,
                          const gsGeometryMap<typename E::Scalar> & C) { return curve_binormal_var1_expr<E>(u,S,C); }

}
}

using namespace gismo;

template <class T>
T curveCoordinate(const gsGeometry<T> & geometry, const T & value, const short_t & dir);

// template <class T>
// T findRoot();

template <class T = real_t>
class cnv_fun : public gsFunction<T>
{
public:
    // Constructor
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

    short_t targetDim() const override { return 3; } // = gismo::expr::cnv().rows()

    void eval_into(const gsMatrix<T>  & u, gsMatrix<T> &result) const override
    {
        // We need to evaluate the cnv expression at point gamma.
        result.resize(targetDim(),u.cols());

        for (index_t k = 0; k!= u.cols(); ++k)
        {
            // Deform the surface with u.col(k)
            gsMultiPatch<T> surf_def = m_surf;
            surf_def.patch(0).coefs() += u.reshapeCol(k,m_surf.patch(0).coefs().rows(),m_surf.patch(0).coefs().cols());

            // Construct expression evaluator
            gsExprEvaluator<T> exprEvaluator;
            auto G_curve = exprEvaluator.getMap(m_curve);
            auto G_surf  = exprEvaluator.getMap(surf_def);

            // Evaluate expression
            result.col(k) = exprEvaluator.eval( gismo::expr::cnv(G_surf, G_curve), m_gamma);
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

    // Make a BSpline curve
    gsKnotVector<> kv_c(0, 1, 2, 3); //start,end,interior knots, start/end multiplicity
    gsMatrix<> coefs_c(5, 2); //ξ,η;..
    coefs_c << 0, 0,
             0.2618, 0.053,
             0.5, 0.5,
             0.738, 0.9465,
             1,1;

    gsBSpline<> curve( kv_c, give(coefs_c));

    // Print the Bspline curve
    gsInfo << "I am a " << curve << "\n";

    // 2D->2D
    gsMultiPatch<> mp_surf;
    mp_surf.addPatch(gsNurbsCreator<>::BSplineSquare());
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
    // auto u_curve = exprAssembler.getSpace(basis_curve); // not needed...
    auto u_surf  = exprAssembler.getSpace(basis_surf,mp_surf.geoDim());
    gsDebugVar(mp_surf.geoDim());
    // Declare expression evaluator from assembler
    gsExprEvaluator<> exprEvaluator(exprAssembler);

    // auto u_surf = exprAssembler.getSpace(basis_surf,mp_surf.geoDim(),G_curve);
    // auto G_surf = exprAssembler.getMap(mp_surf,G_curve);

    gsVector<> gamma(1);
    gamma << 0.5;
    gsInfo<<"Values at point "<<gamma<<"\n";
    gsInfo << "theta           = "<<exprEvaluator.eval( G_curve, gamma ).transpose() << "\n";
    gsInfo << "normal vector   = "<<exprEvaluator.eval( gismo::expr::cnv(G_surf, G_curve), gamma ).transpose() << "\n";
    gsInfo << "tangent vector  = "<<exprEvaluator.eval( gismo::expr::ctv(G_surf, G_curve), gamma ).transpose() << "\n";
    gsInfo << "binormal vector = "<<exprEvaluator.eval( gismo::expr::cbv(G_surf, G_curve), gamma ).transpose() << "\n";
    gsInfo << "variation of tv = "<<exprEvaluator.eval( gismo::expr::ctv_var1(u_surf,G_surf, G_curve), gamma ) << "\n";
    // gsInfo << "variation of nv = "<<exprEvaluator.eval( gismo::expr::cnv_var1(u_surf,G_surf, G_curve), gamma ) << "\n";
    // BREAKS:
    // gsInfo << "variation of bv = "<<exprEvaluator.eval( gismo::expr::cbv_var1(u_surf,G_surf, G_curve), gamma ) << "\n";


    /*
        What we have is `cnv(G_surf,G_curve)`. G_curve is fixed, so it does not change with the displacement of the surface.
        Hence, cnv(G_surf | G_curve). We want to find the variation of cnv w.r.t. G_surf or the displacement vector.

        How do we validate this?
        - We write a "function" that computes cnv given a 'perturbation vector'.
        - If we take the finite difference of that "function", we should get the same result as the variation of cnv.

        How to do this in G+Smo?
        - We write a gsFunction-derived class that takes the role of "function" above
        - Inside this gsFunction, we implement eval_into() that computes the normal vector given a perturbation vector
        - Given eval_into, gsFunction implements finite difference automatically in deriv_into()
     */
    GISMO_ENSURE(mp_surf.nPatches()==1, "The surface must be a single patch");
    gsMatrix<> u(mp_surf.patch(0).coefs().rows()*mp_surf.patch(0).coefs().cols(),1);
    u.setZero();

    // We test our derivative in gamma
    cnv_fun<real_t> cnvFun(mp_surf,mp_curve,gamma);
    gsDebug<<cnvFun.eval(u); // should give the same result as "normal vector   = "<<exprEvaluator.eval( gismo::expr::cnv(G_surf, G_curve), gamma ).transpose() << "\n";
    gsDebug<<cnvFun.deriv(u); // should give the same result as "variation of nv = "<<exprEvaluator.eval( gismo::expr::cnv_var1(u_surf,G_surf, G_curve), gamma ) << "\n";

    //gsDebugVar(exprEvaluator.eval(jac(u),point));

    // gsMatrix<> normal = exprEvaluator.eval( sn(G_surf) );
    // gsMatrix<> tangent = exprEvaluator.eval( tv(G_curve) );
    // gsMatrix<> binormal = normal.cross(tangent);

    // exprEvaluator.integral()


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