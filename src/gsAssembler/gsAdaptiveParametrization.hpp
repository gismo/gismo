/** @file gsSquareDomain.hpp

    @brief Implementation of the gsSquareDomain class.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris
*/

#pragma once


#include <gsCore/gsComposedFunction.h>
#include <gsNurbs/gsSquareDomain.h>
#include <gsAssembler/gsExpressions.h>
#include <gsAssembler/gsExprHelper.h>
#include <gsAssembler/gsExprEvaluator.h>
#include <gsNurbs/gsTensorNurbsBasis.h>

namespace gismo
{

namespace expr
{

// template<class E0, class E1, class E2>
// class ternary_expr : public _expr<ternary_expr<E0, E1, E2> >
// {
//   typename E0::Nested_t _u;
//   typename E1::Nested_t _v;
//   typename E2::Nested_t _w;
//  public:
//   typedef typename E1::Scalar Scalar;

//   explicit ternary_expr(_expr<E0> const &u,
//                         _expr<E1> const &v,
//                         _expr<E2> const &w)
//       :
//       _u(u),
//       _v(v),
//       _w(w) {
//     GISMO_ASSERT(E0::ScalarValued, "Condition must be scalar valued");
//     GISMO_ASSERT((int) E1::ScalarValued == (int) E2::ScalarValued,
//                  "Both v and w must be scalar valued (or not).");
//     GISMO_ASSERT((int) E1::ColBlocks == (int) E2::ColBlocks,
//                  "Both v and w must be colblocks (or not).");
//     GISMO_ASSERT((int) E1::Space == (int) E2::Space,
//                  "Both v and w must be space (or not), but E1::Space = "
//                      << E1::Space << " and E2::Space = " << E2::Space);
//     GISMO_ASSERT(_v.rows() == _w.rows(),
//                  "Rows of v and w differ. _v.rows() = " << _v.rows()
//                                                         << ", _w.rows() = "
//                                                         << _w.rows());
//     GISMO_ASSERT(_v.cols() == _w.cols(),
//                  "Columns of v and w differ. _v.cols() = " << _v.cols()
//                                                            << ", _w.cols() = "
//                                                            << _w.cols());
//     GISMO_ASSERT(_v.rowVar() == _w.rowVar(), "rowVar of v and w differ.");
//     GISMO_ASSERT(_v.colVar() == _w.colVar(), "colVar of v and w differ.");
//   }
//  public:
//   enum {
//     ScalarValued = E1::ScalarValued,
//     ColBlocks = E1::ColBlocks,
//     Space = E1::Space
//   }; // == E2::Space

// //  const Scalar eval(const index_t k) const { return (_u.eval(k) > 0 ? _v.eval
// //  (k) : _w.eval(k)); }

//   const Temporary_t eval(const index_t k) const
//   {
//     return (_u.eval(k) > 0 ? _v.eval(k) : _w.eval(k));
//   }

//   // { res = eval_impl(_u,_v,_w,k); return  res;}

//   index_t rows() const { return _v.rows(); }
//   index_t cols() const { return _v.cols(); }
//   void parse(gsExprHelper<Scalar> &evList) const {
//     _u.parse(evList);
//     _v.parse(evList);
//     _w.parse(evList);
//   }

//   const gsFeSpace<Scalar> &rowVar() const { return _v.rowVar(); }
//   const gsFeSpace<Scalar> &colVar() const { return _v.colVar(); }

// };


// /// Ternary ternary_expr
// template<class E0, class E1, class E2>
// EIGEN_STRONG_INLINE
// ternary_expr<E0, E1, E2> ternary(const E0 &u,
//                                  const E1 &v,
//                                  const E2 &w)
// {
//   return ternary_expr<E0, E1, E2>(u, v, w);
// }

/*
template<class T>
class dmeas_expr : public _expr<dmeas_expr<T> >
{
    typename gsGeometryMap<T>::Nested_t _G;

public:
    enum {Space = 0, ScalarValued = 0, ColBlocks = 0};

    typedef T Scalar;

    mutable gsMatrix<Scalar> grad, hess; // temporary storage for evaluation
    mutable gsMatrix<Scalar> cross;
    mutable gsMatrix<Scalar> res;

    dmeas_expr(const gsGeometryMap<T> & G) : _G(G) { }


    const gsMatrix<Scalar> & eval(const index_t k) const
    {
        res.resize(1,2); // [dMeasure/dxi, dMeasure/deta]
        grad = _G.data().values[1].reshapeCol(k, _G.data().dim.first, _G.data().dim.second).transpose();

        // hess = [dG_x/dxidxi, dG_x/dxideta, dG_x/detadxi, dG_x/detadeta;
        //         dG_y/dxidxi, dG_y/dxideta, dG_y/detadxi, dG_y/detadeta;
        //         dG_z/dxidxi, dG_z/dxideta, dG_z/detadxi, dG_z/detadeta;

        // hess = [dG_x/dxidxi, dG_x/detadeta, dG_x/detadxi
        //         dG_y/dxidxi, dG_y/detadeta, dG_y/detadxi
        //         dG_z/dxidxi, dG_z/detadeta, dG_z/detadxi
        hess.resize(3,3);
        hess.row(0) = _G.data().values[2].col(k).segment(0,3).transpose();
        hess.row(1) = _G.data().values[2].col(k).segment(3,3).transpose();
        hess.row(2) = _G.data().values[2].col(k).segment(6,3).transpose();
        cross = grad.col(0).cross(grad.col(1));

        res(0,0) = ( cross / cross.norm() ).dot( hess.col(0).cross(grad.col(1)) + grad.col(0).cross(hess.col(2)) ); // cross/||cross|| . ( dG/dxixi   x dG/deta + dG/dxi x dG/dxideta  )
        res(0,1) = ( cross / cross.norm() ).dot( hess.col(2).cross(grad.col(0)) + grad.col(1).cross(hess.col(1)) ); // cross/||cross|| . ( dG/detadxi x dG/deta + dG/dxi x dG/detadeta )

        return res;
    }

    index_t cols() const { return _G.source().domainDim(); }

    void parse(gsExprHelper<Scalar> & evList) const
    {
        evList.add(_u);
        _u.data().flags |= NEED_GRAD | NEED_DERIV2;
    }

    const gsFeSpace<Scalar> & rowVar() const {return gsNullExpr<Scalar>::get();}

    const gsFeSpace<Scalar> & colVar() const {return gsNullExpr<Scalar>::get();}

    void print(std::ostream &os) const { os << "????("; _G.print(os); os <<")"; }
private:

    template<class U> static inline
    typename util::enable_if<util::is_same<U,gsComposition<Scalar> >::value,
                             const gsMatrix<Scalar> &>::type
    eval_impl(const U & u, const index_t k)
    {
        return u.eval(k);
    }
};

/// TODO
template<class E> EIGEN_STRONG_INLINE
dmeasure_expr<E> dmeas(const E & u) { return dmeasure_expr<E>(u); }
*/

template<class T>
class dEdSigma_expr : public _expr<dEdSigma_expr<T> >
{
    typename gsGeometryMap<T>::Nested_t _G;

public:
    enum {Space = 0, ScalarValued = 0, ColBlocks = 0};

    typedef T Scalar;

    mutable gsMatrix<Scalar,3,2> grad;
    mutable gsMatrix<Scalar,3,3> hess; // temporary storage for evaluation
    mutable gsMatrix<Scalar,3,2> dJdxi, dJdeta;
    mutable gsVector<Scalar,3> cross;
    mutable Scalar E;
    mutable Scalar meas;
    mutable Scalar dmeasdxi, dmeasdeta;
    mutable gsMatrix<Scalar> res;

    dEdSigma_expr(const gsGeometryMap<T> & G)
    :
    _G(G)
    {
        GISMO_ASSERT( _G.source().domainDim()==2 && _G.source().targetDim()==3,"dEdSigma_expr: The geometry must be a surface");
    }


    const gsMatrix<Scalar> & eval(const index_t k) const
    {
         res.resize(2,1); // [dMeasure/dxi, dMeasure/deta]

        grad = _G.data().values[1].reshapeCol(k, _G.data().dim.first, _G.data().dim.second).transpose();

        // dG/dxi x dG/deta
        cross = grad.col(0).cross(grad.col(1));
        meas  = cross.norm();

        // grad = [dG_x/dxi, dG_x/deta
        //         dG_y/dxi, dG_y/deta
        //         dG_z/dxi, dG_z/deta];

        // hess = [dG_x/dxidxi, dG_x/detadeta, dG_x/detadxi
        //         dG_y/dxidxi, dG_y/detadeta, dG_y/detadxi
        //         dG_z/dxidxi, dG_z/detadeta, dG_z/detadxi

        // dJdxi = [dG_x/dxidxi, dG_x/detadxi
        //          dG_y/dxidxi, dG_y/detadxi
        //          dG_z/dxidxi, dG_z/detadxi];

        // dJdeta= [dG_x/dxideta, dG_x/detadeta
        //          dG_y/dxideta, dG_y/detadeta
        //          dG_z/dxideta, dG_z/detadeta];
        hess.row(0) = _G.data().values[2].col(k).segment(0,3).transpose();
        hess.row(1) = _G.data().values[2].col(k).segment(3,3).transpose();
        hess.row(2) = _G.data().values[2].col(k).segment(6,3).transpose();

        dJdxi.col(0) = hess.col(0);
        dJdxi.col(1) = hess.col(2);

        dJdeta.col(0) = hess.col(2);
        dJdeta.col(1) = hess.col(1);

        dmeasdxi  = ( cross / meas ).dot( hess.col(0).cross(grad.col(1)) + grad.col(0).cross(hess.col(2)) ); // cross/||cross|| . ( dG/dxixi   x dG/deta + dG/dxi x dG/dxideta  )
        dmeasdeta = ( cross / meas ).dot( hess.col(2).cross(grad.col(1)) + grad.col(0).cross(hess.col(1)) ); // cross/||cross|| . ( dG/detadxi x dG/deta + dG/dxi x dG/detadeta )

        E = (grad.transpose()*grad).trace() / meas;

        res(0,0) = ( 2 * (grad.transpose()*dJdxi).trace() )  / meas - E*dmeasdxi /meas;
        res(1,0) = ( 2 * (grad.transpose()*dJdeta).trace() ) / meas - E*dmeasdeta/meas;
        return res;

    }

    index_t rows() const { return 2; }

    index_t cols() const { return 1; }

    void parse(gsExprHelper<Scalar> & evList) const
    {
        evList.add(_G);
        _G.data().flags |= NEED_GRAD | NEED_DERIV2;
    }

    const gsFeSpace<Scalar> & rowVar() const {return gsNullExpr<Scalar>::get();}

    const gsFeSpace<Scalar> & colVar() const {return gsNullExpr<Scalar>::get();}

    void print(std::ostream &os) const { os << "dEdσ("; _G.print(os); os <<")"; }
};

/// TODO
template<class T> EIGEN_STRONG_INLINE
dEdSigma_expr<T> dEdSigma(const gsGeometryMap<T> & G) { return dEdSigma_expr<T>(G); }


template<enum MonitorMode MODE, class E>
class monitor_expr : public _expr<monitor_expr<MODE,E> >
{
public:
    typedef typename E::Scalar Scalar;

private:
    typename E::Nested_t _u;
    typename gsGeometryMap<Scalar>::Nested_t _G;
    const index_t DIM;

    mutable gsMatrix<Scalar> res;
    mutable gsMatrix<Scalar> grad, jac, jacInv;
    mutable gsMatrix<Scalar> ones;
    mutable Scalar m_theta;

public:
    enum{ Space = E::Space, ScalarValued= 0, ColBlocks= 0};

    monitor_expr(const E & u, const gsGeometryMap<Scalar> & G, const Scalar theta)
    :
    _u(u),
    _G(G),
    DIM(_u.source().domainDim()),
    m_theta(theta)
    {
    }

#   define Eigen gsEigen
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
#   undef Eigen

    const gsMatrix<Scalar> & eval(const index_t k) const
    {return eval_impl<MODE>(_u,k); }

    index_t rows() const { return DIM; }

    index_t cols() const { return DIM; }

    void parse(gsExprHelper<Scalar> & evList) const
    {
        parse_impl<MODE>(evList);
    }

    const gsFeSpace<Scalar> & rowVar() const { return _u.rowVar(); }
    const gsFeSpace<Scalar> & colVar() const {return gsNullExpr<Scalar>::get();}
    index_t cardinality_impl() const { return _u.cardinality_impl(); }

    void print(std::ostream &os) const { os << "M("; _u.print(os); os <<")"; }

private:

    template<enum MonitorMode _MODE> inline
    typename util::enable_if<_MODE == ValueBased, void>::type
    parse_impl(gsExprHelper<Scalar> & evList) const
    {
        _G.parse(evList);

        evList.add(_u);
        _u.data().flags |= NEED_VALUE; //for value-based
    }

    template<enum MonitorMode _MODE> inline
    typename util::enable_if<_MODE == GradientBased, void>::type
    parse_impl(gsExprHelper<Scalar> & evList) const
    {
        jacInv_expr<Scalar>(_G).parse(evList);
        _G.parse(evList);

        evList.add(_u);
        _u.data().flags |= NEED_JACOBIAN;
    }

    template<enum MonitorMode _MODE> inline
    typename util::enable_if< _MODE != ValueBased && _MODE != GradientBased , void>::type
    parse_impl(gsExprHelper<Scalar> & evList) const
    {
        GISMO_NO_IMPLEMENTATION;
    }

    template<enum MonitorMode _MODE, class U> inline
    typename util::enable_if< util::is_same<U,gsFeSpace<Scalar> >::value, const gsMatrix<Scalar> & >::type
    eval_impl(const U & u, const index_t k)  const
    {
        GISMO_NO_IMPLEMENTATION;
    }

    template<enum MonitorMode _MODE, class U> inline
    typename util::enable_if< !util::is_same<U,gsFeSpace<Scalar> >::value && _MODE == ValueBased, const gsMatrix<Scalar> & >::type
    eval_impl(const U & u, const index_t k)  const
    {
        ones.resize(DIM,DIM);
        ones.setIdentity();

        // Value-based:
        Scalar uval = _u.data().values[0].col(k).value();
        res = 1.0 / ( math::sqrt(1.0 +m_theta*math::pow(uval,2) ) ) * ones;
        return res;
    }

    template<enum MonitorMode _MODE, class U> inline
    typename util::enable_if< !util::is_same<U,gsFeSpace<Scalar> >::value && _MODE == GradientBased, const gsMatrix<Scalar> & >::type
    eval_impl(const U & u, const index_t k)  const
    {
        ones.resize(DIM,DIM);
        ones.setIdentity();

        // Gradient-based
        jac = _u.data().values[1].col(k).transpose();
        jacInv_expr<Scalar> jacInvExpr = jacInv_expr<Scalar>(_G);
        jacInv = jacInvExpr.eval(k);
        grad = jac * jacInv;
        res = 1.0 / ( math::sqrt(1.0 +m_theta*grad.squaredNorm() ) ) * ones;
        return res;
    }

    template<enum MonitorMode _MODE, class U> inline
    typename util::enable_if< !util::is_same<U,gsFeSpace<Scalar> >::value && _MODE != ValueBased && _MODE != GradientBased, const gsMatrix<Scalar> & >::type
    eval_impl(const U & u, const index_t k)  const
    {
        GISMO_NO_IMPLEMENTATION;
    }

};


template<enum MonitorMode MODE, class E>
EIGEN_STRONG_INLINE
monitor_expr<MODE,E> monitor(const E & u,
                             const gsGeometryMap<typename E::Scalar> & G,
                             const typename E::Scalar theta = 0.1)
                             {
                                return monitor_expr<MODE,E>(u,G,theta);
                             }


} // namespace expr

template<class T, enum MonitorMode MODE>
gsOptMesh<T,MODE>::gsOptMesh()
{}

template<class T, enum MonitorMode MODE>
gsOptMesh<T,MODE>::gsOptMesh(         gsFunction<T> & composition,
                                const gsGeometry<T> & geometry,
                                const gsBasis<T>    * integrationBasis)
:
gsOptMesh(composition,geometry,nullptr,integrationBasis,false)
{}

template<class T, enum MonitorMode MODE>
gsOptMesh<T,MODE>::gsOptMesh(         gsFunction<T> & composition,
                                const gsGeometry<T> & geometry,
                                const gsFunction<T> * fun,
                                const gsBasis<T>    * integrationBasis,
                                const bool            parametric)
:
m_comp(&composition),
m_geom(&geometry),
m_fun(fun),
m_ib(integrationBasis),
m_mb(*m_ib,true),
m_cgeom(*m_comp,geometry),
// THIS DOES NOT WORK FOR PARAMETRIC=FALSE. POINTER IS LOST WHEN ARRIVING IN evalObj
// m_cfun(parametric ? gsComposedFunction<T>(*m_comp,fun) : gsComposedFunction<T>(m_cgeom,fun)),
m_mp(m_cgeom),
m_parametric(parametric)
{
    m_numDesignVars = m_comp->nControls();
    m_curDesign.resize(m_numDesignVars,1);
    m_controls.resize(m_numDesignVars,1);
    m_controls.col(0) = m_comp->getControls();

    m_options.addReal("Smoothing","Smoothing parameter for the monitor function",0.1);
    m_options.addReal("Penalty","Penalty parameter for the monitor function",1e-2);
}

template<class T, enum MonitorMode MODE>
gsFunction<T> & gsOptMesh<T,MODE>::composition()
{
    return *m_comp;
}

template<class T, enum MonitorMode MODE>
gsOptionList & gsOptMesh<T,MODE>::options()
{
    return m_options;
}


// // ORIGINAL
// template<class T, enum MonitorMode MODE>
// T gsOptMesh<T,MODE>::evalObj(const gsAsConstVector<T> &u) const
// {
//     typedef typename gsExprHelper<T>::geometryMap geometryMap; ///< Geometry map type
//     gsExprEvaluator<T> evaluator;

//     evaluator.setIntegrationElements(m_mb); // does not work when in constructor
//     m_comp->setControls(u);

//     // Penalty constant
//     gsConstantFunction<T> pen(m_options.getReal("Penalty"), m_cgeom.domainDim());
//     geometryMap G = evaluator.getMap(m_mp);
//     auto eps = evaluator.getVariable(pen);

//     // auto chi = 0.5 * (jac(G).det() + pow(eps.val() + pow(jac(G).det(), 2.0), 0.5));
//     // auto invJacMat = jac(G).adj()/chi;
//     // auto eta = evaluator.getVariable(m_fun);
//     // return evaluator.integral( (monitor(eta,G).asDiag()*invJacMat).sqNorm()*meas(G));

//     // gsComposedFunction<T> fun(*m_comp,m_fun.function(0));
//     // gsComposedFunction<T> fun(mpG.patch(0),m_fun.function(0));


//     // gsVector<unsigned> np(m_fun.domainDim());
//     // np.setConstant(m_options.getInt("nSamplingPoints"));
//     // gsMatrix<T> grid = gsPointGrid<T>(m_cgeom.support().col(0),m_cgeom.support().col(1),np);
//     if (m_cgeom.domainDim()==m_cgeom.targetDim())
//     {
//         auto detG = jac(G).det();
//         auto chi = 0.5 * (detG + pow(pow(eps.val(),2.0) + pow(detG, 2.0), 0.5));
//         auto invJacMat = jac(G).adj()/chi; // inverse of jacobian matrix with 'determinant' replaced

//         if (m_fun==nullptr)
//         {
//             auto M = 1/detG;
//             return evaluator.integral( (M*invJacMat).sqNorm()*meas(G));
//         }
//         else
//         {
//             // TEMPORARY FIX: SEE COMMENT IN CONSTRUCTOR
//             m_cfun = m_parametric ? gsComposedFunction<T>(*m_comp,*m_fun) : gsComposedFunction<T>(m_cgeom,*m_fun);
//             auto eta = evaluator.getVariable(m_cfun);
//             auto M   = gismo::expr::monitor<MODE>(eta,G,m_options.getReal("Smoothing"));
//             return evaluator.integral( (M*invJacMat).sqNorm()*meas(G));
//         }
//     }
//     else if (m_cgeom.domainDim()<m_cgeom.targetDim())
//     {
//         auto fform = jac(G).tr()*jac(G);
//         // SQUARE ROOT:
//         auto detG = pow(fform.det().val(),0.5).val(); //jacobian determinant for a surface, i.e. the measure
//         // SQUARED:
//         // auto detG = fform.det().val(); //jacobian determinant for a surface, i.e. the measure

//         // OPTION 1
//         // Compute the chi part
//         // auto chiPPart = eps * ((detG.val() - eps.val()).exp());
//         // Ternary operation to compute chi and chip
//         // auto chi = ternary(eps.val() - detG, chiPPart.val(), detG.val());
//         // OPTION 2
//         // auto chi = 0.5 * (detG + pow(pow(eps.val(),2.0) + pow(detG, 2.0), 0.5));
//         // OPTION 3
//         auto chi = detG;
//         // auto chiPPart = -(-2.0 + pow(eps.val(),2) + eps.val())*pow(detG,3)/pow(eps.val(),4) + (- 3.0 + 2.0*pow(eps.val(),2) + 2.0*eps.val())*pow(detG,2)/pow(eps.val(),3) - detG/eps.val() + 1.0/eps.val();
//         // auto chi = ternary(eps.val() - detG, chiPPart.val(), detG.val());
//         // SQUARED
//         //
//         // auto chi = 0.5 * (pow(detG,0.5) + pow(pow(eps.val(),2.0) + detG, 0.5));

//         // SQUARE ROOT:
//         auto invJacMat = fform.sqrt().adj()/chi; // inverse of jacobian matrix with 'determinant' replaced
//         // SQUARED:
//         // auto invJacMat = fform.adj()/chi; // inverse of jacobian matrix with 'determinant' replaced

//         if (m_fun==nullptr)
//         {
//             auto M = 1/detG;
//             return evaluator.integral( M*(invJacMat).sqNorm()*meas(G));
//         }
//         else
//         {
//             // TEMPORARY FIX: SEE COMMENT IN CONSTRUCTOR
//             m_cfun = m_parametric ? gsComposedFunction<T>(*m_comp,*m_fun) : gsComposedFunction<T>(m_cgeom,*m_fun);
//             auto eta = evaluator.getVariable(m_cfun);
//             auto M   = gismo::expr::monitor<MODE>(eta,G,m_options.getReal("Smoothing"));
//             return evaluator.integral( M*(invJacMat).sqNorm()*meas(G));
//         }
//     }
//     else
//         GISMO_ERROR("Domain dimension must be smaller than or equal to the target dimension, but domainDim = "<<m_cgeom.domainDim()<<" and targetDim = "<<m_cgeom.targetDim());
// }


/// @cond
//
// Objective function evaluation.
//
// Map chain:   hat{Omega} --sigma(u;alpha)--> tilde{Omega} --G--> Omega
//
// Jacobians:   J_sigma (dd x dd),  J_g (td x dd),  J_c = J_g * J_sigma (td x dd)
// Metric:      C = J_c^T J_c  (dd x dd)
// Geometry metric: C_g = J_g^T J_g  (dd x dd)
//
// WITHOUT monitor:
//   E(alpha) = int  tr(C^{-1}) / sqrt(det C)  d hat{Omega}
//
// WITH monitor (weight m^2 = 1/(1 + theta * eta^2)):
//   E(alpha) = int  m^2 * tr(C^{-1}) * sqrt(det C)  d hat{Omega}
//
// where eta^2 depends on MonitorMode:
//
// ValueBased:
//   eta = f(xi(alpha))          (monitor value at moving parametric point)
//   m^2 = 1 / (1 + theta * eta^2)
//
// GradientBased:
//   eta^2 = (nabla_xi f)^T C_g^{-1} (nabla_xi f)
//         = || nabla_x f ||^2   (squared physical gradient norm, for td == dd)
//   m^2 = 1 / (1 + theta * eta^2)
//
//   For m_parametric=true:   nabla_xi f = deriv of f w.r.t. xi (direct)
//   For m_parametric=false:  nabla_xi f = J_g^T * nabla_x f  (chain rule)
//
//   The C_g^{-1} metric accounts for the coordinate distortion so that
//   eta^2 measures the physical gradient norm regardless of the
//   parametric representation.
//
/// @endcond
template<class T, enum MonitorMode MODE>
T gsOptMesh<T,MODE>::evalObj(const gsAsConstVector<T> &u) const
{
    const index_t dd = m_comp->domainDim();
    const index_t td = m_geom->targetDim();
    GISMO_ASSERT(dd <= td, "domainDim must be <= targetDim");

    m_comp->setControls(u);

    const T theta = m_options.getReal("Smoothing");
    const bool hasMonitor = (m_fun != nullptr);

    T result = T(0);

    gsExprEvaluator<T> evaluator;
    evaluator.options().setReal("quA",0.0);
    evaluator.options().setInt("quB",1);
    evaluator.setIntegrationElements(m_mb);

    for (unsigned patchInd = 0; patchInd < m_mb.nBases(); ++patchInd)
    {
        gsQuadRule<T> QuRule = gsQuadrature::get(m_mb.basis(patchInd), evaluator.options());
        typename gsBasis<T>::domainIter domIt = m_mb.piece(patchInd).makeDomainIterator();

        gsMatrix<T> uvPoints, xietaPoints;
        gsVector<T> tmpWeights;

        for (; domIt->good(); domIt->next())
        {
            QuRule.mapTo(domIt->lowerCorner(), domIt->upperCorner(),
                         uvPoints, tmpWeights);

            m_comp->eval_into(uvPoints, xietaPoints);

            gsMatrix<T> Jsigma_flat;
            m_comp->deriv_into(uvPoints, Jsigma_flat);

            gsMatrix<T> Jgeom_flat;
            m_geom->deriv_into(xietaPoints, Jgeom_flat);

            gsMatrix<T> monVals, monDerivs_eval;
            if (hasMonitor)
            {
                if (m_parametric)
                {
                    m_fun->eval_into(xietaPoints, monVals);
                    if (MODE == GradientBased)
                        m_fun->deriv_into(xietaPoints, monDerivs_eval);
                }
                else
                {
                    gsMatrix<T> physPoints;
                    m_geom->eval_into(xietaPoints, physPoints);
                    m_fun->eval_into(physPoints, monVals);
                    if (MODE == GradientBased)
                        m_fun->deriv_into(physPoints, monDerivs_eval);
                }
            }

            for (index_t p = 0; p != uvPoints.cols(); ++p)
            {
                // J_sigma (dd x dd), J_g (td x dd), J_c = J_g J_sigma (td x dd)
                gsMatrix<T> Js = Jsigma_flat.col(p).reshaped(dd, dd).transpose();
                gsMatrix<T> Jg = Jgeom_flat.col(p).reshaped(dd, td).transpose();
                gsMatrix<T> Jc = Jg * Js;
                gsMatrix<T> C = Jc.transpose() * Jc;
                gsMatrix<T> Cinv = C.inverse();
                T detG = math::sqrt(C.determinant());

                T integrand;
                if (hasMonitor)
                {
                    T m2;
                    if (MODE == ValueBased)
                    {
                        // m^2 = 1 / (1 + theta * f(xi)^2)
                        T eta = monVals(0, p);
                        m2 = T(1) / (T(1) + theta * eta * eta);
                    }
                    else // GradientBased
                    {
                        // eta^2 = (nabla_xi f)^T C_g^{-1} (nabla_xi f)
                        gsMatrix<T> Cg = Jg.transpose() * Jg;
                        gsMatrix<T> Cg_inv = Cg.inverse();

                        gsMatrix<T> grad_xi_f;
                        if (m_parametric)
                        {
                            grad_xi_f = monDerivs_eval.col(p);
                        }
                        else
                        {
                            gsMatrix<T> grad_x_f = monDerivs_eval.col(p).reshaped(td, 1);
                            grad_xi_f = Jg.transpose() * grad_x_f;
                        }

                        T eta2 = (grad_xi_f.transpose() * Cg_inv * grad_xi_f)(0, 0);
                        m2 = T(1) / (T(1) + theta * eta2);
                    }
                    // E += w * m^2 * tr(C^{-1}) * sqrt(det C)
                    integrand = m2 * Cinv.trace() * detG;
                }
                else
                {
                    // E += w * tr(C^{-1}) / sqrt(det C)
                    integrand = Cinv.trace() / detG;
                }

                result += tmpWeights[p] * integrand;
            }
        }
    }
    return result;
}


/// @cond
//
// Analytical gradient of the objective function w.r.t. the control
// variables alpha of the composition sigma.
//
// ---------------------------------------------------------------------------
// Notation
// ---------------------------------------------------------------------------
//
// alpha_{k,d}  : d-th coordinate of the k-th basis function coefficient
// N_k          : k-th basis function of sigma, evaluated at quadrature point
// nabla_hat N_k: gradient of N_k w.r.t. hat{u} (the reference coordinates)
// xi = sigma(hat{u}; alpha)  : parametric point in tilde{Omega}
// G(xi)        : geometry map  tilde{Omega} -> Omega
// J_sigma      : Jacobian of sigma  (dd x dd)
// J_g          : Jacobian of G      (td x dd)
// J_c = J_g J_sigma : composed Jacobian (td x dd)
// C = J_c^T J_c     : composed metric   (dd x dd)
// C_g = J_g^T J_g   : geometry metric   (dd x dd)
//
// ---------------------------------------------------------------------------
// Kinematic derivatives  (d(...)/d alpha_{k,d})
// ---------------------------------------------------------------------------
//
// d xi_i / d alpha_{k,d}  =  N_k * delta_{id}
//
// d J_sigma / d alpha_{k,d}  =:  dJ_s
//   (dJ_s)_{ij} = delta_{id} * (nabla_hat N_k)_j
//   i.e. only row d is nonzero, equal to nabla_hat N_k transposed
//
// d J_g / d alpha_{k,d}  =:  dJ_g       (td x dd)
//   (dJ_g)_{a,j} = N_k * d^2 G_a / (d xi_j d xi_d)
//   This is the second derivative of the geometry map, contracted with
//   d xi/d alpha = N_k e_d.
//
// d J_c / d alpha_{k,d} = dJ_g * J_s + J_g * dJ_s
// d C   / d alpha_{k,d} = dJ_c^T J_c + J_c^T dJ_c
//
// ---------------------------------------------------------------------------
// Integrand derivative  (no monitor)
// ---------------------------------------------------------------------------
//
// E = tr(C^{-1}) / sqrt(det C)
//
// Using d(tr(C^{-1}))/dalpha = -tr(C^{-1} dC C^{-1})
//   and d(det C)/dalpha = det(C) tr(C^{-1} dC)
//   => d(1/sqrt(det C))/dalpha = -1/(2 sqrt(det C)) tr(C^{-1} dC)
//
// dE = [ -tr(C^{-1} dC C^{-1}) - tr(C^{-1})/2 * tr(C^{-1} dC) ]
//      / sqrt(det C)
//
// ---------------------------------------------------------------------------
// Integrand derivative  (with monitor, both ValueBased & GradientBased)
// ---------------------------------------------------------------------------
//
// E = m^2 * tr(C^{-1}) * sqrt(det C)
//
// dE = dm^2/dalpha * tr(C^{-1}) * sqrt(det C)
//    + m^2 * [ -tr(C^{-1} dC C^{-1}) * sqrt(det C)
//              + tr(C^{-1}) * sqrt(det C)/2 * tr(C^{-1} dC) ]
//
// Note the sign change on the detG term compared to the no-monitor
// case, because the integrand uses sqrt(det C) (not 1/sqrt(det C)).
//
// dm^2/dalpha = (dm^2/d eta^2) * (d eta^2 / d alpha):
//
// ValueBased:
//   eta = f(xi(alpha))
//   m^2 = 1/(1 + theta eta^2)
//   dm^2/deta = -2 theta eta / (1 + theta eta^2)^2
//   deta/dalpha_{k,d} = (nabla_xi f)_d * N_k
//     (chain rule: d f(xi(alpha))/dalpha = nabla_xi f . d xi/dalpha)
//
// GradientBased:
//   eta^2 = (nabla_xi f)^T C_g^{-1} (nabla_xi f)
//   m^2 = 1/(1 + theta eta^2)
//   dm^2/d(eta^2) = -theta / (1 + theta eta^2)^2
//
//   d(eta^2)/dalpha_{k,d} = term1 + term2, where:
//
//   Term 1 (from C_g changing):
//     Using d(M^{-1}) = -M^{-1} dM M^{-1}, let v = C_g^{-1} nabla_xi f:
//     term1 = -v^T dC_g v
//     with dC_g = dJ_g^T J_g + J_g^T dJ_g
//
//   Term 2 (from nabla_xi f changing at the moving point):
//     term2 = 2 v^T d(nabla_xi f)/dalpha
//
//     For m_parametric=true:
//       d(nabla_xi f)/dalpha_{k,d} = N_k * H_f * e_d
//         (H_f is the dd x dd parametric Hessian of f)
//
//     For m_parametric=false:
//       nabla_xi f = J_g^T nabla_x f,  so differentiating:
//       d(nabla_xi f)/dalpha_{k,d}
//         = dJ_g^T * nabla_x f
//         + N_k * J_g^T * H_f * J_g e_d
//       where H_f is the td x td physical-space Hessian
//       and J_g e_d is the d-th column of J_g (= d G/d xi_d).
//
/// @endcond
template<class T, enum MonitorMode MODE>
void gsOptMesh<T,MODE>::gradObj_into ( const gsAsConstVector<T> & u, gsAsVector<T> & result) const
{
    const index_t nc = m_comp->nControls();
    const index_t dd = m_comp->domainDim();
    const index_t td = m_geom->targetDim();
    GISMO_ASSERT(dd <= td, "domainDim must be <= targetDim");

    result.resize(nc);
    result.setZero();
    m_comp->setControls(u);

    const T theta = m_options.getReal("Smoothing");
    const bool hasMonitor = (m_fun != nullptr);

    const gsSquareDomain<T> * sqDomain = static_cast<const gsSquareDomain<T>*>(m_comp);
    const gsBasis<T> & sigmaBasis = sqDomain->domain().basis();
    const gsDofMapper & mapper = sqDomain->mapper();
    const index_t S = dd * (dd + 1) / 2;

    gsExprEvaluator<T> evaluator;
    evaluator.options().setReal("quA",0.0);
    evaluator.options().setInt("quB",1);
    evaluator.setIntegrationElements(m_mb);

    for (unsigned patchInd = 0; patchInd < m_mb.nBases(); ++patchInd)
    {
        gsQuadRule<T> QuRule = gsQuadrature::get(m_mb.basis(patchInd), evaluator.options());
        typename gsBasis<T>::domainIter domIt = m_mb.piece(patchInd).makeDomainIterator();

        gsMatrix<T> uvPoints, xietaPoints;
        gsVector<T> tmpWeights;

        for (; domIt->good(); domIt->next())
        {
            // Domain sequence
            // \hat{\Omega} -> m_comp -> \tilde{\Omega} -> m_geom -> \Omega
            //              --->             m_cgeom    -----------> 

            // Integration points in \hat{\Omega}
            QuRule.mapTo(domIt->lowerCorner(), domIt->upperCorner(),
                         uvPoints, tmpWeights);

            // Integration points in \tilde{\Omega}
            m_comp->eval_into(uvPoints, xietaPoints);

            // dsigma/duv
            gsMatrix<T> Jsigma_flat;
            m_comp->deriv_into(uvPoints, Jsigma_flat);

            // dm_geom/dxieta
            gsMatrix<T> Jgeom_flat;
            m_geom->deriv_into(xietaPoints, Jgeom_flat);

            // d^2m_geom/dxieta^2
            gsMatrix<T> deriv2_geom;
            m_geom->deriv2_into(xietaPoints, deriv2_geom);

            // Data about the basis of sigma
            gsMatrix<index_t> actives;
            sigmaBasis.active_into(uvPoints, actives);

            gsMatrix<T> basisVals;
            sigmaBasis.eval_into(uvPoints, basisVals);

            gsMatrix<T> basisDerivs;
            sigmaBasis.deriv_into(uvPoints, basisDerivs);

            // Monitor function values, gradient, and (for GradientBased) Hessian
            gsMatrix<T> monVals, monDerivs, monDeriv2;
            gsMatrix<T> physPoints_grad;
            if (hasMonitor)
            {
                if (m_parametric)
                {
                    m_fun->eval_into(xietaPoints, monVals);
                    m_fun->deriv_into(xietaPoints, monDerivs);
                    if (MODE == GradientBased)
                        m_fun->deriv2_into(xietaPoints, monDeriv2);
                }
                else
                {
                    m_geom->eval_into(xietaPoints, physPoints_grad);
                    m_fun->eval_into(physPoints_grad, monVals);
                    m_fun->deriv_into(physPoints_grad, monDerivs);
                    if (MODE == GradientBased)
                        m_fun->deriv2_into(physPoints_grad, monDeriv2);
                }
            }

            for (index_t p = 0; p != uvPoints.cols(); ++p)
            {
                // J_sigma (dd x dd), J_g (td x dd), J_c = J_g J_sigma
                gsMatrix<T> Js = Jsigma_flat.col(p).reshaped(dd, dd).transpose();
                gsMatrix<T> Jg = Jgeom_flat.col(p).reshaped(dd, td).transpose();
                gsMatrix<T> Jc = Jg * Js;

                // Composed metric C = J_c^T J_c and related quantities
                gsMatrix<T> C = Jc.transpose() * Jc;
                gsMatrix<T> Cinv = C.inverse();
                T detC = C.determinant();
                T detG = math::sqrt(detC);
                T trCinv = Cinv.trace();

                // Monitor-related per-quadrature-point quantities
                T m2 = T(0), dm2_deta2 = T(0);
                gsMatrix<T> gradMon;   // nabla_xi f  (dd x 1, ValueBased only)
                gsMatrix<T> Cg, Cg_inv;
                gsMatrix<T> grad_xi_f; // nabla_xi f  (dd x 1, GradientBased)
                gsMatrix<T> grad_x_f;  // nabla_x f   (td x 1, physical-space, GradientBased)
                gsMatrix<T> Hess_f;    // Hessian of f (dd x dd or td x td)

                if (hasMonitor)
                {
                    if (MODE == ValueBased)
                    {
                        // m^2 = 1/(1+theta*eta^2), dm^2/deta = -2 theta eta / denom^2
                        T eta = monVals(0, p);
                        T denom = T(1) + theta * eta * eta;
                        m2 = T(1) / denom;
                        T dm2_deta = -T(2) * theta * eta / (denom * denom);

                        // nabla_xi f (parametric gradient of the monitor)
                        if (m_parametric)
                            gradMon = monDerivs.col(p);
                        else
                        {
                            gsMatrix<T> gradMon_phys = monDerivs.col(p).reshaped(td, 1);
                            gradMon = Jg.transpose() * gradMon_phys;
                        }
                        // Store dm^2/deta in dm2_deta2 for use in the per-DOF loop
                        dm2_deta2 = dm2_deta;
                    }
                    else // GradientBased
                    {
                        // C_g = J_g^T J_g  (geometry metric)
                        Cg     = Jg.transpose() * Jg;
                        Cg_inv = Cg.inverse();

                        if (m_parametric)
                        {
                            // nabla_xi f from m_fun's parametric derivatives
                            grad_xi_f = monDerivs.col(p);

                            // Reconstruct dd x dd Hessian from compact storage
                            // deriv2 layout: [d^2f/dxi_0^2, d^2f/dxi_1^2, ..., d^2f/(dxi_0 dxi_1), ...]
                            Hess_f.resize(dd, dd);
                            for (index_t i = 0; i != dd; ++i)
                                for (index_t j = 0; j != dd; ++j)
                                {
                                    index_t lo = math::min(i,j);
                                    index_t hi = math::max(i,j);
                                    index_t hidx = (lo==hi) ? lo : dd + lo*(2*dd-lo-3)/2+hi-1;
                                    Hess_f(i,j) = monDeriv2(hidx, p);
                                }
                        }
                        else
                        {
                            // nabla_xi f = J_g^T nabla_x f  (chain rule)
                            grad_x_f  = monDerivs.col(p).reshaped(td, 1);
                            grad_xi_f = Jg.transpose() * grad_x_f;

                            // Reconstruct td x td physical-space Hessian
                            Hess_f.resize(td, td);
                            for (index_t i = 0; i != td; ++i)
                                for (index_t j = 0; j != td; ++j)
                                {
                                    index_t lo = math::min(i,j);
                                    index_t hi = math::max(i,j);
                                    index_t hidx = (lo==hi) ? lo : td + lo*(2*td-lo-3)/2+hi-1;
                                    Hess_f(i,j) = monDeriv2(hidx, p);
                                }
                        }

                        // eta^2 = nabla_xi_f^T C_g^{-1} nabla_xi_f
                        T eta2 = (grad_xi_f.transpose() * Cg_inv * grad_xi_f)(0,0);
                        T denom = T(1) + theta * eta2;
                        m2          = T(1) / denom;
                        // dm^2/d(eta^2) = -theta / (1 + theta eta^2)^2
                        dm2_deta2   = -theta / (denom * denom);
                    }
                }

                const index_t nActive = actives.rows();
                for (index_t loc = 0; loc != nActive; ++loc)
                {
                    const index_t k = actives(loc, p);
                    T Nk = basisVals(loc, p);

                    gsMatrix<T> gradNk(dd, 1);
                    for (index_t j = 0; j != dd; ++j)
                        gradNk(j) = basisDerivs(loc * dd + j, p);

                    for (index_t d = 0; d != dd; ++d)
                    {
                        if (!mapper.is_free(k, 0, d))
                            continue;
                        index_t ii = mapper.index(k, 0, d);

                        // dJ_s / d alpha_{k,d}:  only row d is nonzero
                        gsMatrix<T> dJs = gsMatrix<T>::Zero(dd, dd);
                        for (index_t j = 0; j != dd; ++j)
                            dJs(d, j) = gradNk(j);

                        // dJ_g / d alpha_{k,d}:
                        //   (dJ_g)_{a,j} = N_k * d^2 G_a / (d xi_j d xi_d)
                        gsMatrix<T> dJg = gsMatrix<T>::Zero(td, dd);
                        for (index_t a = 0; a != td; ++a)
                        {
                            for (index_t j = 0; j != dd; ++j)
                            {
                                index_t lo = math::min(d, j);
                                index_t hi = math::max(d, j);
                                index_t hess_idx = (lo == hi) ? lo : dd + lo * (2 * dd - lo - 3) / 2 + hi - 1;
                                dJg(a, j) = Nk * deriv2_geom(a * S + hess_idx, p);
                            }
                        }

                        // dJ_c = dJ_g J_s + J_g dJ_s,  dC = dJ_c^T J_c + J_c^T dJ_c
                        gsMatrix<T> dJc = dJg * Js + Jg * dJs;
                        gsMatrix<T> dC = dJc.transpose() * Jc + Jc.transpose() * dJc;

                        T trCinvdC     = (Cinv * dC).trace();
                        T trCinvdCCinv = (Cinv * dC * Cinv).trace();

                        T dE;
                        if (hasMonitor)
                        {
                            T dm2_dalpha;

                            if (MODE == ValueBased)
                            {
                                // deta/dalpha = (nabla_xi f)_d * N_k
                                T deta_dalpha = gradMon(d) * Nk;
                                dm2_dalpha = dm2_deta2 * deta_dalpha;
                            }
                            else // GradientBased
                            {
                                // d(eta^2)/dalpha = term1 + term2
                                // (see block comment above for full derivation)

                                // dC_g = dJ_g^T J_g + J_g^T dJ_g
                                gsMatrix<T> dCg = dJg.transpose() * Jg + Jg.transpose() * dJg;

                                // v = C_g^{-1} nabla_xi f
                                gsMatrix<T> v = Cg_inv * grad_xi_f;

                                // Term 1: -v^T dC_g v  (metric contribution)
                                T term1 = -(v.transpose() * dCg * v)(0,0);

                                // Term 2: 2 v^T d(nabla_xi f)/dalpha  (gradient-shift contribution)
                                gsMatrix<T> d_grad_xi_f(dd, 1);
                                if (m_parametric)
                                {
                                    // d(nabla_xi f)/dalpha = N_k * H_f * e_d
                                    d_grad_xi_f = Nk * Hess_f.col(d);
                                }
                                else
                                {
                                    // d(nabla_xi f)/dalpha = dJ_g^T nabla_x f
                                    //                      + N_k J_g^T H_f J_g e_d
                                    gsMatrix<T> Jg_col_d = Jg.col(d);
                                    d_grad_xi_f = dJg.transpose() * grad_x_f
                                                + Nk * Jg.transpose() * (Hess_f * Jg_col_d);
                                }

                                T term2 = T(2) * (v.transpose() * d_grad_xi_f)(0,0);

                                T deta2_dalpha = term1 + term2;
                                dm2_dalpha = dm2_deta2 * deta2_dalpha;
                            }

                            // dE = d(m^2 tr(C^{-1}) sqrt(det C)) / dalpha
                            dE = m2 * (-trCinvdCCinv * detG
                                       + trCinv * detG / T(2) * trCinvdC)
                               + dm2_dalpha * trCinv * detG;
                        }
                        else
                        {
                            // dE = d(tr(C^{-1}) / sqrt(det C)) / dalpha
                            dE = (-trCinvdCCinv
                                  - trCinv / T(2) * trCinvdC) / detG;
                        }

                        result(ii) += tmpWeights[p] * dE;
                    }
                }
            }
        }
    }
}

template<class T, enum MonitorMode MODE>
void gsOptMesh<T,MODE>::gradObj_FD_into( const gsAsConstVector<T> & u, gsAsVector<T> & result) const
{
    const index_t n = u.rows();
    const T h = T(1e-7);
    result.resize(n);

    gsVector<T> uu = u;
    gsAsVector<T> tmp(uu.data(), n);
    gsAsConstVector<T> ctmp(uu.data(), n);

    for (index_t i = 0; i < n; ++i)
    {
        tmp[i] = u[i] + h;
        const T fp = this->evalObj(ctmp);
        tmp[i] = u[i] - h;
        const T fm = this->evalObj(ctmp);
        tmp[i] = u[i];
        result[i] = (fp - fm) / (T(2) * h);
    }
}

/*
TODO:
ADD constructor in .h file
delegate construction to a separate (DIM) templated function
 */
template <class T, enum MonitorMode MODE>
gsAdaptiveParametrization<T,MODE>::gsAdaptiveParametrization(         gsFunction<T>  & composition,
                                                                const gsGeometry<T>  & geometry,
                                                                const gsBasis<T>     & integrationBasis,
                                                                      gsOptimizer<T> & optimizer,
                                                                const bool             parametric)
:
gsAdaptiveParametrization(composition,geometry,nullptr,integrationBasis,optimizer,parametric)
{
}

template <class T, enum MonitorMode MODE>
gsAdaptiveParametrization<T,MODE>::gsAdaptiveParametrization(         gsFunction<T>  & composition,
                                                                const gsGeometry<T>  & geometry,
                                                                const gsFunction<T>  & function,
                                                                const gsBasis<T>     & integrationBasis,
                                                                      gsOptimizer<T> & optimizer,
                                                                const bool             parametric)
:
gsAdaptiveParametrization(composition,geometry,&function,integrationBasis,optimizer,parametric)
{
}

template <class T, enum MonitorMode MODE>
gsAdaptiveParametrization<T,MODE>::gsAdaptiveParametrization(     gsFunction<T>  & composition,
                                                            const gsGeometry<T>  & geometry,
                                                                  gsOptimizer<T> & optimizer,
                                                            const bool             parametric)
:
gsAdaptiveParametrization(composition,geometry,nullptr,geometry.basis(),optimizer,parametric)
{
}

template <class T, enum MonitorMode MODE>
gsAdaptiveParametrization<T,MODE>::gsAdaptiveParametrization(         gsFunction<T>  & composition,
                                                                const gsGeometry<T>  & geometry,
                                                                const gsFunction<T>  & function,
                                                                      gsOptimizer<T> & optimizer,
                                                                const bool             parametric)
:
gsAdaptiveParametrization(composition,geometry,function,geometry.basis(),optimizer,parametric)
{
}

template <class T, enum MonitorMode MODE>
gsAdaptiveParametrization<T,MODE>::gsAdaptiveParametrization(         gsFunction<T>  & composition,
                                                                const gsGeometry<T>  & geometry,
                                                                const gsFunction<T>  * function,
                                                                const gsBasis<T>     & integrationBasis,
                                                                      gsOptimizer<T> & optimizer,
                                                                const bool             parametric)
:
m_comp(composition),
m_geom(geometry),
m_fun(function),
m_optimizer(optimizer),
m_integrationBasis(integrationBasis.clone())//,
{
    m_optProblem = gsOptMesh<T,MODE>(m_comp,m_geom,m_fun,m_integrationBasis.get(),parametric);
    this->defaultOptions();
}

template <class T, enum MonitorMode MODE>
gsAdaptiveParametrization<T,MODE>::gsAdaptiveParametrization(         gsSquareDomain<T> & composition,
                                                                const gsGeometry<T>     & geometry,
                                                                const gsBasis<T>        & integrationBasis,
                                                                      gsOptimizer<T>    & optimizer,
                                                                const bool                parametric)
:
gsAdaptiveParametrization(composition,geometry,nullptr,integrationBasis,optimizer,parametric)
{
}

template <class T, enum MonitorMode MODE>
gsAdaptiveParametrization<T,MODE>::gsAdaptiveParametrization(         gsSquareDomain<T> & composition,
                                                                const gsGeometry<T>     & geometry,
                                                                const gsFunction<T>     & function,
                                                                const gsBasis<T>        & integrationBasis,
                                                                      gsOptimizer<T>    & optimizer,
                                                                const bool                parametric)
:
gsAdaptiveParametrization(composition,geometry,&function,integrationBasis,optimizer,parametric)
{
}

template <class T, enum MonitorMode MODE>
gsAdaptiveParametrization<T,MODE>::gsAdaptiveParametrization(         gsSquareDomain<T> & composition,
                                                                const gsGeometry<T>     & geometry,
                                                                      gsOptimizer<T>    & optimizer,
                                                                const bool                parametric)
:
gsAdaptiveParametrization(composition,geometry,nullptr,geometry.basis(),optimizer,parametric)
{
}

template <class T, enum MonitorMode MODE>
gsAdaptiveParametrization<T,MODE>::gsAdaptiveParametrization(         gsSquareDomain<T> & composition,
                                                                const gsGeometry<T>     & geometry,
                                                                const gsFunction<T>     & function,
                                                                      gsOptimizer<T>    & optimizer,
                                                                const bool                parametric)
:
gsAdaptiveParametrization(composition,geometry,&function,geometry.basis(),optimizer,parametric)
{
}

template <class T, enum MonitorMode MODE>
gsAdaptiveParametrization<T,MODE>::gsAdaptiveParametrization(         gsSquareDomain<T> & composition,
                                                                const gsGeometry<T>     & geometry,
                                                                const gsFunction<T>     * function,
                                                                const gsBasis<T>        & integrationBasis,
                                                                      gsOptimizer<T>    & optimizer,
                                                                const bool                parametric)
:
m_comp(composition),
m_geom(geometry),
m_fun(function),
m_optimizer(optimizer)
{
    GISMO_ASSERT((dynamic_cast<const gsTensorBSplineBasis<2,T> *>(&composition.domain().basis())),"The composition must be a tensor B-spline basis");
    const gsTensorBSplineBasis<2,T> & comp_tbasis = static_cast<const gsTensorBSplineBasis<2,T> &>(composition.domain().basis());
    if (const gsTensorBSplineBasis<2,T> * tbasis = dynamic_cast<const gsTensorBSplineBasis<2,T> *>(&integrationBasis))
    {
        gsTensorBSplineBasis<2,T> ibasis = makeIntegrationBasis(*tbasis,comp_tbasis);
        m_integrationBasis = memory::make_unique(new gsTensorBSplineBasis<2,T>(ibasis));
    }
    else if (const gsTensorNurbsBasis<2,T> * nbasis = dynamic_cast<const gsTensorNurbsBasis<2,T> *>(&integrationBasis))
    {
        gsTensorNurbsBasis<2,T> ibasis = makeIntegrationBasis(nbasis->source(),comp_tbasis);
        m_integrationBasis = memory::make_unique(new gsTensorNurbsBasis<2,T>(ibasis));
    }
    else
        GISMO_ERROR("The integration basis must be either a tensor B-spline or a tensor NURBS basis");
    m_optProblem = gsOptMesh<T,MODE>(m_comp,m_geom,m_fun,m_integrationBasis.get(),parametric);
    this->defaultOptions();
}

template <class T, enum MonitorMode MODE>
template <short_t _d>
gsTensorBSplineBasis<_d,T> gsAdaptiveParametrization<T,MODE>::makeIntegrationBasis(const gsTensorBSplineBasis<_d,T> & basis1,
                                                                                  const gsTensorBSplineBasis<_d,T> & basis2)
{
    gsTensorBSplineBasis<_d,T> ibasis(basis1);
    // Integration basis: parent basis with knots of composition basis inserted, and the degree is the sum of the two degrees (?)
    index_t targetDegree;
    for (size_t d = 0; d!=_d; d++)
    {
        // 1. Insert interior knots of composition basis
        for (typename gsKnotVector<T>::uiterator it = std::next(basis2.knots(d).ubegin());
                                                    it!= std::prev(basis2.knots(d).uend());
                                                    ++it)
            {
                if (ibasis.knots(d).has(*it))
                    continue;
                ibasis.insertKnot(*it,d);
            }
        // 2. Increase the degree
        targetDegree = ibasis.degree(d) * basis2.degree(d);
        ibasis.degreeIncrease(targetDegree-ibasis.degree(d),d);

    }
    return ibasis;
}

template <class T, enum MonitorMode MODE>
gsOptionList & gsAdaptiveParametrization<T,MODE>::options()
{
    return m_options;
}

template <class T, enum MonitorMode MODE>
void gsAdaptiveParametrization<T,MODE>::defaultOptions()
{
    m_options.addReal("Penalty","Penalization coefficient for Jacobian determinant [default=0.01]",1e-2);
    m_options.addReal("Smoothing","Smoothing parameter in the monitor function [default=0.1]",0.1);
    m_options.addInt("Mode","0: Relocate based on f [default]; 1: Relocate based on grad(f)",0);
}

template <class T, enum MonitorMode MODE>
void gsAdaptiveParametrization<T,MODE>::solve()
{
    // Set the optimization problem
    m_optProblem.options().setReal("Smoothing",m_options.getReal("Smoothing"));
    m_optProblem.options().setReal("Penalty",m_options.getReal("Penalty"));
    m_optimizer.setProblem(&m_optProblem);
    // Solve the optimization problem
    gsVector<T> controls = m_optProblem.composition().getControls();
    m_optimizer.solve(controls);
    controls = m_optimizer.currentDesign();
    m_optProblem.composition().setControls(controls);
    gsInfo<<"Finished with objective value: "<<m_optimizer.objective()<<std::endl;
}


} // namespace gismo

