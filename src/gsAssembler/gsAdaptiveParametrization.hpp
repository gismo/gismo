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
        res(0,0) = dmeasdxi;
        res(1,0) = dmeasdeta;
        return res;

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


// NEW CODE
template<class T, enum MonitorMode MODE>
T gsOptMesh<T,MODE>::evalObj(const gsAsConstVector<T> &u) const
{
    typedef typename gsExprHelper<T>::geometryMap geometryMap; ///< Geometry map type
    gsExprEvaluator<T> evaluator;
    evaluator.options().setReal("quA",0.0);
    evaluator.options().setInt("quB",1);
    evaluator.setIntegrationElements(m_mb); // does not work when in constructor
    m_comp->setControls(u);

    // Penalty constant
    gsConstantFunction<T> pen(m_options.getReal("Penalty"), m_cgeom.domainDim());
    geometryMap G = evaluator.getMap(m_mp);
    auto eps = evaluator.getVariable(pen);
    if (m_cgeom.domainDim()==m_cgeom.targetDim())
    {
        auto detG = jac(G).det();
        auto chi = 0.5 * (detG + pow(pow(eps.val(),2.0) + pow(detG, 2.0), 0.5));
        auto invJacMat = jac(G).adj()/chi; // inverse of jacobian matrix with 'determinant' replaced

        if (m_fun==nullptr)
        {
            auto M = 1/detG;
            return evaluator.integral( (M*invJacMat).sqNorm()*meas(G));
        }
        else
        {
            // TEMPORARY FIX: SEE COMMENT IN CONSTRUCTOR
            m_cfun = m_parametric ? gsComposedFunction<T>(*m_comp,*m_fun) : gsComposedFunction<T>(m_cgeom,*m_fun);
            auto eta = evaluator.getVariable(m_cfun);
            auto M   = gismo::expr::monitor<MODE>(eta,G,m_options.getReal("Smoothing"));
            return evaluator.integral( (M*invJacMat).sqNorm()*meas(G));
        }
    }
    else if (m_cgeom.domainDim()<m_cgeom.targetDim())
    {
        auto fform = jac(G).tr()*jac(G);
        if (m_fun==nullptr)
        {
//          auto detG = pow(fform.det().val(),0.5).val(); //jacobian determinant for a surface, i.e. the measure
//          auto M = 1/detG;
//          return evaluator.integral( M*fform.trace()/pow(meas(G),2)*meas(G) );
//            return evaluator.integral( fform.trace()/pow(meas(G),2) );

            typename gsBasis<T>::domainIter domIt;
            gsQuadRule<T> QuRule;  // Quadrature rule
            gsMatrix<T> uvPoints; // Quadrature points in (u,v) space
            gsMatrix<T> xietaPoints; // Quadrature points in (xi,eta) space
            gsVector<T> tmpWeights;

            T result = 0.0;
            for (unsigned patchInd=0; patchInd < m_mb.nBases(); ++patchInd)
            {
                // Quadrature rule
                evaluator.options().setReal("quA",0.0);
                evaluator.options().setInt("quB",1);
                QuRule =  gsQuadrature::get(m_mb.basis(patchInd), evaluator.options());

                // Initialize domain element iterator
                domIt = m_mb.piece(patchInd).makeDomainIterator();
                for (; domIt->good(); domIt->next() )
                {
                    // Map the Quadrature rule to the element
                    QuRule.mapTo( domIt->lowerCorner(), domIt->upperCorner(),
                                 uvPoints, tmpWeights);
//                    gsDebugVar(tmpWeights);

                    // Evaluate the geometry map at the quadrature points
//                    m_comp->eval_into(uvPoints, xietaPoints);
                    for (index_t p = 0; p!=uvPoints.cols(); ++p)
                    {
                      gsDebugVar(evaluator.eval( meas(G), uvPoints.col(p) ).value());
                      result += evaluator.eval( meas(G), uvPoints.col(p) ).value() * tmpWeights(p);
//                      result += evaluator.eval( fform.trace()/pow(meas(G),1), uvPoints.col(p) ).value() * tmpWeights(p);
                    }
                }
            }


//          return evaluator.integral( fform.trace()/meas(G) );
            return result;
        }
        else
        {
            // TEMPORARY FIX: SEE COMMENT IN CONSTRUCTOR
            m_cfun = m_parametric ? gsComposedFunction<T>(*m_comp,*m_fun) : gsComposedFunction<T>(m_cgeom,*m_fun);
            auto eta = evaluator.getVariable(m_cfun);
            auto M   = gismo::expr::monitor<MODE>(eta,G,m_options.getReal("Smoothing"));
            return evaluator.integral( M*fform.trace()/meas(G)*meas(G));
        }
    }
    else
        GISMO_ERROR("Domain dimension must be smaller than or equal to the target dimension, but domainDim = "<<m_cgeom.domainDim()<<" and targetDim = "<<m_cgeom.targetDim());
}


// NEW CODE
template<class T, enum MonitorMode MODE>
void gsOptMesh<T,MODE>::gradObj_into ( const gsAsConstVector<T> & u, gsAsVector<T> & result) const
{
    result.resize(m_comp->nControls());
    result.setZero();

    typedef typename gsExprHelper<T>::geometryMap geometryMap;
    m_comp->setControls(u);

    gsExprEvaluator<T> evaluator_geom;
    gsMultiBasis<T> mmbasis(*m_geom);
    evaluator_geom.setIntegrationElements(mmbasis);
    geometryMap G_geom = evaluator_geom.getMap(*m_geom);

    const gsSquareDomain<T> * sqDomain = static_cast<const gsSquareDomain<T>*>(m_comp);

    typename gsBasis<T>::domainIter domIt;
    gsQuadRule<T> QuRule;
    gsMatrix<T> uvPoints, xietaPoints;
    gsVector<T> tmpWeights;
    gsMatrix<T> dSigmadAlpha, dDetJsigmadAlpha;

    for (unsigned patchInd=0; patchInd < m_mb.nBases(); ++patchInd)
    {
        gsExprEvaluator<T> evaluator;
        evaluator.options().setReal("quA",0.0);
        evaluator.options().setInt("quB",1);
        evaluator.setIntegrationElements(m_mb);
        QuRule = gsQuadrature::get(m_mb.basis(patchInd), evaluator.options());

        domIt = m_mb.piece(patchInd).makeDomainIterator();
        for (; domIt->good(); domIt->next() )
        {
            QuRule.mapTo( domIt->lowerCorner(), domIt->upperCorner(),
                         uvPoints, tmpWeights);

            m_comp->eval_into(uvPoints, xietaPoints);
            m_comp->control_deriv_into(uvPoints, dSigmadAlpha);
            sqDomain->control_jacobian_deriv_into(uvPoints, dDetJsigmadAlpha);

            for (index_t p = 0; p!=uvPoints.cols(); ++p)
            {
                gsAsMatrix<T> dSigmadAlphaMatrix(dSigmadAlpha.col(p).data(), m_comp->nControls(), 2);

                gsMatrix<T> jacobianSigma = m_comp->deriv(uvPoints.col(p));
                T detJsigma = jacobianSigma(0) * jacobianSigma(3) - jacobianSigma(1) * jacobianSigma(2);
                T absDetJsigma = math::abs(detJsigma);
                T signDetJsigma = (detJsigma > 0) ? T(1) : T(-1);

                gsMatrix<T> dEdSigma_val = evaluator_geom.eval(dEdSigma(G_geom), xietaPoints.col(p));

                T measG_at_sigma = evaluator_geom.eval(meas(G_geom), xietaPoints.col(p))(0);

                gsMatrix<T> term1 = dSigmadAlphaMatrix * dEdSigma_val * absDetJsigma;

                gsVector<T> term2 = measG_at_sigma * signDetJsigma * dDetJsigmadAlpha.col(p);

                result += tmpWeights[p] * (term1 + term2);
            }
        }
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

