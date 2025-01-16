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
                                const gsBasis<T>    * integrationBasis,
                                const bool            parametric)
:
gsOptMesh(composition,geometry,nullptr,integrationBasis,parametric)
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
m_mb(*m_ib),
m_cgeom(*m_comp,geometry),
// THIS DOES NOT WORK FOR PARAMETRIC=FALSE. POINTER IS LOST WHEN ARRIVING IN evalObj
// m_cfun(parametric ? gsComposedFunction<T>(*m_comp,fun) : gsComposedFunction<T>(m_cgeom,fun)),
m_mp(m_cgeom),
m_parametric(parametric)
{
    m_numDesignVars = m_comp->nControls();
    m_curDesign.resize(m_numDesignVars,1);
    m_controls.resize(m_numDesignVars,1);
    for (index_t k=0; k!=m_numDesignVars; k++)
        m_controls(k,0) = m_comp->control(k);

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

template<class T, enum MonitorMode MODE>
T gsOptMesh<T,MODE>::evalObj(const gsAsConstVector<T> &u) const
{
    typedef typename gsExprHelper<T>::geometryMap geometryMap; ///< Geometry map type
    gsExprEvaluator<T> evaluator;

    evaluator.setIntegrationElements(m_mb); // does not work when in constructor
    for (index_t k=0; k!=m_numDesignVars; k++)
        m_comp->control(k) = u[k];

    // Penalty constant
    gsConstantFunction<T> pen(m_options.getReal("Penalty"), m_cgeom.domainDim());

    geometryMap G = evaluator.getMap(m_mp);
    auto eps = evaluator.getVariable(pen);

    // auto chi = 0.5 * (jac(G).det() + pow(eps.val() + pow(jac(G).det(), 2.0), 0.5));
    // auto invJacMat = jac(G).adj()/chi;
    // auto eta = evaluator.getVariable(m_fun);
    // return evaluator.integral( (monitor(eta,G).asDiag()*invJacMat).sqNorm()*meas(G));

    // gsComposedFunction<T> fun(*m_comp,m_fun.function(0));
    // gsComposedFunction<T> fun(mpG.patch(0),m_fun.function(0));


    // gsVector<unsigned> np(m_fun.domainDim());
    // np.setConstant(m_options.getInt("nSamplingPoints"));
    // gsMatrix<T> grid = gsPointGrid<T>(m_cgeom.support().col(0),m_cgeom.support().col(1),np);
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
        auto detG = pow(fform.det().val(),0.5); //jacobian determinant for a surface, i.e. the measure
        // Compute the chi part
        auto chiPPart = eps * ((detG.val() - eps.val()).exp());
        // Ternary operation to compute chi and chip
        auto chi = ternary(eps.val() - detG, chiPPart.val(), detG.val());
        auto invJacMat = fform.sqrt().adj()/chi; // inverse of jacobian matrix with 'determinant' replaced

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
    else
        GISMO_ERROR("Domain dimension must be smaller than or equal to the target dimension, but domainDim = "<<m_cgeom.domainDim()<<" and targetDim = "<<m_cgeom.targetDim());
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
    if (const gsTensorBSplineBasis<2,T> * tbasis = dynamic_cast<const gsTensorBSplineBasis<2,T> *>(&integrationBasis))
    {
        gsTensorBSplineBasis<2,T> ibasis(*tbasis);
        // Integration basis: parent basis with knots of composition basis inserted, and the degree is the sum of the two degrees (?)
        index_t targetDegree;
        const gsTensorBSplineBasis<2,T> & comp_tbasis = dynamic_cast<const gsTensorBSplineBasis<2,T> &>(composition.domain().basis());
        for (size_t d = 0; d!=2; d++)
        {
            // 1. Insert interior knots of composition basis
            for (typename gsKnotVector<real_t>::uiterator it = std::next(comp_tbasis.knots(d).ubegin());
                                                        it!= std::prev(comp_tbasis.knots(d).uend());
                                                        ++it)
                {
                    if (ibasis.knots(d).has(*it))
                        continue;
                    ibasis.insertKnot(*it,d);
                }
            // 2. Increase the degree
            targetDegree = ibasis.degree(d) * comp_tbasis.degree(d);
            ibasis.degreeIncrease(targetDegree-ibasis.degree(d),d);

            m_integrationBasis = memory::make_unique(new gsTensorBSplineBasis<2,T>(ibasis));
        }
    }
    m_optProblem = gsOptMesh<T,MODE>(m_comp,m_geom,m_fun,m_integrationBasis.get(),parametric);
    this->defaultOptions();
}

template <class T, enum MonitorMode MODE>
gsOptionList & gsAdaptiveParametrization<T,MODE>::options()
{
    return m_optProblem.options();
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
    gsVector<T> controls(m_optProblem.composition().nControls());
    for (size_t k=0; k!=m_optProblem.composition().nControls(); k++)
        controls[k] = m_optProblem.composition().control(k);

    m_optimizer.solve(controls);
    controls = m_optimizer.currentDesign();
    // Write the optimized solution back to the composition
    for (index_t k=0; k!=controls.rows(); k++)
        m_optProblem.composition().control(k) = controls[k];
}


} // namespace gismo

