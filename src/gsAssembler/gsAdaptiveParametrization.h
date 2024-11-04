/** @file gsAdaptiveParametrization.h

    @brief Provides class for adaptive parametrization (r-adaptivity).

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst (TU Delft, 2019-)
*/

#pragma once

#include <gismo.h>
// #include <gsCore/gsComposedFunction.h>
#include <gsNurbs/gsSquareDomain.h>
#include <gsHLBFGS/gsHLBFGS.h>
#include <gsModeling/gsBarrierCore.h>
#include <gsOptimizer/gsGradientDescent.h>
#include <gsOptim/gsOptim.h>

enum MonitorMode
{
    // We can add the option "Runtime" where the mode is determined at runtime using an option
    ValueBased = 0,
    GradientBased = 1
};

namespace gismo
{
namespace expr
{

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


template<typename T, enum MonitorMode MODE>
class gsOptMesh : public gsOptProblem<T>
{
    using Base = gsOptProblem<T>;

private:
    typedef typename gsExprAssembler<T>::geometryMap geometryMap;
    typedef typename gsExprAssembler<T>::space space;
    typedef typename gsExprAssembler<T>::solution solution;

public:

    gsOptMesh() {}

    gsOptMesh(        gsFunction<T> & composition,
                const gsGeometry<T> & geometry,
                const gsFunction<T> & fun,
                const gsBasis<T>    * integrationBasis,
                const bool            parametric)
    :
    m_comp(&composition),
    m_ib(integrationBasis),
    m_mb(*m_ib),
    m_cgeom(*m_comp,geometry),
    m_cfun(parametric ? gsComposedFunction<T>(*m_comp,fun) : gsComposedFunction<T>(m_cgeom,fun)),
    m_mp(m_cgeom)
    {
        m_numDesignVars = m_comp->nControls();
        m_curDesign.resize(m_numDesignVars,1);
        m_controls.resize(m_numDesignVars,1);
        for (index_t k=0; k!=m_numDesignVars; k++)
            m_controls(k,0) = m_comp->control(k);

        m_options.addReal("Smoothing","Smoothing parameter for the monitor function",0.1);
        m_options.addReal("Penalty","Penalty parameter for the monitor function",1e-2);
    }

    gsFunction<T> & composition() { return *m_comp; }

    gsOptionList & options() { return m_options; }

    /// Evaluates the objective function at the given point u.
    T evalObj(const gsAsConstVector<T> &u) const override
    {
        m_evaluator.setIntegrationElements(m_mb); // does not work when in constructor
        for (index_t k=0; k!=m_numDesignVars; k++)
            m_comp->control(k) = u[k];

        // Penalty constant
        gsConstantFunction<T> pen(m_options.getReal("Penalty"), m_cgeom.domainDim());

        geometryMap G = m_evaluator.getMap(m_mp);
        auto eps = m_evaluator.getVariable(pen);

        // auto chi = 0.5 * (jac(G).det() + pow(eps.val() + pow(jac(G).det(), 2.0), 0.5));
        // auto invJacMat = jac(G).adj()/chi;
        // auto eta = m_evaluator.getVariable(m_fun);
        // return m_evaluator.integral( (monitor(eta,G).asDiag()*invJacMat).sqNorm()*meas(G));

        // gsComposedFunction<T> fun(*m_comp,m_fun.function(0));
        // gsComposedFunction<T> fun(mpG.patch(0),m_fun.function(0));


        // gsVector<unsigned> np(m_fun.domainDim());
        // np.setConstant(m_options.getInt("nSamplingPoints"));
        // gsMatrix<T> grid = gsPointGrid<T>(m_cgeom.support().col(0),m_cgeom.support().col(1),np);
        if (m_cgeom.domainDim()==m_cgeom.targetDim())
        {
            auto jacG = jac(G).det();
            auto chi = 0.5 * (jacG + pow(pow(eps.val(),2.0) + pow(jacG, 2.0), 0.5));
            auto invJacMat = jac(G).adj()/chi; // inverse of jacobian matrix with 'determinant' replaced

            // auto eta = m_evaluator.getVariable(fun);
            auto eta = m_evaluator.getVariable(m_cfun);
            // auto eta = m_evaluator.getVariable(m_fun,G);

            auto M   = gismo::expr::monitor<MODE>(eta,G,m_options.getReal("Smoothing"));
            return m_evaluator.integral( (M*invJacMat).sqNorm()*meas(G));
        }
        else
        {
            // auto fform = jac(G).tr()*jac(G);
            // auto jacG = sqrt(fform.det()); //jacobian determinant for a surface, i.e. the measure
            // // auto chi = 0.5 * (jacG + pow(pow(eps.val(),2.0) + pow(jacG, 2.0), 0.5));

            // // Compute the chi part
            // auto chiPPart = eps * ((jacG - eps.val()).exp());

            // // Ternary operation to compute chi and chip
            // auto chi = ternary(eps.val() - jacG, chiPPart.val(), jacG.val());

            // auto invJacMat = fform.adj()/chi; // inverse of jacobian matrix with 'determinant' replaced
            // auto eta = m_evaluator.getVariable(fun);
            // return m_evaluator.integral( (monitor(eta,G)*invJacMat).sqNorm()*meas(G));

            GISMO_ERROR("The dimension of target domain should be 2 or 3.");
            return 0;
        }
    }

    // /// Computes the gradient of the objective function at the given point u
    // // and stores it in result.
    // void gradObj_into(const gsAsConstVector<T> &u,
    //                 gsAsVector<T> &result) const override;

protected:

    // From gsOptProblem
    using Base::m_curDesign;
    using Base::m_numDesignVars;

    gsFunction<T>             * m_comp;
    const gsBasis<T>          * m_ib;
    gsMultiBasis<T>             m_mb;
    gsComposedGeometry<T>       m_cgeom;
    gsComposedFunction<T>       m_cfun;
    gsMultiPatch<T>             m_mp;


    // Controls of the composition
    // NOTE: Different from m_curDesign, since m_curDesign updates every time
    gsMatrix<T> m_controls;

    gsOptionList                m_options;

    mutable gsExprEvaluator<T>  m_evaluator;

};


/*
    NOTES:
    * Give the integration basis?? The integration elements should be the finest of the elements of the composition and the basis/geometry
    * There should be a rule of thumb for the number of integration points. Given a composition of degree p and a geometry of degree q, the number of integration points should be at least p*q+1?? or p+q+1??
 */
template<class T, enum MonitorMode MODE=MonitorMode::ValueBased>
class gsAdaptiveParametrization
{
protected:

public:

    /**
     * @brief Constructs a gsAdaptiveParametrization object.
     *
     * This constructor takes any composition and integration basis among other inputs, and computes integrals based on the integration rule on the integration basis.
     *
     * @param composition       a \a gsFunction object representing the composition of the parametrization
     * @param geometry          a \a gsGeometry object representing the geometry
     * @param function          a \a gsFunction object representing the indicator function
     * @param integrationBasis  a \a gsBasis object used to define the integration points
     * @param optimizer         a \a gsOptimizer object used to solve the optimization problem
     * @param parametric        a boolean indicating whether the composition \a function is defined in the parametric domain (default = true)
     */
    gsAdaptiveParametrization(        gsFunction<T>  & composition,
                                const gsGeometry<T>  & geometry,
                                const gsFunction<T>  & function,
                                const gsBasis<T>     & integrationBasis,
                                      gsOptimizer<T> & optimizer,
                                const bool             parametric=true)
    :
    m_optimizer(optimizer),
    m_integrationBasis(integrationBasis.clone())//,
    {
        m_optProblem = gsOptMesh<T,MODE>(composition,geometry,function,m_integrationBasis.get(),parametric);
        this->defaultOptions();
    }

    /**
     * @brief Constructs a gsAdaptiveParametrization object.
     *
     * This constructor takes any composition among other inputs, and computes integrals based on the integration rule on the integration basis.
     *
     * @param composition       a \a gsFunction object representing the composition of the parametrization
     * @param geometry          a \a gsGeometry object representing the geometry
     * @param function          a \a gsFunction object representing the indicator function
     * @param optimizer         a \a gsOptimizer object used to solve the optimization problem
     * @param parametric        a boolean indicating whether the composition \a function is defined in the parametric domain (default = true)
     */
    gsAdaptiveParametrization(        gsFunction<T>  & composition,
                                const gsGeometry<T>  & geometry,
                                const gsFunction<T>  & function,
                                      gsOptimizer<T> & optimizer,
                                const bool             parametric=true)
    :
    gsAdaptiveParametrization(composition,geometry,function,geometry.basis(),optimizer,parametric)
    {
    }

    /**
     * @brief Constructs a gsAdaptiveParametrization object.
     *
     * This constructor takes a \ref gsSquareDomain as composition and integration basis among other inputs, and computes integrals based on the union of the integration basis and the composition
     *
     * @param composition       a \a gsSquareDomain object representing the composition of the parametrization. This object will be exactly integrated by inserting its integration points in the integration basis.
     * @param geometry          a \a gsGeometry object representing the geometry
     * @param function          a \a gsFunction object representing the indicator function
     * @param integrationBasis  a \a gsBasis object used to define the integration points. The integration points from the composition will be added to this basis.
     * @param optimizer         a \a gsOptimizer object used to solve the optimization problem
     * @param parametric        a boolean indicating whether the composition \a function is defined in the parametric domain (default = true)
     */
    template <short_t DIM>
    gsAdaptiveParametrization(        gsSquareDomain<DIM,T> & composition,
                                const gsGeometry<T>         & geometry,
                                const gsFunction<T>         & function,
                                const gsBasis<T>            & integrationBasis,
                                      gsOptimizer<T>        & optimizer,
                                const bool                    parametric=true)
    :
    m_optimizer(optimizer)
    {
        if (const gsTensorBSplineBasis<DIM,T> * tbasis = dynamic_cast<const gsTensorBSplineBasis<DIM,T> *>(&integrationBasis))
        {
            gsTensorBSplineBasis<DIM,T> ibasis(*tbasis);
            // Integration basis: parent basis with knots of composition basis inserted, and the degree is the sum of the two degrees (?)
            index_t targetDegree;
            const gsTensorBSplineBasis<DIM,T> & comp_tbasis = dynamic_cast<const gsTensorBSplineBasis<DIM,T> &>(composition.domain().basis());
            for (size_t d = 0; d!=DIM; d++)
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

                m_integrationBasis = memory::make_unique(new gsTensorBSplineBasis<DIM,T>(ibasis));
            }
        }
        m_optProblem = gsOptMesh<T,MODE>(composition,geometry,function,m_integrationBasis.get(),parametric);
        this->defaultOptions();
    }

    /**
     * @brief Constructs a gsAdaptiveParametrization object.
     *
     * This constructor takes a \ref gsSquareDomain as composition among other inputs, and computes integrals based on the union of the integration basis and the composition
     * The integration basis is set to the geometry's basis.
     *
     * @param composition       a \a gsSquareDomain object representing the composition of the parametrization. This object will be exactly integrated by inserting its integration points in the integration basis.
     * @param geometry          a \a gsGeometry object representing the geometry
     * @param function          a \a gsFunction object representing the indicator function
     * @param optimizer         a \a gsOptimizer object used to solve the optimization problem
     * @param parametric        a boolean indicating whether the composition \a function is defined in the parametric domain (default = true)
     */
    template <short_t DIM>
    gsAdaptiveParametrization(        gsSquareDomain<DIM,T> & composition,
                                const gsGeometry<T>         & geometry,
                                const gsFunction<T>         & function,
                                      gsOptimizer<T>        & optimizer,
                                const bool                    parametric=true)
    :
    gsAdaptiveParametrization(composition,geometry,function,geometry.basis(),optimizer,parametric)
    {
    }

    gsOptionList & options() { return m_options; }

    void defaultOptions()
    {
        m_options.addReal("Penalty","Penalization coefficient for Jacobian determinant [default=0.01]",1e-2);
        m_options.addReal("Smoothing","Smoothing parameter in the monitor function [default=0.1]",0.1);
        m_options.addInt("Mode","0: Relocate based on f [default]; 1: Relocate based on grad(f)",0);
    }

    void solve()
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

protected:

    gsOptMesh<T,MODE>           m_optProblem;
    gsOptimizer<T>            & m_optimizer;
    typename gsBasis<T>::uPtr   m_integrationBasis;
    gsOptionList                m_options;
};


} // namespace gismo

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsAdaptiveParametrization.hpp)
#endif