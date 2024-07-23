/** @file composed_domain_L2.cpp

    @brief Tutorial on how to use expression assembler to solve the Poisson equation

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris
*/

//! [Include namespace]
#include <gismo.h>
#include <gsNurbs/gsSquareDomain.h>
#include <gsHLBFGS/gsHLBFGS.h>
#include <gsModeling/gsBarrierCore.h>

//! [Include namespace]

namespace gismo{
namespace expr{

template<class E>
class monitor_expr : public _expr<monitor_expr<E> >
{
public:
    typedef typename E::Scalar Scalar;

    mutable gsMatrix<Scalar> res;
    mutable gsMatrix<Scalar> grad;
    mutable gsVector<Scalar> ones;

private:
    typename E::Nested_t _u;

public:
    enum{ Space = E::Space, ScalarValued= 0, ColBlocks= 0};

    monitor_expr(const E & u) : _u(u)
    {
    }

#   define Eigen gsEigen
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
#   undef Eigen

    const gsMatrix<Scalar> & eval(const index_t k) const
    {return eval_impl(_u,k); }

    index_t rows() const { return _u.source().domainDim(); }

    index_t cols() const { return 1; }

    void parse(gsExprHelper<Scalar> & evList) const
    {
        evList.add(_u);
        _u.data().flags |= NEED_GRAD;
    }

    const gsFeSpace<Scalar> & rowVar() const { return _u.rowVar(); }
    const gsFeSpace<Scalar> & colVar() const {return gsNullExpr<Scalar>::get();}
    index_t cardinality_impl() const { return _u.cardinality_impl(); }

    void print(std::ostream &os) const { os << "M("; _u.print(os); os <<")"; }

private:
    template<class U> inline
    typename util::enable_if< util::is_same<U,gsFeSpace<Scalar> >::value, const gsMatrix<Scalar> & >::type
    eval_impl(const U & u, const index_t k)  const
    {
        // const index_t A = _u.cardinality()/_u.dim(); // _u.data().actives.rows()
        // res.resize(_u.cardinality(), cols()); // rows()*

        // normal = _G.data().normal(k);// not normalized to unit length
        // normal.normalize();
        // bGrads = _u.data().values[1].col(k);
        // cJac = _G.data().values[1].reshapeCol(k, _G.data().dim.first, _G.data().dim.second).transpose();
        // const Scalar measure =  (cJac.col3d(0).cross( cJac.col3d(1) )).norm();

        // for (index_t d = 0; d!= cols(); ++d) // for all basis function components
        // {
        //     const short_t s = d*A;
        //     for (index_t j = 0; j!= A; ++j) // for all actives
        //     {
        //         // Jac(u) ~ Jac(G) with alternating signs ?..
        //         m_v.noalias() = (vecFun(d, bGrads.at(2*j  ) ).cross( cJac.col(1).template head<3>() )
        //                       - vecFun(d, bGrads.at(2*j+1) ).cross( cJac.col(0).template head<3>() )) / measure;

        //         // ---------------  First variation of the normal
        //         // res.row(s+j).noalias() = (m_v - ( normal.dot(m_v) ) * normal).transpose();
        //         res.row(s+j).noalias() = (m_v - ( normal*m_v.transpose() ) * normal).transpose(); // outer-product version
        //     }
        // }
        // return res;
    }

    template<class U> inline
    typename util::enable_if< util::is_same<U,gsFeSolution<Scalar> >::value, const gsMatrix<Scalar> & >::type
    eval_impl(const U & u, const index_t k)  const
    {
        GISMO_ASSERT(1==_u.data().actives.cols(), "Single actives expected");
        grad_expr<gsFeSolution<Scalar>> sGrad =  grad_expr<gsFeSolution<Scalar>>(_u);
        grad = sGrad.eval(k);
        ones.resize( _u.source().domainDim());
        ones.setOnes();
        ones[0] = 4;
        // res = 1.0 / ( math::sqrt(1.0 + grad.dot(grad) ) ) * ones;
        res = ones;
        return res;
    }

    template<class U> inline
    typename util::enable_if< !util::is_same<U,gsFeSolution<Scalar> >::value, const gsMatrix<Scalar> & >::type
    eval_impl(const U & u, const index_t k)  const
    {
        grad_expr<U> vGrad = grad_expr<U>(u);
        grad = vGrad.eval(k);
        ones.resize( _u.source().domainDim());
        ones.setOnes();

        res.resize( _u.source().domainDim(), 1);

        res = 1.0 / ( math::sqrt(1.0 + grad.norm()*grad.norm() ) ) * ones;
        // res(0,0) = 1.0 / ( math::sqrt(1.0 + (grad*grad.transpose()).value() ) );
        // res(1,0) = 1;
        return res;
    }
};

template<class E> EIGEN_STRONG_INLINE
monitor_expr<E> m(const E & u) { return monitor_expr<E>(u); }


}
}

using namespace gismo;


template<typename T = real_t>
class gsOptMesh : public gsOptProblem<T> {

private:
    typedef typename gsExprAssembler<T>::geometryMap geometryMap;
    typedef typename gsExprAssembler<T>::space space;
    typedef typename gsExprAssembler<T>::solution solution;

public:
    gsOptMesh(  gsFunction<T> & composition,
                const gsGeometry<T> & geometry,
                const gsFunction<T> & fun)
    :
    m_comp(&composition),
    m_geom(geometry),
    m_fun(fun)
    {}

    /// Evaluates the objective function at the given point u.
    T evalObj(const gsAsConstVector<T> &u) const override
    {
        for (index_t k=0; k!=u.rows(); k++)
            m_comp->control(k) = u[k];

        gsComposedGeometry<T> cgeom(*m_comp,m_geom);

        gsMultiPatch<> mp;
        mp.addPatch(cgeom);
        gsMultiBasis<> mb(m_geom.basis());

        m_evaluator.setIntegrationElements(mb);
        geometryMap G = m_evaluator.getMap(mp);

        auto invJacMat = jac(G).inv();
        auto eta = m_evaluator.getVariable(m_fun,G); // This G might not be needed for a spline-based error field
        return m_evaluator.integral((m(eta).asDiag()*invJacMat).sqNorm()*meas(G));
    }

    // /// Computes the gradient of the objective function at the given point u
    // // and stores it in result.
    // void gradObj_into(const gsAsConstVector<T> &u,
    //                 gsAsVector<T> &result) const override;

protected:
    gsFunction<T> * m_comp;
    const gsGeometry<T> & m_geom;
    const gsFunction<T> & m_fun;

    mutable gsExprEvaluator<T> m_evaluator;
};

template<typename T = real_t>
class gsOptMesh2 : public gsOptProblem<T> {

private:
    typedef typename gsExprAssembler<T>::geometryMap geometryMap;
    typedef typename gsExprAssembler<T>::space space;
    typedef typename gsExprAssembler<T>::solution solution;

public:
    gsOptMesh2( gsGeometry<T> & geometry,
                const gsFunction<T> & fun)
    :
    m_geom(geometry),
    m_fun(fun)
    {
        // Mapper storing control points
        m_mapper = gsDofMapper(m_geom.basis(),m_geom.targetDim());

        gsBoxTopology topology(m_geom.domainDim(),1);
        topology.addAutoBoundaries();
        // gsMatrix<index_t> boundary = m_geom.basis().allBoundary();
        for (typename gsBoxTopology::biterator it = topology.bBegin(); it != topology.bEnd(); ++it)
        {
            gsMatrix<index_t> boundary = m_geom.basis().boundary(*it);
            for (index_t d = 0; d!=m_geom.targetDim(); d++)
                m_mapper.markBoundary(0,boundary,d);
        }
        m_mapper.finalize();

    }

    gsDofMapper mapper() { return m_mapper;}

    /// Evaluates the objective function at the given point u.
    T evalObj(const gsAsConstVector<T> &u) const override
    {
        for (short_t d=0; d!=m_geom.domainDim(); d++)
            for (index_t k=0; k!=m_geom.coefs().rows(); k++)
                if (m_mapper.is_free(k,0,d))
                    m_geom.coefs()(k,d) = u[m_mapper.index(k, 0, d)];

        gsMultiPatch<> mp;
        mp.addPatch(m_geom);
        gsMultiBasis<> mb(m_geom.basis());

        m_evaluator.setIntegrationElements(mb);
        geometryMap G = m_evaluator.getMap(mp);

        auto invJacMat = jac(G).inv();
        auto eta = m_evaluator.getVariable(m_fun,G);
        return m_evaluator.integral((m(eta).asDiag()*invJacMat).sqNorm() * meas(G));
    }

    // /// Computes the gradient of the objective function at the given point u
    // // and stores it in result.
    // void gradObj_into(const gsAsConstVector<T> &u,
    //                 gsAsVector<T> &result) const override;

protected:

    gsFunction<T> * m_comp;
    gsGeometry<T> & m_geom;
    const gsFunction<T> & m_fun;

    gsDofMapper m_mapper;

    gsMatrix<T> m_free;

    mutable gsExprEvaluator<T> m_evaluator;
};


int main(int arg, char *argv[])
{
    //! [Parse command line]
    bool plot = false;
    index_t numRefine  = 2;
    index_t numElevate = 0;

    gsCmdLine cmd("Tutorial on solving a Poisson problem.");
    cmd.addInt( "e", "degreeElevation",
                "Number of degree elevation steps to perform before solving (0: equalize degree in all directions)", numElevate );
    cmd.addInt( "r", "uniformRefine", "Number of Uniform h-refinement loops",  numRefine );
    cmd.addSwitch("plot", "Create a ParaView visualization file with the solution", plot);

    try { cmd.getValues(arg,argv); } catch (int rv) { return rv; }
    //! [Parse command line]


    gsKnotVector<> kv({0,0,0,0.25,0.5,0.75,1,1,1},2);
    gsTensorBSplineBasis<2> tbasis(kv,kv);

    gsMatrix<> coefs = tbasis.anchors().transpose();
    // coefs.row(coefs.rows()-1) *= 1.2;
    gsTensorBSpline<2> tbspline(tbasis,coefs);

    gsSquareDomain<2,real_t> domain(1,2);

    gsComposedGeometry<real_t> cspline(domain,tbspline);


    gsMultiPatch<> mp;
    mp.addPatch(cspline);
    gsMultiBasis<> mb(mp);


/*
        std::string k = "100";
        std::string w = "0.1";
        std::string l = "0.2";
        gsFunctionExpr<> fun("-(1/2 + tanh("+k+"*(y - 0.5 + "+w+"/2))/2)*(1/2 + tanh("+k+"*(0.5 + "+w+"/2 - y))/2)*(-"+k+"^2*tanh("+k+"*(x - 1 + "+l+"))*(1 - tanh("+k+"*(x - 1 + "+l+"))^2) - "+k+"^2*tanh("+k+"*("+l+" - x))*(1 - tanh("+k+"*("+l+" - x))^2)) + "+k+"^2*tanh("+k+"*(x - 0.5 + "+w+"/2))*(1 - tanh("+k+"*(x - 0.5 + "+w+"/2))^2)*(1/2 + tanh("+k+"*(0.5 + "+w+"/2 - x))/2)*(1 + tanh("+k+"*(y - 1 + "+l+"))/2 + tanh("+k+"*("+l+" - y))/2) + "+k+"^2*(1 - tanh("+k+"*(x - 0.5 + "+w+"/2))^2)*(1 - tanh("+k+"*(0.5 + "+w+"/2 - x))^2)*(1 + tanh("+k+"*(y - 1 + "+l+"))/2 + tanh("+k+"*("+l+" - y))/2)/2 + (1/2 + tanh("+k+"*(x - 0.5 + "+w+"/2))/2)*"+k+"^2*tanh("+k+"*(0.5 + "+w+"/2 - x))*(1 - tanh("+k+"*(0.5 + "+w+"/2 - x))^2)*(1 + tanh("+k+"*(y - 1 + "+l+"))/2 + tanh("+k+"*("+l+" - y))/2) + "+k+"^2*tanh("+k+"*(y - 0.5 + "+w+"/2))*(1 - tanh("+k+"*(y - 0.5 + "+w+"/2))^2)*(1/2 + tanh("+k+"*(0.5 + "+w+"/2 - y))/2)*(1 + tanh("+k+"*(x - 1 + "+l+"))/2 + tanh("+k+"*("+l+" - x))/2) + "+k+"^2*(1 - tanh("+k+"*(y - 0.5 + "+w+"/2))^2)*(1 - tanh("+k+"*(0.5 + "+w+"/2 - y))^2)*(1 + tanh("+k+"*(x - 1 + "+l+"))/2 + tanh("+k+"*("+l+" - x))/2)/2 + (1/2 + tanh("+k+"*(y - 0.5 + "+w+"/2))/2)*"+k+"^2*tanh("+k+"*(0.5 + "+w+"/2 - y))*(1 - tanh("+k+"*(0.5 + "+w+"/2 - y))^2)*(1 + tanh("+k+"*(x - 1 + "+l+"))/2 + tanh("+k+"*("+l+" - x))/2) - (1/2 + tanh("+k+"*(x - 0.5 + "+w+"/2))/2)*(1/2 + tanh("+k+"*(0.5 + "+w+"/2 - x))/2)*(-"+k+"^2*tanh("+k+"*(y - 1 + "+l+"))*(1 - tanh("+k+"*(y - 1 + "+l+"))^2) - "+k+"^2*tanh("+k+"*("+l+" - y))*(1 - tanh("+k+"*("+l+" - y))^2))",2);
*/


    gsFunctionExpr<> fun("tanh(30*((x)^2+(y)^2-1.0/16.0))",2);
    gsWriteParaview(mp,fun,"fun");

    // gsWriteParaview(mp,domain,"domain",1000,true,true);


    gsExprEvaluator<> ev;
    ev.setIntegrationElements(mb);
    gsExprEvaluator<>::geometryMap G = ev.getMap(mp);

    auto invJacMat = jac(G).inv();
    auto u = ev.getVariable(fun);


    // gsMatrix<> pts(2,3);
    // pts.col(0).setConstant(0.3);
    // pts.col(1).setConstant(0.6);
    // pts.col(2).setConstant(0.9);
    // gsDebugVar(ev.eval(m(u).asDiag(),pts));


    ev.writeParaview(m(u),G,"m");


    gsDebugVar(ev.integral((m(u).asDiag()*invJacMat).sqNorm() * meas(G)));

    // gsOptMesh<> optMesh(domain,tbspline,fun);
    gsOptMesh2<> optMesh(tbspline,fun);
    gsMatrix<> initial = convertMultiPatchToFreeVector(mp,optMesh.mapper());
    initial *= 0.8;

    // gsVector<> controls(domain.nControls());
    // for (size_t k=0; k!=domain.nControls(); k++)
    //     controls[k] = domain.control(k)*0.8;

    // gsDebugVar(controls.transpose());

    // gsDebugVar(optMesh.evalObj(gsAsConstVector<>(controls.data(),controls.rows())));
    // gsDebugVar(domain.nControls());
    gsHLBFGS<real_t> optimizer(&optMesh);
    optimizer.options().setInt("MaxIterations",200);
    optimizer.options().setInt("Verbose",2);
    // gsDebugVar(optimizer.currentDesign().transpose());
    // optimizer.solve(controls);
    optimizer.solve(initial);



    gsVector<> solution = optimizer.currentDesign();
    gsDebugVar(solution.transpose());
    convertFreeVectorToMultiPatch(solution,optMesh.mapper(),mp);
    gsWriteParaview(mp,"mp",1000);


    // gsDebugVar(solution.transpose());
    // for (size_t k=0; k!=solution.rows(); k++)
    //     domain.control(k) = solution[k];



    // gsWriteParaview(cspline,"cspline",1000);
    // gsWriteParaview(domain,"domain",1000);


    // for (index_t k=0; k!=domain.nControls(); k++)
    //     domain.control(k) *= 0.5;

    return 0;
}// end main
