/** @file monitor_composed_r-adaptivity.cpp

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
        //ones[0] = 4;
        res = 1.0 / ( math::sqrt(1.0 + grad.dot(grad) ) ) * ones;
        //res = ones;
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
        
        auto EfoldoverFree = (1.e-2 - jac(G).det()).ppartval();

        return m_evaluator.integral( (m(eta).asDiag()*invJacMat).sqNorm()*meas(G)
                                        + 1e4 * EfoldoverFree *meas(G)
                                        );
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


int main(int arg, char *argv[])
{
    //! [Parse command line]
    bool plot = false;
    index_t numRefine  = 2;
    index_t numElevate = 0;
    index_t maxIt = 100;
    real_t tol_g = 5e-5;

    gsCmdLine cmd("Tutorial on solving a Poisson problem.");
    cmd.addInt( "e", "degreeElevation",
                "Number of degree elevation steps to perform before solving (0: equalize degree in all directions)", numElevate );
    cmd.addInt( "r", "uniformRefine", "Number of Uniform h-refinement loops",  numRefine );
    cmd.addSwitch("plot", "Create a ParaView visualization file with the solution", plot);
    cmd.addReal("g", "tolG", "relative tol", tol_g);
    cmd.addInt( "i", "maxIt", "max num iteterations",  maxIt );


    try { cmd.getValues(arg,argv); } catch (int rv) { return rv; }
    //! [Parse command line]


    gsKnotVector<> kv({0,0,0,0.25,0.5,0.75,1,1,1},2);
    gsTensorBSplineBasis<2> tbasis(kv,kv);

    gsMatrix<> coefs = tbasis.anchors().transpose();
    // coefs.row(coefs.rows()-1) *= 1.2;
    gsTensorBSpline<2> tbspline(tbasis,coefs);

    gsSquareDomain<2,real_t> domain(numElevate,numRefine);
    domain.options().addSwitch("Slide","",true);

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

    gsOptMesh<> optMesh(domain,tbspline,fun);
    gsVector<> controls(domain.nControls());
     for (size_t k=0; k!=domain.nControls(); k++)
        controls[k] = domain.control(k);
         // controls[k] = domain.control(k) * 0.9;


    // gsDebugVar(controls.transpose());
    
    gsHLBFGS<real_t> optimizer(&optMesh);
    optimizer.options().setInt("MaxIterations",maxIt);
    optimizer.options().setInt("Verbose",2);
    optimizer.options().setReal("tolRelG",tol_g);
    // gsDebugVar(optimizer.currentDesign().transpose());
    optimizer.solve(controls);
    



    gsVector<> solution = optimizer.currentDesign();
    gsDebugVar(solution.transpose());
    

    for (size_t k=0; k!=solution.rows(); k++)
        domain.control(k) = solution[k];



    gsWriteParaview(cspline,"cspline",1000);
    gsWriteParaview(domain,"domain",1000);


    

    return 0;
}// end main
