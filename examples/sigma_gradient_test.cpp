/** @file compositions_example.cpp

    @brief Tutorial on how to generate a composed basis and a composed geometry

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris
*/

//! [Include namespace]
#include <gismo.h>
#include <gsNurbs/gsSquareDomain.h>
#include <gsAssembler/gsAdaptiveParametrization.h>

namespace gismo
{
namespace expr
{

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


}
}

using namespace gismo;
//! [Include namespace]

template <class T>
class SigmaFunction : public gsFunction<T>
{

public:
    typedef gsFunction<T> Base;

    SigmaFunction(const gsSquareDomain<T> & domain, gsVector<T> & pt)
    :
    m_domain(domain),
    m_pt(pt)
    {
    }


    void eval_into(const gsMatrix<T> & controls, gsMatrix<T> & result) const override
    {
        result.resize(m_domain.targetDim(), controls.cols());
        gsMatrix<T> tmp;
        for (index_t k=0; k!=controls.cols(); k++)
        {
            m_domain.setControls(controls.col(k));
            m_domain.eval_into(m_pt, tmp);
            result.col(k) = tmp;
        }
    }

    void an_deriv_into(const gsMatrix<T> & controls, gsMatrix<T> & result) const
    {
        result.resize(m_domain.targetDim()*m_domain.nControls(), controls.cols());
        gsMatrix<T> tmp;
        for (index_t k=0; k!=controls.cols(); k++)
        {
            m_domain.setControls(controls.col(k));
            m_domain.control_deriv_into(m_pt, tmp);
            result.col(k) = tmp;
        }
    }

    short_t domainDim() const override
    {
        return m_domain.nControls();
    }

    short_t targetDim() const override
    {
        return m_domain.targetDim();
    }

protected:
    mutable gsSquareDomain<T> m_domain; ///< The multi-patch to evaluate
    gsVector<T> m_pt; ///< The point at which to evaluate the function
};

template <class T>
class Efunction : public gsFunction<T>
{

public:
    typedef gsFunction<T> Base;

    Efunction(const gsMultiPatch<T> & mp)
    :
    m_mp(mp),
    m_mb(mp)
    {
        m_evaluator.setIntegrationElements(m_mb);
    }

    void eval_into(const gsMatrix<T> & u, gsMatrix<T> & result) const override
    {
        result.resize(1,u.cols());
        auto G = m_evaluator.getMap(m_mp);
        auto fform = jac(G).tr()*jac(G);
        auto detG = pow(fform.det().val(),0.5).val(); //jacobian determinant for a surface, i.e. the measure
        auto M = 1/detG;
        for (index_t k=0; k!=u.cols(); k++)
            // ONLY MEAS
            result(0,k) = m_evaluator.eval( pow(meas(G),1),u.col(k)).value();
            // ONLY FFORM
            // result(0,k) = m_evaluator.eval( fform.trace(),u.col(k)).value();
            // FULL EXPR
            // result(0,k) = m_evaluator.eval( M*fform.trace()/pow(meas(G),2),u.col(k)).value();
    }

    void an_deriv_into(const gsMatrix<T> & u, gsMatrix<T> & result) const
    {
        result.resize(2,u.cols());
        auto G = m_evaluator.getMap(m_mp);
        for (index_t k=0; k!=u.cols(); ++k) {
          result.col(k) = m_evaluator.eval( dEdSigma(G),u.col(k) );
        }
    }

    short_t domainDim() const override
    {
        return m_mp.domainDim();
    }

    short_t targetDim() const override
    {
        return 1;
    }

protected:
    const gsMultiPatch<T> & m_mp; ///< The multi-patch to evaluate
    const gsMultiBasis<T> m_mb;
    mutable gsExprEvaluator<T> m_evaluator;
};


int main(int argc, char *argv[])
{
    // Construct the composition
    gsGeometry<>::uPtr composition = gsNurbsCreator<>::BSplineSquareDeg(2);
    // composition->coefs().row(4).setConstant(0.9);
    gsInfo<<"The composition has basis:\n"<<composition->basis()<<"\n";
    gsInfo<<"The control points of the map are:\n"<<composition->coefs()<<"\n";

    gsSquareDomain<real_t> domain(*composition);

    // Plot the composition
    gsWriteParaview(*composition,"composition");

    gsDebugVar(composition->coefs());
    // Construct a basis
//    gsKnotVector<> kv1({0,0,0,1./3.,2./3.,1,1,1},2);
//    gsKnotVector<> kv2({0,0,0,0.25,0.50,0.75,1,1,1},2);
    gsKnotVector<> kv1({0,0,0,1,1,1}, 2);
    gsKnotVector<> kv2({0,0,0,0.50,1,1,1}, 2);
    gsTensorBSplineBasis<2> tbasis(kv1, kv2);

    // Composte it
    gsComposedBasis<> cbasis(*composition, tbasis);

    // Plot the basis
    gsWriteParaview(cbasis, "basis", 1000);
    // .. and its mesh
    gsMesh<> mesh(tbasis);
    cbasis.mapMesh(mesh);
    gsWriteParaview(mesh, "mesh");

    // Construct a random geometry
    gsMatrix<> coefs(cbasis.size(),3);
    coefs.leftCols(2) = cbasis.anchors().transpose();
    for (index_t i=0; i<coefs.rows(); ++i)
    {
        coefs(i,2) = math::sin(coefs(i,0)*M_PI) * math::cos(coefs(i,1)*M_PI);
    }
//    coefs.col(2).setRandom();
    // coefs.col(2).setZero();
    gsDebugVar(coefs);

    // Make the geometries (composed and non-composed)
    gsGeometry<>::uPtr geom = tbasis.makeGeometry(coefs);
    gsGeometry<>::uPtr cgeom= cbasis.makeGeometry(coefs);

    // Plot the geometries (composed and non-composed)
    gsWriteParaview(*geom,"geom",1000,true);
    gsWriteParaview(*cgeom,"cgeom",1000,true);


    gsOptMesh<real_t,MonitorMode::ValueBased> opt(domain,*geom,&tbasis);
    gsVector<real_t> controls = gsVector<real_t>::LinSpaced(domain.nControls(),0.,1.);

    gsInfo<<"Objective function for controls: "<<controls.transpose()<<"\n";
    gsInfo<<opt.evalObj(gsAsConstVector<real_t>(controls.data(), controls.size()))<<"\n";

    gsVector<real_t> grad(controls.size());
    gsAsVector<real_t> asgrad(grad.data(),grad.rows());
    opt.gradObj_into(gsAsConstVector<real_t>(controls.data(), controls.size()),asgrad);

    gsInfo<<"Gradient of the objective function for controls: "<<controls.transpose()<<"\n";
    gsInfo<<asgrad.transpose()<<"\n";
    gsDebugVar(asgrad);

    // numerical gradient test
    const real_t epsilon = 1e-6;
    gsVector<real_t> grad_fd(controls.size());

    for (index_t i = 0; i < controls.size(); ++i)
    {
      gsVector<real_t> controls_plus = controls;
      gsVector<real_t> controls_minus = controls;

      controls_plus(i) += epsilon;
      controls_minus(i) -= epsilon;

      real_t f = opt.evalObj(gsAsConstVector<real_t>(controls.data(), controls.size()));
      real_t f_plus = opt.evalObj(gsAsConstVector<real_t>(controls_plus.data(), controls_plus.size()));
      real_t f_minus = opt.evalObj(gsAsConstVector<real_t>(controls_minus.data(), controls_minus.size()));

      grad_fd(i) = (f_plus - f_minus) / (2 * epsilon);
//      grad_fd(i) = (f_plus - f) / (epsilon);
    }

    gsInfo << "Analytical gradient:\n" << grad.transpose() << "\n";
    gsInfo << "Finite difference gradient:\n" << grad_fd.transpose() << "\n";
    gsInfo << "Difference:\n" << (grad - grad_fd).transpose() << "\n";

    // // Get the control derivative
    // gsVector<> pt(2);
    // pt.setConstant(0.5);
    // gsMatrix<real_t> dSigmadAlpha;
    // domain.control_deriv_into(pt, dSigmadAlpha);
    // gsDebug<< "dSigma/dAlpha at "<<pt.transpose()<<" is:\n"<<dSigmadAlpha<<"\n";

    // // Create the Sigma function
    // SigmaFunction<real_t> sigma(domain,pt);
    // gsInfo<<"Sigma function: R^"<<sigma.domainDim()<<" -> R^"<<sigma.targetDim()<<"\n";
    // gsInfo<<"Controls: "<<controls.transpose()<<"\n";
    // gsInfo<<"Sigma function at "<<pt.transpose()<<" is:\n";
    // gsMatrix<real_t> testPt(2,1);
    // testPt << 0.5, 0.5;
    // gsDebugVar(testPt);
    // gsMatrix<real_t> result;
    // sigma.eval_into(testPt, result);
    // gsInfo<<result.transpose()<<"\n";

    // gsInfo<<"Its derivative is:\n";
    // sigma.deriv_into(testPt, result);
    // gsInfo<<result.transpose()<<"\n";

    // gsInfo<<"Its analytical derivative is:\n";
    // sigma.an_deriv_into(testPt, result);
    // gsInfo<<result.transpose()<<"\n";



    // gsMultiPatch<> mp(*geom);
    // Efunction<real_t> Efun(mp);
    // gsInfo<<"E function: R^"<<Efun.domainDim()<<" -> R^"<<Efun.targetDim()<<"\n";
    // gsInfo<<"E function at "<<pt.transpose()<<" is:\n";
    // Efun.eval_into(pt,result);
    // gsInfo<<result<<"\n";

    // gsInfo<<"Its derivative is:\n";
    // Efun.deriv_into(pt,result);
    // gsInfo<<result<<"\n";

    // gsInfo<<"Its analytical derivative is:\n";
    // Efun.an_deriv_into(pt,result);
    // gsInfo<<" Its analytical derivative is = " << result<<"\n";

    ///////// derivative of objective function
    // auto domIt =


    return EXIT_SUCCESS;

}// end main
