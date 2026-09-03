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
   composition->coef(4,0) = 0.1;
   composition->coef(4,1) = 0.1;
    gsDebugVar(composition->coefs());
    // Construct a basis
//    gsKnotVector<> kv1({0,0,0,1./3.,2./3.,1,1,1},2);
//    gsKnotVector<> kv2({0,0,0,0.25,0.50,0.75,1,1,1},2);
    gsKnotVector<> kv1({0,0,0,1,1,1}, 2);
//    gsKnotVector<> kv2({0,0,0,1,1,1}, 2);
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

    // Define point
    gsMatrix<> pt(2,1);
    pt.setConstant(0.5);

    // Basic data
    gsInfo<<"-------------- Basic data --------------\n";

    // Make the geometries (composed and non-composed)
    gsGeometry<>::uPtr geom = tbasis.makeGeometry(coefs);
    gsGeometry<>::uPtr cgeom= cbasis.makeGeometry(coefs);

    // Plot the geometries (composed and non-composed)
    gsWriteParaview( *geom,  "geom", 1000, true);
    gsWriteParaview(*cgeom, "cgeom", 1000, true);

    gsGeometry<>::uPtr σ = composition->clone();

    // Derivatives
    gsInfo<<"-------------- Derivatives --------------\n";
    gsMatrix<> dcG, dG, dσ;
    // dGdσ
    geom->deriv_into(pt,dG);
    dG.resize(geom->domainDim(), geom->targetDim());
    // dσ/dξ
    σ->deriv_into(pt,dσ);
    dσ.resize(σ->domainDim(), σ->targetDim());
    // exact:
    // cgeom->deriv_into(pt,dGdσ);
    // dGdσ.resize(geom->domainDim(), cgeom->targetDim());

    dcG = dσ*dG;

    gsInfo<<"Derivatives at point "<<pt.transpose()<<":\n";
    gsInfo<<" * dG/dσ: \n"<<dG<<"\n";
    gsInfo<<" * dσ/dξ: \n"<<dσ<<"\n";
    gsInfo<<" * dG/dξ = dσ/dξ * dG/dσ: \n"<<dcG<<"\n";

    // Measure
    gsInfo<<"-------------- Measure --------------\n";
    gsMatrix<> d2G;
    geom->deriv2_into(pt,d2G);
    std::vector<gsMatrix<>> hessians(3); // hessX, hessY, hessZ
    std::vector<gsMatrix<>> deriv2(4); // dξξ, dηη, dξη, dηξ

    gsDebugVar(d2G);

    // Hessian of G_x
    hessians[0].resize(geom->domainDim(), geom->domainDim());
    hessians[0](0,0) = d2G(0,0); // d2Gx/dξ2
    hessians[0](1,1) = d2G(1,0); // d2Gx/dη2
    hessians[0](1,0) =             // d2Gx/dξdη
    hessians[0](0,1) = d2G(2,0); // d2Gx/dηdξ

    // Hessian of G_y
    hessians[1].resize(geom->domainDim(), geom->domainDim());
    hessians[1](0,0) = d2G(0+3,0); // d2Gy/dξ2
    hessians[1](1,1) = d2G(1+3,0); // d2Gy/dη2
    hessians[1](1,0) =               // d2Gy/dξdη
    hessians[1](0,1) = d2G(2+3,0); // d2Gy/dηdξ

    // Hessian of G_z
    hessians[2].resize(geom->domainDim(), geom->domainDim());
    hessians[2](0,0) = d2G(0+6,0); // d2Gz/dξ2
    hessians[2](1,1) = d2G(1+6,0); // d2Gz/dη2
    hessians[2](1,0) =               // d2Gz/dξdη
    hessians[2](0,1) = d2G(2+6,0); // d2Gz/dηdξ

    // ξ ξ derivative
    std::for_each(deriv2.begin(), deriv2.end(), [&geom](gsMatrix<> & m){m.resize(1,geom->targetDim());});
    deriv2[0](0,0) = d2G(0,0);
    deriv2[0](0,1) = d2G(0+3,0);
    deriv2[0](0,2) = d2G(0+6,0);
    // η η derivative
    deriv2[1](0,0) = d2G(1,0);
    deriv2[1](0,1) = d2G(1+3,0);
    deriv2[1](0,2) = d2G(1+6,0);
    // ξ η derivative
    // η ξ derivative
    deriv2[2](0,0) = deriv2[3](0,0) = d2G(2,0);
    deriv2[2](0,1) = deriv2[3](0,1) = d2G(2+3,0);
    deriv2[2](0,2) = deriv2[3](0,2) = d2G(2+6,0);

    gsVector<real_t,3> dGdξ = dcG.row(0);
    gsVector<real_t,3> dGdη = dcG.row(1);
    gsVector<real_t,3> d2Gdξξ = deriv2[0].row(0);
    gsVector<real_t,3> d2Gdηη = deriv2[1].row(0);
    gsVector<real_t,3> d2Gdξη = deriv2[2].row(0);
    gsVector<real_t,3> d2Gdηξ = deriv2[3].row(0);

    gsVector<real_t,3> fform = dGdξ.cross(dGdη);
    real_t meas = fform.norm();
    gsInfo<<"Measure at point "<<pt.transpose()<<": "<<meas<<"\n";

    gsMatrix<> dmeasdξ(2,1);
    dmeasdξ(0,0) = ( fform / meas ).dot(d2Gdξξ.cross(dGdη) + dGdξ.cross(d2Gdξη) ); // fform/||fform|| . ( d2G/dξ2   x dG/deta + dG/dξ x d2G/dξdeta  )
    dmeasdξ(1,0) = ( fform / meas ).dot(d2Gdηξ.cross(dGdη) + dGdξ.cross(d2Gdηη) ); // fform/||fform|| . ( d2G/dξ2   x dG/deta + dG/dξ x d2G/dξdeta  )
    gsInfo<<"∇measure at point "<<pt.transpose()<<": "<<dmeasdξ.transpose()<<"\n";

    gsMatrix<> dσdα;
    domain.control_deriv_into(pt, dσdα);
    gsInfo<<"dσ/dα at point "<<pt.transpose()<<": "<<dσdα.transpose()<<"\n";

    gsInfo<<"-------------- J matrix --------------\n";
    gsMatrix<> J(geom->targetDim(), geom->domainDim());
    J.col(0) = dGdξ;
    J.col(1) = dGdη;
    gsInfo<<"J = [dG/dξ, dG/dη] = \n"<<J<<"\n";

    gsMatrix<> dJdξ(geom->targetDim(), geom->domainDim());
    dJdξ.col(0) = d2Gdξξ;
    dJdξ.col(1) = d2Gdξη;
    gsInfo<<"dJ/dξ = [d²G/dξ², d²G/dξdη] = \n"<<dJdξ<<"\n";

    gsMatrix<> dJdη(geom->targetDim(), geom->domainDim());
    dJdη.col(0) = d2Gdηξ;
    dJdη.col(1) = d2Gdηη;
    gsInfo<<"dJ/dη = [d²G/dηdξ, d²G/dη²] = \n"<<dJdη<<"\n";

    return EXIT_SUCCESS;

}// end main
