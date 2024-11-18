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

    mutable gsMatrix<Scalar> theta, dtheta, bGrads, res;
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
        // NEEDED FOR ARTIFICIAL ADDITION OF EXTRA ZERO FOR 2D
        res.resize(_u.cardinality(),cols()); // cols==3

        const index_t A = _u.cardinality()/_u.dim(); // _u.data().actives.rows()
        // The parametric coordinates of the surface are the evaluation of the curve at point k
        theta = _C.data().values[0].col(k);
        dtheta = _C.data().values[1].reshapeCol(k,1,_C.data().dim.second).transpose(); // VERIFY

        // Take the gradient of the basis functions
        bGrads = _u.data().values[1].col(k);
        for (index_t d = 0; d!= cols(); ++d) // for all basis function components
        {
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

    index_t cols() const { return _u.dim(); }

    index_t cardinality_impl() const { return _u.cardinality_impl(); }

    void parse(gsExprHelper<Scalar> & evList) const
    {
        evList.add(_u);
        _u.data().flags |= NEED_ACTIVE | NEED_GRAD; // need actives for cardinality

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

    index_t rows() const { return _C.source().targetDim(); }

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

    index_t rows() const { return _C.source().targetDim(); }

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


template<class T> EIGEN_STRONG_INLINE
curve_tangent_expr<T> ctv(const gsGeometryMap<T> & S,
                          const gsGeometryMap<T> & C) { return curve_tangent_expr<T>(S,C); }

template<class T> EIGEN_STRONG_INLINE
curve_normal_expr<T>  cnv(const gsGeometryMap<T> & S,
                          const gsGeometryMap<T> & C) { return curve_normal_expr<T>(S,C); }

template<class T> EIGEN_STRONG_INLINE
curve_binormal_expr<T>  cbv(const gsGeometryMap<T> & S,
                            const gsGeometryMap<T> & C) { return curve_binormal_expr<T>(S,C); }


/*
Next:
- Binormal = n x t
- Variation of the normal vector
- Variation of the tangent vector
- Variation of the binormal vector
 */


}
}




using namespace gismo;

template <class T>
T curveCoordinate(const gsGeometry<T> & geometry, const T & value, const short_t & dir);

// template <class T>
// T findRoot();

int main(int argc, char *argv[])
{

    gsInfo<<"To do list:\n";
    gsInfo<<"- Parse expressions as compositions, i.e. pre-compute the coordinates of the curve and pre-compute the surface at this point.\n";
    gsInfo<<"- Parse dependent expressions (e.g., the use of curve_normal_expr in curve_tangent_expr in the curve_binormal_expr)\n";
    gsInfo<<"- Make the expressions work for 2D and 3D surfaces using enable_if\n";

    // bool plot = false; // If set to true, paraview file is generated and launched on exit
    // bool trim = false; // If set to true, trim/merge operations are displayed
    // bool intersect = false; // If set to true, intersection example is displayed

    gsCmdLine cmd("TODO");
    // cmd.addSwitch("plot", "Plot result in ParaView format", plot);
    // cmd.addSwitch("trim", "Basic trim/merge operations", trim);
    // cmd.addSwitch("intersect", "Intersection operations", intersect);

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
    gsInfo<<"The surface (R^"<<mp_surf.patch(0).domainDim()<<" -> R^"<<mp_surf.patch(0).targetDim()<<") is:\n"<<mp_surf.patch(0)<<"\n";
    // 2D->2D
    gsMultiPatch<> mp_curve;
    mp_curve.addPatch(curve);
    gsInfo<<"The curve (R^"<<mp_curve.patch(0).domainDim()<<" -> R^"<<mp_curve.patch(0).targetDim()<<") is:\n"<<mp_curve.patch(0)<<"\n";

    gsMultiBasis<> basis_curve(mp_curve);

    gsExprEvaluator<> ev;
    ev.setIntegrationElements(basis_curve); // ONLY NEEDED WHEN CALLING ev.integral
    auto G_curve = ev.getMap(mp_curve);
    auto G_surf = ev.getMap(mp_surf);

    gsVector<> gamma(1);
    gamma << 0.5;
    gsInfo<<"Values at point "<<gamma<<"\n";
    gsInfo << "theta           = "<<ev.eval( G_curve, gamma ).transpose() << "\n";
    gsInfo << "normal vector   = "<<ev.eval( gismo::expr::cnv(G_surf, G_curve), gamma ).transpose() << "\n";
    gsInfo << "tangent vector  = "<<ev.eval( gismo::expr::ctv(G_surf, G_curve), gamma ).transpose() << "\n";
    gsInfo << "binormal vector = "<<ev.eval( gismo::expr::cbv(G_surf, G_curve), gamma ).transpose() << "\n";

    gsDebugVar(ev.eval(jac(u),point));

    // gsMatrix<> normal = ev.eval( sn(G_surf) );
    // gsMatrix<> tangent = ev.eval( tv(G_curve) );
    // gsMatrix<> binormal = normal.cross(tangent);

    // ev.integral()

    // gsInfo<<curveCoordinate(coefs);

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