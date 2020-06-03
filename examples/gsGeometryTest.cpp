/** @file geometry_example.cpp

    @brief Tutorial on gsGeometry abstract class.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): J. Speh
*/

#include <iostream>
#include <gismo.h>
#include <typeinfo>

#  define MatExprType  auto

namespace gismo{
namespace expr{

template<class T>
class otangent_expr : public _expr<otangent_expr<T> >
{
    typename gsGeometryMap<T>::Nested_t _G;

public:
    typedef T Scalar;

    otangent_expr(const gsGeometryMap<T> & G) : _G(G) { }

    mutable gsVector<T,3> onormal, normal;

    MatExprType eval(const index_t k) const
    {
        GISMO_ASSERT(_G.data().dim.second==3,"Domain dimension should be 3, is "<<_G.data().dim.second);

        gsMatrix<T> Jacobian = _G.data().jacobian(k);

        onormal = _G.data().outNormal(k);
        normal =  _G.data().normal(k);

        return normal.cross(onormal);

        // gsMatrix<T> result(_G.data().dim.second,1);
        // result.setZero();
        // return result;
        // if (_G.data().dim.second!=3)
        //     return normal.head(_G.data().dim.second);
        // else
        //     return normal;
    }

    index_t rows() const { return _G.data().dim.second; }
    index_t cols() const { return 1; }

    const gsFeSpace<T> & rowVar() const {return gsNullExpr<T>::get();}
    const gsFeSpace<T> & colVar() const {return gsNullExpr<T>::get();}

    static constexpr bool rowSpan() {return false;}
    static bool colSpan() {return false;}

    void setFlag() const { _G.data().flags |= NEED_OUTER_NORMAL | NEED_NORMAL | NEED_DERIV; }

    void parse(gsSortedVector<const gsFunctionSet<Scalar>*> & evList) const
    {
        //GISMO_ASSERT(NULL!=m_fd, "FeVariable: FuncData member not registered");
        evList.push_sorted_unique(&_G.source());
        _G.data().flags |= NEED_OUTER_NORMAL | NEED_NORMAL | NEED_DERIV;
    }

    // Normalized to unit length
    normalized_expr<otangent_expr<T> > normalized()
    { return normalized_expr<otangent_expr<T> >(*this); }

    void print(std::ostream &os) const { os << "tv2("; _G.print(os); os <<")"; }
};

template<class E>
class tvar1_expr : public _expr<tvar1_expr<E> >
{
public:
    typedef typename E::Scalar Scalar;

private:
    typename E::Nested_t _u;
    typename gsGeometryMap<Scalar>::Nested_t _G;

public:
    enum{ Space = E::Space };

    tvar1_expr(const E & u, const gsGeometryMap<Scalar> & G) : _u(u), _G(G) { }

    mutable gsVector<Scalar,3> onormal, tangent, dtan;
    mutable gsVector<Scalar> tmp;
    mutable gsMatrix<Scalar> bGrads, cJac, res;
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    // helper function
    static inline gsVector<Scalar,3> vecFun(index_t pos, Scalar val)
    {
        gsVector<Scalar,3> result = gsVector<Scalar,3>::Zero();
        result[pos] = val;
        return result;
    }

    MatExprType eval(const index_t k) const
    {
        GISMO_ASSERT(_G.data().dim.second==3,"Domain dimension should be 3, is "<<_G.data().dim.second);

        const index_t A = _u.cardinality()/_u.targetDim();
        res.resize(A*_u.targetDim(), cols()); // rows()*
        cJac = _G.data().jacobian(k);

        onormal = _G.data().outNormal(k);
        tmp = cJac.transpose() * onormal;
        Scalar tol = 1e-8;

        gsDebugVar(onormal);
        /*
            We can check which column of the Jacobian corresponds to the outer normal vector or to the tangent.
            The tangent is a covariant vector and hence the column of the Jacobian should be equal to the tangent.
            The normal is a contravariant vector and hence the corresponding column of the Jacobian times the outward normal should give 1. We use this property.
        */
        index_t colIndex;
        if ( (abs(tmp.at(0)) < tol) && (abs(tmp.at(1)) > 1-tol ) )         // then the normal is vector 2 and the tangent vector 1
            colIndex = 0;
        else if ( (abs(tmp.at(1)) < tol) && (abs(tmp.at(0)) > 1-tol ) )     // then the normal is vector 1 and the tangent vector 2
            colIndex = 1;
        else                    // then the normal is unknown??
            gsInfo<<"warning: choice unknown\n";

        tangent = cJac.col(colIndex);

        // Now we will compute the derivatives of the basis functions
        bGrads = _u.data().values[1].col(k);
        for (index_t d = 0; d!= cols(); ++d) // for all basis function components
        {
            const short_t s = d*A;
            for (index_t j = 0; j!= A; ++j) // for all actives
            {
                // The tangent vector is in column colIndex in cJac and thus in 2*j+colIndex in bGrads.
                // Furthermore, as basis function for dimension d, it has a nonzero in entry d, and zeros elsewhere
                dtan = vecFun(d, bGrads.at(2*j+colIndex));
                res.row(s+j).noalias() = (1 / tangent.norm() * ( dtan - ( tangent * dtan ) * tangent / (tangent.norm() * tangent.norm()) )).transpose();
            }
        }
        return res;
    }

    index_t rows() const { return 1; }
    index_t cols() const { return _u.dim(); }

    const gsFeSpace<Scalar> & rowVar() const { return _u.rowVar(); }
    const gsFeSpace<Scalar> & colVar() const {return gsNullExpr<Scalar>::get();}
    index_t cardinality_impl() const { return _u.cardinality_impl(); }

    static constexpr bool rowSpan() {return E::rowSpan(); }
    static constexpr bool colSpan() {return false;}

    void setFlag() const
    {
        _u.data().flags |= NEED_GRAD | NEED_ACTIVE;
        _G.data().flags |= NEED_NORMAL | NEED_DERIV | NEED_MEASURE;
    }

    void parse(gsSortedVector<const gsFunctionSet<Scalar>*> & evList) const
    {
        evList.push_sorted_unique(&_u.source());
        _u.data().flags |= NEED_GRAD;
        _G.data().flags |= NEED_NORMAL | NEED_DERIV | NEED_MEASURE;
    }

    void print(std::ostream &os) const { os << "tvar("; _G.print(os); os <<")"; }
};

template<class E>
class ovar1_expr : public _expr<ovar1_expr<E> >
{
public:
    typedef typename E::Scalar Scalar;

private:
    typename E::Nested_t _u;
    typename gsGeometryMap<Scalar>::Nested_t _G;

public:
    enum{ Space = E::Space };

    ovar1_expr(const E & u, const gsGeometryMap<Scalar> & G) : _u(u), _G(G) { }

    mutable gsVector<Scalar,3> onormal, tangent, normal, dtan, tvar, snvar, mv;
    mutable gsVector<Scalar> tmp;
    mutable gsMatrix<Scalar> bGrads, cJac, res;
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    // helper function
    static inline gsVector<Scalar,3> vecFun(index_t pos, Scalar val)
    {
        gsVector<Scalar,3> result = gsVector<Scalar,3>::Zero();
        result[pos] = val;
        return result;
    }

    MatExprType eval(const index_t k) const
    {
        GISMO_ASSERT(_G.data().dim.second==3,"Domain dimension should be 3, is "<<_G.data().dim.second);
        const index_t A = _u.cardinality()/_u.targetDim();
        res.resize(A*_u.targetDim(), cols()); // rows()*

        cJac = _G.data().jacobian(k);

        onormal = _G.data().outNormal(k);
        normal  =  _G.data().normal(k);
        tmp = cJac.transpose() * onormal;
        Scalar tol = 1e-8;

        gsDebugVar(onormal);
        /*
            We can check which column of the Jacobian corresponds to the outer normal vector or to the tangent.
            The tangent is a covariant vector and hence the column of the Jacobian should be equal to the tangent.
            The normal is a contravariant vector and hence the corresponding column of the Jacobian times the outward normal should give 1. We use this property.
        */
        index_t colIndex;
        if ( (abs(tmp.at(0)) < tol) && (abs(tmp.at(1)) > 1-tol ) )         // then the normal is vector 2 and the tangent vector 1
            colIndex = 0;
        else if ( (abs(tmp.at(1)) < tol) && (abs(tmp.at(0)) > 1-tol ) )     // then the normal is vector 1 and the tangent vector 2
            colIndex = 1;
        else                    // then the normal is unknown??
            gsInfo<<"warning: choice unknown\n";

        tangent = cJac.col(colIndex);

        // Now we will compute the derivatives of the basis functions
        bGrads = _u.data().values[1].col(k);

        const Scalar measure =  _G.data().measures.at(k);
        for (index_t d = 0; d!= cols(); ++d) // for all basis function components
        {
            const short_t s = d*A;
            for (index_t j = 0; j!= A; ++j) // for all actives
            {
                // The tangent vector is in column colIndex in cJac and thus in 2*j+colIndex in bGrads.
                // Furthermore, as basis function for dimension d, it has a nonzero in entry d, and zeros elsewhere
                dtan = vecFun(d, bGrads.at(2*j+colIndex));
                tvar.noalias() = 1 / tangent.norm() * ( dtan - ( tangent.dot(dtan) ) * tangent / (tangent.norm() * tangent.norm()) );

                // Jac(u) ~ Jac(G) with alternating signs ?..
                mv.noalias() = (vecFun(d, bGrads.at(2*j  ) ).cross( cJac.col(1).template head<3>() )
                              - vecFun(d, bGrads.at(2*j+1) ).cross( cJac.col(0).template head<3>() )) / measure;

                // ---------------  First variation of the normal
                snvar.noalias() = mv - ( normal.dot(mv) ) * normal;

                res.row(s+j).noalias() = tvar.cross(normal) + tangent.cross(snvar);
            }
        }
        return res;
    }

    index_t rows() const { return 1; }
    index_t cols() const { return _u.dim(); }

    const gsFeSpace<Scalar> & rowVar() const { return _u.rowVar(); }
    const gsFeSpace<Scalar> & colVar() const {return gsNullExpr<Scalar>::get();}
    index_t cardinality_impl() const { return _u.cardinality_impl(); }

    static constexpr bool rowSpan() {return E::rowSpan(); }
    static constexpr bool colSpan() {return false;}

    void setFlag() const
    {
        _u.data().flags |= NEED_GRAD | NEED_ACTIVE;
        _G.data().flags |= NEED_NORMAL | NEED_OUTER_NORMAL | NEED_DERIV | NEED_MEASURE;
    }

    void parse(gsSortedVector<const gsFunctionSet<Scalar>*> & evList) const
    {
        evList.push_sorted_unique(&_u.source());
        _u.data().flags |= NEED_GRAD;
        _G.data().flags |= NEED_NORMAL | NEED_OUTER_NORMAL | NEED_DERIV | NEED_MEASURE;
    }

    void print(std::ostream &os) const { os << "ovar("; _G.print(os); os <<")"; }
};

template<class E1, class E2, class E3>
class ovar2_expr : public _expr<ovar2_expr<E1,E2,E3> >
{
public:
    typedef typename E1::Scalar Scalar;

private:
    typename E1::Nested_t _u;
    typename E2::Nested_t _v;
    typename gsGeometryMap<Scalar>::Nested_t _G;
    typename E3::Nested_t _C;

public:
    enum{ Space = E1::Space };

    ovar2_expr(const E1 & u, const E2 & v, const gsGeometryMap<Scalar> & G, _expr<E3> const& C) : _u(u), _v(v), _G(G), _C(C) { }

    mutable gsVector<Scalar,3> onormal, normal, tangent, dtanu, dtanv, tvaru, tvarv, tvar2,
                            mu, mv, muv, mu_der, snvaru, snvarv, snvar2, nvar2;
    mutable gsVector<Scalar> tmp;
    mutable gsMatrix<Scalar> uGrads, vGrads, cJac, res, eC;

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    // helper function
    static inline gsVector<Scalar,3> vecFun(index_t pos, Scalar val)
    {
        gsVector<Scalar,3> result = gsVector<Scalar,3>::Zero();
        result[pos] = val;
        return result;
    }

    MatExprType eval(const index_t k) const
    {
        GISMO_ASSERT(_G.data().dim.second==3,"Domain dimension should be 3, is "<<_G.data().dim.second);
        GISMO_ASSERT(_C.cols()*_C.rows()==3, "Size of vector is incorrect");
        res.resize(_u.cardinality(), _u.cardinality());

        eC = _C.eval(k);
        eC.resize(1,3);

        cJac = _G.data().jacobian(k);

        onormal = _G.data().outNormal(k);
        normal =  _G.data().normal(k);
        tmp = cJac.transpose() * onormal;
        Scalar tol = 1e-8;

        /*
            We can check which column of the Jacobian corresponds to the outer normal vector or to the tangent.
            The tangent is a covariant vector and hence the column of the Jacobian should be equal to the tangent.
            The normal is a contravariant vector and hence the corresponding column of the Jacobian times the outward normal should give 1. We use this property.
        */
        index_t colIndex;
        if ( (abs(tmp.at(0)) < tol) && (abs(tmp.at(1)) > 1-tol ) )         // then the normal is vector 2 and the tangent vector 1
            colIndex = 0;
        else if ( (abs(tmp.at(1)) < tol) && (abs(tmp.at(0)) > 1-tol ) )     // then the normal is vector 1 and the tangent vector 2
            colIndex = 1;
        else                    // then the normal is unknown??
            gsInfo<<"warning: choice unknown\n";

        tangent = cJac.col(colIndex);

        // Now we will compute the derivatives of the basis functions
        uGrads = _u.data().values[1].col(k);
        vGrads = _v.data().values[1].col(k);
        const index_t cardU = _u.data().values[0].rows(); // number of actives per component of u
        const index_t cardV = _v.data().values[0].rows(); // number of actives per component of v

        const Scalar measure =  _G.data().measures.at(k);

        for (index_t j = 0; j!= cardU; ++j) // for all basis functions u (1)
        {
            for (index_t i = 0; i!= cardV; ++i) // for all basis functions v (1)
            {
                for (index_t d = 0; d!= _u.dim(); ++d) // for all basis functions u (2)
                {
                    const short_t s = d*cardU;

                    // first variation of the tangent (colvector)
                    dtanu = vecFun(d, uGrads.at(2*j+colIndex));
                    tvaru = 1 / tangent.norm() * ( dtanu - ( tangent.dot(dtanu) ) * tangent / (tangent.norm() * tangent.norm()) );

                    // first variation of the surface normal (colvector)
                    mu.noalias() = ( vecFun(d, uGrads.at(2*j  ) ).cross( cJac.col(1).template head<3>() )
                                    -vecFun(d, uGrads.at(2*j+1) ).cross( cJac.col(0).template head<3>() ))
                                    / measure;
                    snvaru.noalias() = (mu - ( normal.dot(mu) ) * normal);

                    for (index_t c = 0; c!= _v.dim(); ++c) // for all basis functions v (2)
                    {
                        const short_t r = c*cardV;

                        // first variation of the tangent (colvector)
                        dtanv = vecFun(c, vGrads.at(2*i+colIndex));
                        tvarv = 1 / tangent.norm() * ( dtanv - ( tangent.dot(dtanv) ) * tangent / (tangent.norm() * tangent.norm()) );

                        // first variation of the surface normal (colvector)
                        mv.noalias() = ( vecFun(c, vGrads.at(2*i  ) ).cross( cJac.col(1).template head<3>() )
                                        -vecFun(c, vGrads.at(2*i+1) ).cross( cJac.col(0).template head<3>() ))
                                        / measure;
                        snvarv.noalias() = (mv - ( normal.dot(mv) ) * normal);


                        // Second variation of the tangent (colvector)
                        tvar2 = 1 / tangent.norm() * ( tvarv.dot(dtanu) * tangent )
                                + 1 / (tangent.norm()*tangent.norm())
                                * ( 2*( (tangent.dot(dtanu))*(tangent.dot(dtanv))*tangent )
                                    - ( tangent.dot(dtanu)*dtanv ) - ( tangent.dot(dtanv)*dtanu) );

                        // Second variation of the surface normal (colvector)
                        muv.noalias() = ( vecFun(d, uGrads.at(2*j  ) ).cross( vecFun(c, vGrads.at(2*i+1) ) )
                                         +vecFun(c, vGrads.at(2*i  ) ).cross( vecFun(d, uGrads.at(2*j+1) ) ))
                                          / measure;

                        mu_der.noalias() = (muv - ( normal.dot(mv) ) * mu);

                        snvar2 = mu_der - (mu.dot(snvarv) + normal.dot(mu_der) ) * normal - (normal.dot(mu) ) * snvarv;

                        // Second variation of the outer normal (colvector)
                        nvar2 = tvar2.cross(normal) + tvaru.cross(snvarv) + tvarv.cross(snvaru) + tangent.cross(snvar2);

                        res(s + j, r + i ) = (eC * nvar2)(0,0);
                    }
                }
            }
        }
        return res;
    }

    index_t cardinality_impl() const { return 1; }

    index_t rows() const { return 1; }

    index_t cols() const { return 1; }

    void setFlag() const
    {
        _u.data().flags |= NEED_GRAD;
        _v.data().flags |= NEED_GRAD;
        _G.data().flags |= NEED_NORMAL | NEED_OUTER_NORMAL | NEED_DERIV | NEED_MEASURE;
        _C.setFlag();
    }

    void parse(gsSortedVector<const gsFunctionSet<Scalar>*> & evList) const
    {
        evList.push_sorted_unique(&_u.source());
        _u.data().flags |= NEED_GRAD;
        _v.data().flags |= NEED_GRAD;
        _G.data().flags |= NEED_NORMAL | NEED_OUTER_NORMAL | NEED_DERIV | NEED_MEASURE;
        _C.setFlag();
    }

    const gsFeSpace<Scalar> & rowVar() const { return _u.rowVar(); }
    const gsFeSpace<Scalar> & colVar() const { return _v.rowVar(); }

    static constexpr bool rowSpan() {return true; }
    static constexpr bool colSpan() {return true; }

    void print(std::ostream &os) const { os << "nvar2("; _u.print(os); os <<")"; }
};


template<class T> EIGEN_STRONG_INLINE
otangent_expr<T> otangent(const gsGeometryMap<T> & u) { return otangent_expr<T>(u); }

template<class E> EIGEN_STRONG_INLINE
tvar1_expr<E> tvar1(const E & u, const gsGeometryMap<typename E::Scalar> & G) { return tvar1_expr<E>(u, G); }

template<class E> EIGEN_STRONG_INLINE
ovar1_expr<E> ovar1(const E & u, const gsGeometryMap<typename E::Scalar> & G) { return ovar1_expr<E>(u, G); }

template<class E1, class E2, class E3> EIGEN_STRONG_INLINE
ovar2_expr<E1,E2,E3> ovar2(const E1 & u, const E2 & v, const gsGeometryMap<typename E1::Scalar> & G, const E3 & C)
{ return ovar2_expr<E1,E2,E3>(u,v, G, C); }
}
}

using namespace gismo;


// Returns the string with the size of a matrix.
template <typename T>
std::string size(const gsMatrix<T>& matrix)
{
    std::string result = "(" + util::to_string(matrix.rows()) + " x " +
        util::to_string(matrix.cols()) + ")";

    return result;
}

template <class T>
void evaluateFunction(gsExprEvaluator<T> ev, auto expression, gsVector<T> pt);


int main(int argc, char* argv[])
{
    bool plot = false;
    index_t numRefine  = 1;
    index_t numElevate = 1;
    std::string input("surfaces/simple.xml");
    std::string output("");

    gsCmdLine cmd("Tutorial on gsGeometry class.");
    cmd.addPlainString("filename", "G+Smo input geometry file.", input);
    cmd.addString("o", "output", "Name of the output file", output);
    cmd.addInt( "e", "degreeElevation",
                "Number of degree elevation steps to perform before solving (0: equalize degree in all directions)", numElevate );
    cmd.addInt( "r", "uniformRefine", "Number of Uniform h-refinement steps to perform before solving",  numRefine );
    cmd.addSwitch("plot", "Create a ParaView visualization file with the solution", plot);
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    // ======================================================================
    // reading the geometry
    // ======================================================================

    gsMultiPatch<> mp;
    gsReadFile<>(input, mp);

    // p-refine
    if (numElevate!=0)
        mp.degreeElevate(numElevate);

    // h-refine
    for (int r =0; r < numRefine; ++r)
        mp.uniformRefine();

    if (plot)
        gsWriteParaview<>( mp, "mp", 1000, true);

    if (mp.domainDim()!=3)
        mp.embed(3);

    gsMultiBasis<> dbasis(mp);

    // ======================================================================
    // printing some information about the basis
    // ======================================================================


    // ----------------------------------------------------------------------
    // printing the geometry
    // ----------------------------------------------------------------------

    //! [printing the geometry]
    gsInfo << "The file contains: \n" << mp << "\n";

    // G+Smo geometries contains basis and coefficients
    gsInfo << "\nBasis 0: \n" << dbasis << "\n";

    // const gsMatrix<>& coefs = mp.patch(0).coefs();
    // gsInfo << "\nCoefficients 0: \n" << coefs << "\n" << "\n";
    //! [printing the geometry]


    // ======================================================================
    // test normals etc
    // ======================================================================
    typedef gsExprAssembler<>::geometryMap geometryMap;
    typedef gsExprAssembler<>::space       space;

    gsExprAssembler<> A;
    A.setIntegrationElements(dbasis);
    gsExprEvaluator<> ev(A);
    geometryMap G = A.getMap(mp);

    space u = A.getSpace(dbasis, 3);

    gsVector<> pt(2);
    pt << 0,0.5;

    auto onor   = nv(G);
    auto otan   = otangent(G);
    auto tanv   = tv(G);
    auto snor   = sn(G);
    auto tvar   = tvar1(u,G);
    auto ovar   = ovar1(u,G);
    auto ovar2f = ovar2(u,u,G,sn(G));
    evaluateFunction(ev, G, pt); // evaluates an expression on a point
    evaluateFunction(ev, onor, pt); // evaluates an expression on a point
    evaluateFunction(ev, snor, pt); // evaluates an expression on a point
    evaluateFunction(ev, otan, pt); // evaluates an expression on a point
    evaluateFunction(ev, jac(G), pt); // evaluates an expression on a point
    // evaluateFunction(ev, tanv, pt); // evaluates an expression on a point
    // evaluateFunction(ev, tanv.tr()*onor, pt); // evaluates an expression on a point
    evaluateFunction(ev, otan.tr()*onor, pt); // evaluates an expression on a point
    evaluateFunction(ev, snor.tr()*onor, pt); // evaluates an expression on a point
    evaluateFunction(ev, snor.tr()*otan, pt); // evaluates an expression on a point

    evaluateFunction(ev, ovar2f, pt); // evaluates an expression on a point
    return 0;
}

template <class T>
void evaluateFunction(gsExprEvaluator<T> ev, auto expression, gsVector<T> pt)
{
    gsMatrix<T> evresult = ev.eval( expression,pt );
    gsInfo << "Eval on point ("<<pt.at(0)<<" , "<<pt.at(1)<<") :\n"<< evresult;
    gsInfo << "\nEnd ("<< evresult.rows()<< " x "<<evresult.cols()<<")\n";
};


