/** @file flow-over-heated-plate.cpp

    @brief Heat equation participant for the PreCICE example "flow over heated plate"

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst
*/

#include <gismo.h>

namespace gismo
{
namespace expr
{
    /*
  Expression for the gradient of a finite element variable

  Transposed gradient vectors are returned as a matrix
*/
template<class E>
class subgrad_expr : public _expr<subgrad_expr<E> >
{
    typename E::Nested_t _u;
public:
    enum {Space = E::Space, ScalarValued= 0, ColBlocks= 0};

    typedef typename E::Scalar Scalar;
    mutable gsMatrix<Scalar> tmp;
    short_t _start, _end;

    subgrad_expr(const E & u, short_t startOffset = 0, short_t endOffset = 0)
    :
    _u(u),
    _start(startOffset),
    _end(endOffset)
    {
        GISMO_ASSERT(1==u.dim(),"grad(.) requires 1D variable, use jac(.) instead.");
        GISMO_ASSERT(0<=startOffset && 0<=endOffset,"grad(.) requires non-negative offsets.");
        GISMO_ASSERT(startOffset+endOffset<=u.source().domainDim(),"grad(.) requires startOffset+endOffset<=u.source().domainDim().");
        GISMO_ASSERT(startOffset<u.source().domainDim() && endOffset<u.source().domainDim(),"grad(.) requires startOffset<endOffset<u.sourceDomainDim().");
        // GISMO_ASSERT(startOffset<u.source().domainDim()-
    }

    const gsMatrix<Scalar> & eval(const index_t k) const
    {
        // Assumes: derivatives are in _u.data().values[1]
        // gsExprHelper acounts for compositions/physical expressions
        // so that derivs are directly computed
        tmp = _u.data().values[1].reshapeCol(k, _u.source().domainDim(), cardinality_impl()).transpose()(gsEigen::all, gsEigen::seq(_start,gsEigen::last-_end));
        return tmp;
    }

    index_t rows() const { return 1 /*==u.dim()*/; }

    index_t cols() const { return _u.source().domainDim()-_start-_end; }

    index_t cardinality_impl() const
    { return _u.data().values[1].rows() /  _u.source().domainDim(); }

    void parse(gsExprHelper<Scalar> & evList) const
    {
        evList.add(_u);
        _u.data().flags |= NEED_GRAD;
    }

    const gsFeSpace<Scalar> & rowVar() const { return _u.rowVar(); }
    const gsFeSpace<Scalar> & colVar() const
    {return gsNullExpr<Scalar>::get();}

    void print(std::ostream &os) const { os << "\u2207("; _u.print(os); os <<")"; }
private:

    template<class U> static inline
    typename util::enable_if<util::is_same<U,gsComposition<Scalar> >::value,
                             const gsMatrix<Scalar> &>::type
    eval_impl(const U & u, const index_t k)
    {
        return u.eval(k);
    }
};

/*
  \brief Expression for the gradient of a finite element variable

  Transposed gradient vectors are returned as a matrix.
  This specialization is for a gsFeSolution object
*/

template<class T>
class subgrad_expr<gsFeSolution<T> > : public _expr<subgrad_expr<gsFeSolution<T> > >
{
protected:
    const gsFeSolution<T> _u;

public:
    typedef T Scalar;
    enum {Space = 0, ScalarValued= 0, ColBlocks= 0};
    short_t _start, _end;

    explicit subgrad_expr(const gsFeSolution<T> & u, short_t startOffset = 0, short_t endOffset = 0)
    :
    _u(u),
    _start(startOffset),
    _end(endOffset)
    {
        GISMO_ASSERT(1==u.dim(),"grad(.) requires 1D variable, use jac(.) instead.");
        GISMO_ASSERT(0<=startOffset && 0<=endOffset,"grad(.) requires non-negative offsets.");
        GISMO_ASSERT(startOffset+endOffset<=u.parDim(),"grad(.) requires startOffset+endOffset<=u.parDim().");
        GISMO_ASSERT(startOffset<u.parDim() && endOffset<u.parDim(),"grad(.) requires startOffset<endOffset<u.sourceDomainDim().");
    }

    mutable gsMatrix<T> res;
    mutable gsMatrix<T> der;
    const gsMatrix<T> & eval(index_t k) const
    {
        GISMO_ASSERT(1==_u.data().actives.cols(), "Single actives expected");

        res.setZero(_u.dim(), _u.parDim()-_start-_end);
        const gsDofMapper & map = _u.mapper();
        for (index_t c = 0; c!= _u.dim(); c++)
        {
            for (index_t i = 0; i!=_u.data().actives.size(); ++i)
            {
                const index_t ii = map.index(_u.data().actives.at(i), _u.data().patchId,c);
                der = _u.data().values[1].col(k)
                            .segment(i*_u.parDim(), _u.parDim()).transpose()
                                        ( 0,gsEigen::seq(_start,gsEigen::last-_end) );
                if ( map.is_free_index(ii) ) // DoF value is in the solVector
                    res.row(c) += _u.coefs().at(ii) * der;
                else
                    res.row(c) += _u.fixedPart().at( map.global_to_bindex(ii) ) * der;
            }
        }
        return res;
    }

    index_t rows() const {return _u.dim();}
    index_t cols() const {return _u.parDim()-_start-_end; }

    const gsFeSpace<Scalar> & rowVar() const
    {return gsNullExpr<Scalar>::get();}
    const gsFeSpace<Scalar> & colVar() const
    {return gsNullExpr<Scalar>::get();}

    void parse(gsExprHelper<Scalar> & evList) const
    {
        _u.parse(evList);                         // add symbol
        evList.add(_u.space());
        _u.data().flags |= NEED_GRAD|NEED_ACTIVE; // define flags
    }

    void print(std::ostream &os) const { os << "\u2207(s)"; }
};


/*
  Expression for the Jacobian matrix of a geometry map
*/
template<class T>
class subjac_expr : public _expr<subjac_expr<T> >
{
    typename gsGeometryMap<T>::Nested_t _G;

public:
    typedef T Scalar;
    enum {Space = 0, ScalarValued= 0, ColBlocks= 0};
    short_t _start, _end;

    subjac_expr(const gsGeometryMap<T> & G, short_t startOffset, short_t endOffset)
    :
    _G(G),
    _start(startOffset),
    _end(endOffset)
    {
        GISMO_ASSERT(0<=startOffset && 0<=endOffset,"grad(.) requires non-negative offsets.");
        GISMO_ASSERT(startOffset+endOffset<=_G.source().domainDim(),"grad(.) requires startOffset+endOffset<=_G.source().domainDim().");
        GISMO_ASSERT(startOffset<_G.source().domainDim() && endOffset<_G.source().domainDim(),"grad(.) requires startOffset<endOffset<u.sourceDomainDim().");
        // GISMO_ASSERT(startOffset<_G.source().domainDim()-
    }

    auto eval(const index_t k) const
    {
        // TarDim x ParDim
        return _G.data().values[1]
            .reshapeCol(k, _G.data().dim.first, _G.data().dim.second).transpose()(gsEigen::all, gsEigen::seq(_start,gsEigen::last-_end));
    }

    index_t rows() const { return _G.source().targetDim(); }

    index_t cols() const { return _G.source().domainDim()-_start-_end; }

    static const gsFeSpace<Scalar> & rowVar() { return gsNullExpr<Scalar>::get(); }
    static const gsFeSpace<Scalar> & colVar() { return gsNullExpr<Scalar>::get(); }

    void parse(gsExprHelper<Scalar> & evList) const
    {
        evList.add(_G);
        _G.data().flags |= NEED_DERIV;
    }

    meas_expr<T> absDet() const
    {
        GISMO_ASSERT(rows() == cols(), "The Jacobian matrix is not square");
        return meas_expr<T>(_G);
    }

    jacInv_expr<T> inv() const
    {
        GISMO_ASSERT(rows() == cols(), "The Jacobian matrix is not square");
        return jacInv_expr<T>(_G);
    }

    /// The generalized Jacobian matrix inverse, i.e.: (J^t J)^{-t} J^t
    jacInv_expr<T> ginv() const { return jacInv_expr<T>(_G); }

    void print(std::ostream &os) const { os << "\u2207("; _G.print(os); os <<")"; }
};

/**
   Expression for the Laplacian of a finite element variable
*/
template<class E>
class sublapl_expr : public _expr<sublapl_expr<E> >
{
    typename E::Nested_t _u;

public:
    typedef typename E::Scalar Scalar;
    enum {Space = E::Space, ScalarValued= 0, ColBlocks= 0};
    short_t _start, _end;
    mutable gsMatrix<Scalar> der2, hess, result;

    sublapl_expr(const E & u, short_t startOffset, short_t endOffset)
    :
    _u(u),
    _start(startOffset),
    _end(endOffset)
    {
        GISMO_ASSERT(1==u.dim(),"lapl(.) requires 1D variable, use lapl(.) instead.");
        GISMO_ASSERT(0<=startOffset && 0<=endOffset,"lapl(.) requires non-negative offsets.");
        GISMO_ASSERT(startOffset+endOffset<=u.source().domainDim(),"lapl(.) requires startOffset+endOffset<=u.source().domainDim().");
        GISMO_ASSERT(startOffset<u.source().domainDim() && endOffset<u.source().domainDim(),"lapl(.) requires startOffset<endOffset<u.sourceDomainDim().");
    }

    const gsMatrix<Scalar> eval(const index_t k) const
    {
        index_t numActs = _u.data().values[0].rows();
        index_t numDers = _u.parDim() * (_u.parDim() + 1) / 2;
        result.resize(numActs,1);
        for (index_t i = 0; i!=numActs; ++i)
        {
            der2 = _u.data().values[2].block(i*numDers,k,_u.parDim(),1); // this only takes d11, d22, d33 part. For all the derivatives [d11, d22, d33, d12, d13, d23]: col.block(i*numDers,k,numDers,1)
            result(i,0) = der2(gsEigen::seq(_start,gsEigen::last-_end),0).sum();
        }
        return result;
    }

    index_t rows() const { return 1; }
    index_t cols() const { return 1; }

    index_t cardinality_impl() const { return _u.cardinality_impl(); }

    void parse(gsExprHelper<Scalar> & evList) const
    {
        evList.add(_u);
        _u.data().flags |= NEED_DERIV2;
    }

    const gsFeSpace<Scalar> & rowVar() const { return _u.rowVar(); }
    const gsFeSpace<Scalar> & colVar() const { return gsNullExpr<Scalar>::get(); }

    void print(std::ostream &os) const { os << "\u2206("; _u.print(os); os <<")"; } //or \u0394
};

/*
  Expression for the Laplacian of a finite element solution
*/
template<class T>
class sublapl_expr<gsFeSolution<T> > : public _expr<sublapl_expr<gsFeSolution<T> > >
{
protected:
    const gsFeSolution<T> _u;

public:
    typedef T Scalar;
    enum {Space = 0, ScalarValued= 0, ColBlocks= 0};
    short_t _start, _end;
    mutable gsMatrix<Scalar> der2, hess, result;

    sublapl_expr(const gsFeSolution<T> & u, short_t startOffset, short_t endOffset)
    :
    _u(u),
    _start(startOffset),
    _end(endOffset)
    {
        GISMO_ASSERT(1==u.dim(),"lapl(.) requires 1D variable, use lapl(.) instead.");
        GISMO_ASSERT(0<=startOffset && 0<=endOffset,"lapl(.) requires non-negative offsets.");
        GISMO_ASSERT(startOffset+endOffset<=u.parDim(),"lapl(.) requires startOffset+endOffset<=u.parDim().");
        GISMO_ASSERT(startOffset<u.parDim() && endOffset<u.parDim(),"lapl(.) requires startOffset<endOffset<u.sourceDomainDim().");
    }

    mutable gsMatrix<T> res;
    const gsMatrix<T> & eval(const index_t k) const
    {
        GISMO_ASSERT(1==_u.data().actives.cols(), "Single actives expected");

        res.setZero(_u.dim(), 1); //  scalar, but per component
        const gsDofMapper & map = _u.mapper();

        index_t numActs = _u.data().values[0].rows();
        index_t numDers = _u.parDim() * (_u.parDim() + 1) / 2;
        gsMatrix<T> deriv2;

        for (index_t c = 0; c!= _u.dim(); c++)
            for (index_t i = 0; i!=numActs; ++i)
            {
                const index_t ii = map.index(_u.data().actives.at(i), _u.data().patchId,c);
                deriv2 = _u.data().values[2].block(i*numDers,k,_u.parDim(),1); // this only takes d11, d22, d33 part. For all the derivatives [d11, d22, d33, d12, d13, d23]: col.block(i*numDers,k,numDers,1)
                deriv2 = deriv2(gsEigen::seq(_start,gsEigen::last-_end),0);
                if ( map.is_free_index(ii) ) // DoF value is in the solVector
                    res.at(c) += _u.coefs().at(ii) * deriv2.sum();
                else
                    res.at(c) +=_u.fixedPart().at( map.global_to_bindex(ii) ) * deriv2.sum();
            }
        return res;
    }

    index_t rows() const { return _u.dim(); }
    index_t cols() const { return 1; }

    void parse(gsExprHelper<Scalar> & evList) const
    {
        evList.add(_u.space());
        _u.data().flags |= NEED_ACTIVE | NEED_DERIV2;
    }

    const gsFeSpace<Scalar> & rowVar() const { return gsNullExpr<Scalar>::get(); }
    const gsFeSpace<Scalar> & colVar() const { return gsNullExpr<Scalar>::get(); }

    void print(std::ostream &os) const { os << "\u2206(s)"; }
};

/// The gradient of a variable
template<class E> EIGEN_STRONG_INLINE
subgrad_expr<E> subgrad(const E & u, short_t startOffset, short_t endOffset)
{ return subgrad_expr<E>(u,startOffset,endOffset); }

/// The gradient of a finite element solution
template<class T> EIGEN_STRONG_INLINE
subgrad_expr<gsFeSolution<T> > subgrad(const gsFeSolution<T> & u, short_t startOffset, short_t endOffset)
{ return subgrad_expr<gsFeSolution<T> >(u,startOffset,endOffset); }

/// The Jacobian matrix of a geometry map
template<class T> EIGEN_STRONG_INLINE
subjac_expr<T> subjac(const gsGeometryMap<T> & G, short_t startOffset, short_t endOffset)
{return subjac_expr<T>(G,startOffset,endOffset);}

template<class E> EIGEN_STRONG_INLINE
sublapl_expr<E> sublapl(const symbol_expr<E> & u, short_t startOffset, short_t endOffset)
{ return sublapl_expr<E>(u,startOffset,endOffset); }

template<class T> EIGEN_STRONG_INLINE
sublapl_expr<gsFeSolution<T> > sublapl(const gsFeSolution<T> & u, short_t startOffset, short_t endOffset)
{ return sublapl_expr<gsFeSolution<T> >(u,startOffset,endOffset); }


}
}

using namespace gismo;

int main(int argc, char *argv[])
{
    //! [Parse command line]
    bool plot = false;
    index_t numRefine  = 5;
    index_t numElevate = 0;
    index_t numRefineT = 0;
    index_t numElevateT= 0;
    index_t steps = 10;
    bool last{false}, export_b64{false};

    gsCmdLine cmd("Tutorial on solving a Poisson problem.");
    cmd.addInt( "e", "degreeElevation",
                "Number of degree elevation steps to perform before solving (0: equalize degree in all directions)", numElevate );
    cmd.addInt( "r", "uniformRefine", "Number of Uniform h-refinement loops",  numRefine );
    cmd.addInt( "R", "uniformRefineT", "Number of Uniform h-refinement loops in time",  numRefineT );
    cmd.addInt( "E", "degreeElevationT",
                "Number of degree elevation steps to perform in time before solving (0: equalize degree in all directions)", numElevateT );
    cmd.addInt( "N", "steps", "Number of time steps", steps );
    cmd.addSwitch("last", "Solve solely for the last level of h-refinement",
                  last);
    cmd.addSwitch(
        "plot", "Create a ParaView visualization file with the solution", plot);
    cmd.addSwitch("binary", "Use B64 encoding for Paraview", export_b64);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }
    //! [Parse command line]

    //! [Read input file]

    real_t t_max = 10;
    real_t alpha = 5e-3;
    real_t lambda = 1.0;

    // directions 0 and 1 are spatial
    // direction 2 is time
    gsMultiPatch<> mp;
    mp.addPatch(*gsNurbsCreator<>::BSplineCube());
    // mp.patch(0).coefs().col(2)*=t_max;

    std::string r = "10*sqrt((x-0.5)^2+(y-0.5)^2)";
    gsFunctionExpr<> initial("if("+r+" <= 1,(1-("+r+")^2)^2,0)",3);
    // gsFunctionExpr<> initial("(1-("+r+")^2)^2",3);
    gsInfo<<"Initial condition function "<< initial << "\n";
    // front: initial condition
    gsBoundaryConditions<> bc;
    bc.addCondition(boundary::front, condition_type::dirichlet, &initial);
    bc.addCondition(boundary::north, condition_type::dirichlet, 0);
    bc.addCondition(boundary::east , condition_type::dirichlet, 0);
    bc.addCondition(boundary::south, condition_type::dirichlet, 0);
    bc.addCondition(boundary::west , condition_type::dirichlet, 0);
    gsInfo<<"Boundary conditions:\n"<< bc <<"\n";
    bc.setGeoMap(mp);


    mp.degreeElevate(numElevate,0);
    mp.degreeElevate(numElevate,1);
    mp.degreeElevate(numElevateT,2);
    // h-refine each basis
    for (int r =0; r < numRefine; ++r)
    {
        mp.uniformRefine(1,1,0);
        mp.uniformRefine(1,1,1);
    }
    for (int r =0; r < numRefineT; ++r)
        mp.uniformRefine(1,1,2);

    //! [Refinement]
    gsMultiBasis<> dbasis(mp, true);//true: poly-splines (not NURBS)


    // Set heat conduction coefficient
    typedef gsExprAssembler<>::geometryMap geometryMap;
    typedef gsExprAssembler<>::space       space;
    typedef gsExprAssembler<>::solution    solution;

    gsExprAssembler<> A(1,1);

    gsInfo<<"Active options:\n"<< A.options() <<"\n";

    A.setIntegrationElements(dbasis);

    gsExprEvaluator<> ev(A);

    // Set the geometry map
    geometryMap G = A.getMap(mp);

    // Set the discretization space
    space w = A.getSpace(dbasis);
    w.setup(bc, dirichlet::l2Projection, 0);

    auto dwdt = subgrad(w,2,0);
    auto dwdX = subgrad(w,0,1);
    auto d2wdX= sublapl(w,0,1);
    // auto d2wdX= lapl(w);

    gsMatrix<> C;
    auto c  = A.getSolution(w, C);
    // auto dc = A.getSolution(u, dC);

    gsWarn<<"How to handle ijac?\n";
    auto dcdt = subgrad(c,2,0);
    auto dcdX = subgrad(c,0,1);
    auto d2cdX= sublapl(c,0,1);

    // Advection coefficient (Burgers)
    gsVector<> Fvec(2);
    // Fvec.setConstant(0.5);
    Fvec<<0,0;
    gsConstantFunction<> Ffun(Fvec,3);
    auto Fc = A.getCoeff(Ffun, G);
    // auto F  = Fc.tr();
    // auto dF = 0.0*w*Fc.tr();

    // Diffusion coefficient (constant)
    auto D = alpha;

    auto forcing      = 0*w;
    auto residual     = w*dcdt.tr() - dwdX*Fc*c*c  + D*dwdX*dcdX.tr();
    auto jacobian     = w*dwdt.tr() - dwdX*Fc*2*c*w.tr() + D*dwdX*dwdX.tr();

    // auto forcing      = -dwdX*Fc;
    // auto residual     = w*dcdt.tr() - dwdX*Fc + D*dwdX*dcdX.tr();
    // auto jacobian     = w*dwdt.tr() + D*dwdX*dwdX.tr();

    // Define linear solver (install SuperLUMT-devel)
#ifdef GISMO_WITH_PARDISO
    gsSparseSolver<>::PardisoLDLT solver;
#else
    gsSparseSolver<>::QR solver;
#endif

    gsMatrix<> dC;
    gsMatrix<> Q;

    A.initSystem();
    C.setZero(A.numDofs(),1);

    // Linear solve
    A.assemble(jacobian * meas(G));
    A.assemble(forcing  * meas(G));
    Q = A.rhs();
    solver.compute(A.matrix());
    C = solver.solve(Q);
    ev.writeParaview(c,G,"initial");

    // w.setup(bc, dirichlet::homogeneous, 0);
    index_t maxIter = 10;
    for (index_t iter = 0; iter < maxIter; ++iter)
    {
        A.clearMatrix();
        A.clearRhs();

        // Assemble residual
        A.assemble(residual * meas(G));
        Q = A.rhs();
        if (A.rhs().norm() < 1e-6 || iter == maxIter-1)
        {
            gsInfo<<"Converged with residual = "<<Q.norm()<<"\n";
            break;
        }

        gsInfo<<"Iteration "<<iter<<": residual = "<<Q.norm()<<", previous update = "<<dC.norm()<<"\n";

        // Assemble Jacobian
        // A.assembleJacobian(residual * meas(G), c);
        A.assemble(jacobian * meas(G));

        // Solve
        solver.compute(A.matrix());
        dC = solver.solve(-Q);
        C += dC;
    }

    w.setup(bc, dirichlet::l2Projection, 0);

    if (plot)
    {
        gsInfo<<"Plotting in Paraview...\n";

        gsParaviewCollection collection("ParaviewOutput/solution", &ev);
        collection.options().setSwitch("plotElements", true);
        collection.options().setSwitch("base64", export_b64);
        collection.options().setInt("plotElements.resolution", 16);
        collection.options().setInt("numPoints",1e5);
        collection.newTimeStep(&mp);
        collection.addField(c,"numerical solution");
        collection.saveTimeStep();
        collection.save();
    }
    return  EXIT_SUCCESS;
}
