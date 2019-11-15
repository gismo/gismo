/** @file gsExprIntegral_test.cpp

    @brief Testing integral computation using the expression evaluator

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris
*/

#include <gismo.h>

#include <gsAssembler/gsExprEvaluator.h>

using namespace gismo;

/*
    gsIntegrantZ has as input
        -
        -
        -

*/
template<class T>
class gsIntegrantZ : public gismo::gsFunction<T>
{
    protected:
        const typename gsFunction<T>::Ptr _fun;
        mutable gsMatrix<T> tmp;
        gsMatrix<T> _surfPts;
    public:
        /// Shared pointer for gsIntegrantZ
        typedef memory::shared_ptr< gsIntegrantZ > Ptr;

        /// Unique pointer for gsIntegrantZ
        typedef memory::unique_ptr< gsIntegrantZ > uPtr;

    // copy constructor
    explicit gsIntegrantZ(const gsIntegrantZ &other) : _fun(other._fun), _surfPts(other._surfPts) {}

    gsIntegrantZ(const gsFunction<T> & fun) : _fun(fun.clone()) { }

    GISMO_CLONE_FUNCTION(gsIntegrantZ)

    void setPoint(const gsMatrix<T>& surfPts) { _surfPts = surfPts; }

    gsMatrix<T> point() { return _surfPts; }

    short_t domainDim() const {return 1;}

    short_t targetDim() const {return _fun->targetDim();}

    void eval_into(const gsMatrix<T>& u, gsMatrix<T>& result) const
    {
        GISMO_ASSERT(u.rows()==1,
            "The number of rows for the 1D coordinate is not 1 but " <<u.rows()<<"!"  );
        GISMO_ASSERT(_fun->domainDim()==_surfPts.rows() + 1,
            "The domain dimensions do not match! fun.domainDim() != surfPts.rows() + 1! (" <<_fun->domainDim()<<"!="<<_surfPts.rows() + 1<<" )"  );
        GISMO_ASSERT(_surfPts.cols()==1,
            "Multiple ("<<_surfPts.cols()<<") parametric points given, accepts only 1... " <<"!"     );

        index_t m = _surfPts.rows();
        index_t N = u.cols();

        tmp.resize(m + 1,N);
        tmp.topRows(m)=_surfPts.replicate(1,N);
        tmp.bottomRows(1) = u;

        _fun->eval_into(tmp,result);
    }

    std::ostream &print(std::ostream &os) const
      { os << "gsIntegrantZ ( " << _fun << " )"; return os; };
};

template <class T>
class gsIntegrateZ : public gismo::gsFunction<T>
{
    protected:
        const typename gsFunction<T>::Ptr _fun;//, _t;
        mutable gsMatrix<T> tmp;
        gsMatrix<T> _surfPts;
        T _t;
    public:
        /// Shared pointer for gsIntegrateZ
        typedef memory::shared_ptr< gsIntegrateZ > Ptr;

        /// Unique pointer for gsIntegrateZ
        typedef memory::unique_ptr< gsIntegrateZ > uPtr;

    // copy constructor
    explicit gsIntegrateZ(const gsIntegrateZ &other)
    : _fun(other._fun), _t(other._t), _surfPts(other._surfPts) {}

    gsIntegrateZ(const gsFunction<T> & fun, T thickness)
    : _fun(fun.clone()), _t(thickness)
    {
        // gsConstantFunction _t(thickness, 3);
    }

    // gsIntegrateZ(const gsFunction<T> & fun, gsFunction<T> & thickFun)
    // : _fun(fun.clone()), _t(thickFun.clone()) { }

    GISMO_CLONE_FUNCTION(gsIntegrateZ)

    short_t domainDim() const {return 1;}

    short_t targetDim() const {return _fun->targetDim();}

    void setPoint(const gsMatrix<T>& surfPts) { _surfPts = surfPts; }

    // u are z-coordinates only!
    void eval_into(const gsMatrix<T>& u, gsMatrix<T>& result) const
    {
        // // Compute the thickness
        // _tmp.points = u;
        // _t.eval_into(_tmp.values[0], thickMat);

        result.resize(_fun->targetDim(),u.cols());

        // Define integrator for the z direction
        gsExprEvaluator<real_t> ev;

        // Define integration interval
        int k = 1; // interior knots
        int p = 1; // B-spline order

        // Make 1D domain with a basis
        gsKnotVector<> KV(-_t/2.0, _t/2.0, k, p+1);
        gsMultiBasis<> basis;
        basis.addBasis(gsBSplineBasis<>::make(KV));

        // Set integration elements along basis
        ev.setIntegrationElements(basis);

        // Define integrant variables
        typedef gsExprEvaluator<real_t>::variable    variable;
        variable    integrant = ev.getVariable(*_fun, 1);

        for (index_t i=0; i!=_fun->targetDim(); ++i)
            for (index_t j = 0; j != u.cols(); ++j)
            {
                //thickness integral for all components i.
                ev.integral(integrant.tr()[i]);
                result(i,j) = ev.value();
            }
    }

    std::ostream &print(std::ostream &os) const
      { os << "gsIntegrateZ ( " << _fun << " )"; return os; };
};

template <class T>
class gsIntegrate : public gismo::gsFunction<T>
{
    protected:
        const typename gsFunction<T>::Ptr _fun, _t;
        mutable gsMatrix<T> tmp, thickMat;
        gsMatrix<T> _surfPts;
        // T _t;
    public:
        /// Shared pointer for gsIntegrate
        typedef memory::shared_ptr< gsIntegrate > Ptr;

        /// Unique pointer for gsIntegrate
        typedef memory::unique_ptr< gsIntegrate > uPtr;

    // copy constructor
    explicit gsIntegrate(const gsIntegrate &other)
    : _fun(other._fun), _t(other._t), _surfPts(other._surfPts) {}

    // gsIntegrate(const gsFunction<T> & fun, T thickness)
    // : _fun(fun.clone()), _t(thickness)
    // {
    //     gsConstantFunction _t(thickness, 3);
    // }

    gsIntegrate(const gsFunction<T> & fun, gsFunction<T> & thickFun)
    : _fun(fun.clone()), _t(thickFun.clone()) { }

    GISMO_CLONE_FUNCTION(gsIntegrate)

    short_t domainDim() const {return 2;}

    short_t targetDim() const {return _fun->targetDim();}

    // u are xy-coordinates only; domainDim=2
    void eval_into(const gsMatrix<T>& u, gsMatrix<T>& result) const
    {
        // Compute the thickness
        _t->eval_into(u, thickMat);

        result.resize(_fun->targetDim(),u.cols());

        // Define integrator for the z direction
        gsExprEvaluator<real_t> ev;

        // Define integrant variables
        typedef gsExprEvaluator<real_t>::variable    variable;
        gsIntegrantZ integrant(*_fun);

        T tHalf;
        for (index_t i=0; i!=_fun->targetDim(); ++i)
            for (index_t j = 0; j != u.cols(); ++j)
            {
                // this part is quite sloppy since a multi-basis is created every iteration..
                tHalf = thickMat(0,j)/2.0;

                gsKnotVector<> KV(-tHalf, tHalf, 2, 2);
                gsMultiBasis<> basis;
                basis.addBasis(gsBSplineBasis<>::make(KV));

                // Set integration elements along basis
                ev.setIntegrationElements(basis);

                variable intfun = ev.getVariable(integrant, 1);

                // set new integration point
                integrant.setPoint(u.col(j));

                //thickness integral for all components i.
                // ev.eval(intfun.coord(i));
                ev.integral(intfun.tr()[i]);
                result(i,j) = ev.value();
            }
    }

    std::ostream &print(std::ostream &os) const
      { os << "gsIntegrateZ ( " << _fun << " )"; return os; };
};

template <class T>
class gsMaterialMatrix : public gismo::gsFunction<T>
{
  // Computes the material matrix for different material models
  //
protected:
    const gsFunctionSet<T> * _mp;
    const gsFunction<T> * _YoungsModulus;
    const gsFunction<T> * _PoissonRatio;
    mutable gsMapData<T> _tmp;
    mutable gsMatrix<real_t,3,3> F0;
    mutable gsMatrix<T> Emat,Nmat;
    mutable real_t lambda, mu, E, nu, C_constant;

public:
    /// Shared pointer for gsMaterialMatrix
    typedef memory::shared_ptr< gsMaterialMatrix > Ptr;

    /// Unique pointer for gsMaterialMatrix
    typedef memory::unique_ptr< gsMaterialMatrix > uPtr;

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    gsMaterialMatrix(const gsFunctionSet<T> & mp, const gsFunction<T> & YoungsModulus,
                   const gsFunction<T> & PoissonRatio) :
    _mp(&mp), _YoungsModulus(&YoungsModulus), _PoissonRatio(&PoissonRatio), _mm_piece(nullptr)
    {
        _tmp.flags = NEED_JACOBIAN | NEED_NORMAL | NEED_VALUE;
    }

    ~gsMaterialMatrix() { delete _mm_piece; }

    GISMO_CLONE_FUNCTION(gsMaterialMatrix)

    short_t domainDim() const {return 2;}

    short_t targetDim() const {return 9;}

    mutable gsMaterialMatrix<T> * _mm_piece; // todo: improve the way pieces are accessed

    const gsFunction<T> & piece(const index_t k) const
    {
        delete _mm_piece;
        _mm_piece = new gsMaterialMatrix(_mp->piece(k), *_YoungsModulus, *_PoissonRatio);
        return *_mm_piece;
    }

    //class .. matMatrix_z
    // should contain eval_into(thickness variable)

    // Input is parametric coordinates of the surface \a mp
    void eval_into(const gsMatrix<T>& u, gsMatrix<T>& result) const
    {
        // NOTE 1: if the input \a u is considered to be in physical coordinates
        // then we first need to invert the points to parameter space
        // _mp.patch(0).invertPoints(u, _tmp.points, 1e-8)
        // otherwise we just use the input paramteric points
        gsDebugVar(u);
        _tmp.points = u;
        gsDebugVar(_mp->piece(0));
        gsDebugVar(_tmp.points);

        static_cast<const gsFunction<T>*>(_mp)->computeMap(_tmp);

        // NOTE 2: in the case that parametric value is needed it suffices
        // to evaluate Youngs modulus and Poisson's ratio at
        // \a u instead of _tmp.values[0].
        _YoungsModulus->eval_into(_tmp.values[0], Emat);
        _PoissonRatio->eval_into(_tmp.values[0], Nmat);

        result.resize( targetDim() , u.cols() );
        for( index_t i=0; i< u.cols(); ++i )
        {
            gsAsMatrix<T, Dynamic, Dynamic> C = result.reshapeCol(i,3,3);

            F0.leftCols(2) = _tmp.jacobian(i);
            F0.col(2)      = _tmp.normal(i).normalized();
            F0 = F0.inverse();
            F0 = F0 * F0.transpose(); //3x3

            // Evaluate material properties on the quadrature point
            E = Emat(0,i);
            nu = Nmat(0,i);
            lambda = E * nu / ( (1. + nu)*(1.-2.*nu)) ;
            mu     = E / (2.*(1. + nu)) ;

            C_constant = 4*lambda*mu/(lambda+2*mu);

            C(0,0) = C_constant*F0(0,0)*F0(0,0) + 2*mu*(2*F0(0,0)*F0(0,0));
            C(1,1) = C_constant*F0(1,1)*F0(1,1) + 2*mu*(2*F0(1,1)*F0(1,1));
            C(2,2) = C_constant*F0(0,1)*F0(0,1) + 2*mu*(F0(0,0)*F0(1,1) + F0(0,1)*F0(0,1));
            C(1,0) =
            C(0,1) = C_constant*F0(0,0)*F0(1,1) + 2*mu*(2*F0(0,1)*F0(0,1));
            C(2,0) =
            C(0,2) = C_constant*F0(0,0)*F0(0,1) + 2*mu*(2*F0(0,0)*F0(0,1));
            C(2,1) = C(1,2) = C_constant*F0(0,1)*F0(1,1) + 2*mu*(2*F0(0,1)*F0(1,1));

            //gsDebugVar(C);
        }
    }

    // std::ostream &print(std::ostream &os) const
    //   { os << "gsMaterialMatrix "; return os; };
};



int main(int argc, char *argv[])
{
    gsCmdLine cmd("Testing expression evaluator.");
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    gsMultiPatch<> mp;

    mp.addPatch( gsNurbsCreator<>::BSplineSquare(1) ); // degree
    mp.addAutoBoundaries();
    mp.embed(3);

    gsMultiBasis<> b(mp);
    b.uniformRefine(1);


    // Initiate the expression evaluator
    gsExprEvaluator<real_t> ev;

    // Set the parameter mesh as the integration mesh
    ev.setIntegrationElements(b);


    /*
        test gsIntegrantZ function
    */

    gsVector<> pt1D(1); pt1D.setConstant(0.25);
    gsVector<> pt2D(2); pt2D.setConstant(0.25);
    gsVector<> pt3D(3); pt3D.setConstant(0.25);

    gsMatrix<> points(2,11);
    points<<0,1,2,3,4,5,6,7,8,9,10,
            0,1,2,3,4,5,6,7,8,9,10;

    gsFunctionExpr<> fun("1*x","2*y","x*y*z^2",3);

    gsMatrix<> result;
    fun.eval_into(pt3D,result);
    gsInfo<<"Evaluation of the following function on point (x,y) = ("<<pt3D.at(0)<<","<<pt3D.at(1)<<") and z coordinate "<<pt3D.at(2)<<"\n";
    gsInfo<<fun<<"\n";
    gsInfo<<"result = "<<result.transpose()<<"\n";

    gsInfo<<"Evaluation of the following function on point (x,y) = ("<<pt2D.at(0)<<","<<pt2D.at(1)<<") and z coordinate "<<pt1D.at(0)<<"\n";
    gsInfo<<fun<<"\n";
    gsIntegrantZ fun2(fun);
    fun2.setPoint(pt2D); // if changes to be applied
    fun2.eval_into(pt1D,result);
    gsInfo<<"result = "<<result.transpose()<<"\n";

    pt2D.setConstant(0.1);
    gsInfo<<"Evaluation of the following function on point (x,y) = ("<<pt2D.at(0)<<","<<pt2D.at(1)<<") and z coordinate "<<pt1D.at(0)<<"\n";
    gsInfo<<fun<<"\n";
    fun2.setPoint(pt2D); // if changes to be applied
    fun2.eval_into(pt1D,result);
    gsInfo<<"result = "<<result.transpose()<<"\n";

    /*
        test gsIntegrateZ function
    */
    gsFunctionExpr<> fun3("1","x","x^2","x^3","x^4","x^5","x^6","x^7","x^8",1);

    real_t bound = 1.0;
    gsIntegrateZ<real_t> integrator(fun3,bound);
    integrator.eval_into(pt1D,result);

    gsInfo<<"Integration of the following function from "<<-bound/2.0<<" to "<<bound/2.0<<": \n";
    gsInfo<<fun3<<"\n";
    gsInfo<<"Result: "<<result.transpose()<<"\n";

    /*
        test gsIntegrateZ function
    */
    bound = 1.0;
    gsIntegrateZ integrator2(fun2,bound);
    integrator2.eval_into(pt1D,result);

    gsInfo<<"Integration of the third component of the following function from "<<-bound/2.0<<" to "<<bound/2.0<<" on point (x,y) = ("<<pt2D.at(0)<<","<<pt2D.at(1)<<") \n";
    gsInfo<<fun<<"\n";
    gsInfo<<"Result: "<<result.transpose()<<"\n";


    /*
        test gsIntegrate function
    */
    gsConstantFunction<> thickFun(bound,2);
    gsIntegrate integrate(fun,thickFun);
    integrate.eval_into(points,result);

    gsInfo<<"Integration of the third component of the following function from "<<-bound/2.0<<" to "<<bound/2.0<<"\n";
    gsInfo<<fun<<"\n";
    gsInfo<<"on points (x,y) = \n";
    gsInfo<<points.transpose()<<"\n";
    gsInfo<<"Result: \n"<<result.transpose()<<"\n";


    /*
        Integrate now a material matrix point by point
    */
    real_t E_modulus = 1.0;
    real_t PoissonRatio = 0.0;
    gsFunctionExpr<> E(std::to_string(E_modulus),3);
    gsFunctionExpr<> nu(std::to_string(PoissonRatio),3);
    gsMaterialMatrix materialMat(mp, E, nu);
    // gsIntegrate integrateMM(materialMat,thickFun);

    materialMat.eval_into(pt2D, result);
    // materialMat.eval_into(points, result);
    // integrateMM.eval_into(points,result);

    gsInfo<<"Result: \n"<<result.transpose()<<"\n";

    return EXIT_SUCCESS;
}
