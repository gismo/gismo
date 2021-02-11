/** @file quasiinterpolation_example.cpp

    @brief Different Quasi-Interpolation Schemes.

    Different quasi-interpolation-methods of a function.
    Some of the implemented methods are based on
    "Tom Lyche and Knut Morken. Spline methods draft. 2011"
    http://www.uio.no/studier/emner/matnat/ifi/INF-MAT5340/v11/undervisningsmateriale/

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): M. Haberleitner, A. Mantzaflaris, H. Verhelst
**/

#include <gismo.h>
#include <gsSolver/gsSolverUtils.h>

using namespace gismo;

//samples 'numVals' equally distributed points in the parameter
//domain, evaluate the function and

//the quasi-interpolant at those points and sum up the squared
//distances. return the squareroot of the sum as error
template<typename T>
T computeError(const gsFunction<T> &fun, const gsFunction<T> &spl, const int &numVals)
{
    gsGridIterator<T,CUBE> pt(spl.support(), numVals);
    gsMatrix<T> funExact, funApprox;
    T error = 0;
    for(index_t c = 0; pt; ++pt, ++c)
    {
        funExact  = fun.eval(*pt);
        funApprox = spl.eval(*pt);
        error    += (funExact - funApprox).squaredNorm();
    }
    return math::sqrt(error);
}

//compute a quasi-interpolant of some 'type' and print the error to the original function
//then (uniformly) refine the interpolant 'numRef' times and print the error each time
//type...1-> Schoenberg, 2->Taylor, 3->evaluation based with special implementations (based on degree),
//       0-> evaluation based with general formula
template<typename T>
bool errorAnalysis(const gsFunction<T> &fun, const gsBasis<T> & bbasis, int type, int numRef)       //ToDo: Ratio-test, return true or false..
{
    typename gsBasis<T>::uPtr basis = bbasis.clone();
    gsMatrix<T> coefs;
    real_t error, expConvRate=-1.0;// prevError=0.0, ratio=-1.0;
    std::vector<real_t> error_list, h_list;

    gsMatrix<> ab = basis->support();

    int deg = basis->minDegree();

    for(int i=0; i<numRef; i++)
    {
        switch(type)
        {
        case(1): gsQuasiInterpolate<T>::Schoenberg(*basis, fun, coefs); expConvRate = 2.0; break;
        case(2): gsQuasiInterpolate<T>::Taylor(*basis, fun, deg, coefs); expConvRate = (deg+1); break;
        case(3): gsQuasiInterpolate<T>::EvalBased(*basis, fun, true, coefs); expConvRate = (deg+1); break;
        case(0): gsQuasiInterpolate<T>::EvalBased(*basis, fun, false, coefs); expConvRate = (deg+1); break;
        case(4): gsQuasiInterpolate<T>::localIntpl(*basis, fun, coefs); expConvRate = (deg+1); break;
        default: GISMO_ERROR("invalid option");
        }

        error = computeError(fun, *basis->makeGeometry(give(coefs)), 100);
        error_list.push_back(error);

        // int numInnerKnots = basis->knots().size() - 2*(basis->degree()+1);
        // real_t h = (b-a)/(numInnerKnots+1);
        real_t h = math::pow((ab.col(1)-ab.col(0)).prod()/basis->numElements(),
                             1.0/basis->domainDim());
        h_list.push_back(h);
        basis->uniformRefine();
    }
    //fit a linear curve to the errors vs meshsize and use the slope of the curve as convergence rate
    real_t convRateAvg = gsSolverUtils<>::convergenceRateLS(error_list,h_list);

    //if the computed convergence rate is at least 85% of the expected one, we consider this test as passed
    real_t tolerance = 0.85;

    gsInfo<<"The convergence rate is "<< convRateAvg << " (expected rate = "<< expConvRate<< ")\n";
    gsInfo<< (convRateAvg > expConvRate*tolerance ? "OK" : "Not OK") <<"\n";
    return (convRateAvg > expConvRate*tolerance);
}


//compute a quasi-interpolant of some 'type' and print the error to the original function
//check if the error at 100 uniformly distributed points in the domain is "small enough"
//polynomials up the the degree fun.degree() should be reconstructed
//type...1-> Schoenberg, 2->Taylor, 3->evaluation based with special implementations (based on degree),
//       0-> evaluation based with general formula
//degPoly... degree of the polynomial to approximate
template<typename T>
bool polyReconstructionTaylor(const gsFunction<T> &fun, gsBSplineBasis<T> basis, int degPoly)
{
    real_t error = 0, expError = 0;

    T a = basis.knot(0);
    T b = basis.knot(basis.knots().size()-1);

    int deg = basis.degree();

    if(degPoly <= deg)  //set the expected error values
        expError = 1.0/math::pow(10.0,5);   //only 10^-5 because gsFunctionExpr.eval_into is not more accurate
    else
        expError = math::pow(10.0,8);    //some high value, for when we don't expect an exact approximation
                                        //(when degPoly is higher than the degree of the spline function we use to approximate)

    gsMatrix<T> coefs;
    gsQuasiInterpolate<T>::Taylor(basis, fun, deg, coefs);
    gsBSpline<T> result(basis,coefs);

    error = computeError(fun, result, 100);
    int numInnerKnots = basis.knots().size() - 2*(basis.degree()+1);
    real_t h = (b-a)/(numInnerKnots+1);
    gsInfo<<"h = "<< h <<",  err = "<<error <<"\n";

    gsInfo<< (error < expError ? "OK" : "Not OK") <<"\n";
    return (error < expError);

}

// A sinus(x)
template<typename T>
class gsMySinus: public gismo::gsFunction<T>
{
public:
    /// Shared pointer for gsMySinus
    typedef memory::shared_ptr< gsMySinus > Ptr;

    /// Unique pointer for gsMySinus
    typedef memory::unique_ptr< gsMySinus > uPtr;

    GISMO_CLONE_FUNCTION(gsMySinus)

    void eval_into(const gsMatrix<T> &u, gsMatrix<T> &result) const
    {
        result.resize(1,u.size());
        for(int i=0; i<u.size(); i++)
        {
            T val = math::sin(u.at(i));
            result(i) = val;
        }
    }


    void evalAllDers_into(const gsMatrix<T> &u, int n, std::vector<gsMatrix<T> > &result) const
    {
        gsMatrix<T> sin(1,u.size());
        gsMatrix<T> cos(1,u.size());
        for(int i=0; i<u.size(); i++)
        {
            T val1 = math::sin(u.at(i));
            sin(i) = val1;
            T val2 = math::cos(u.at(i));
            cos(i) = val2;
        }
        int type;
        for(int i=0; i<=n; i++)
        {
            type = i % 4;
            switch(type)
            {
            case(0): result.push_back(sin); break;
            case(1): result.push_back(cos); break;
            case(2): result.push_back(-sin); break;
            case(3): result.push_back(-cos); break;
            default: GISMO_ERROR("wrong type");
            }
        }

    }

    short_t domainDim() const {return 1;}

    /// Prints the object as a string.
    std::ostream &print(std::ostream &os) const
    { os << "MySinus(x) = sin(x)"; return os; }
};

//
template<typename T>
class gsMyCircle: public gismo::gsFunction<T>
{
public:
    /// Shared pointer for gsMyCircle
    typedef memory::shared_ptr< gsMyCircle > Ptr;

    /// Unique pointer for gsMyCircle
    typedef memory::unique_ptr< gsMyCircle > uPtr;

    GISMO_CLONE_FUNCTION(gsMyCircle)

    void eval_into(const gsMatrix<T> &u, gsMatrix<T> &result) const
    {
        result.resize(2,u.size());
        for(int i=0; i<u.size(); i++)
        {
            T valX = math::cos(2*EIGEN_PI*u.at(i));
            T valY = math::sin(2*EIGEN_PI*u.at(i));
            result(0,i) = valX;
            result(1,i) = valY;
        }
    }

    void evalAllDers_into(const gsMatrix<T> &u, int n, std::vector<gsMatrix<T> > &result) const
    {
        gsMatrix<T> sincos(2,u.size());
        gsMatrix<T> cossin(2,u.size());
        for(int i=0; i<u.size(); i++)
        {
            sincos(0,i) =
            cossin(1,i) = math::cos(2*EIGEN_PI*u.at(i));
            sincos(1,i) =
            cossin(0,i) = math::sin(2*EIGEN_PI*u.at(i));

        }
        int type;
        for(int i=0; i<=n; i++)
        {
            type = i % 4;
            switch(type)
            {
            case(0): result.push_back(std::pow(2*EIGEN_PI,i)*sincos); break;
            case(1): result.push_back(std::pow(2*EIGEN_PI,i)*cossin); break;
            case(2): result.push_back(-std::pow(2*EIGEN_PI,i)*sincos); break;
            case(3): result.push_back(-std::pow(2*EIGEN_PI,i)*cossin); break;
            default: GISMO_ERROR("wrong type");
            }
        }
    }

    short_t domainDim() const {return 1;}

    short_t targetDim() const {return 2;}

};


bool qi_1D()
{
    gsInfo<<"\n******** Running QI-1D ********\n";

    gsMySinus<real_t> mySinus;
    gsFunctionExpr<> myPolyLin("50*x-28",1);
    gsFunctionExpr<> myPolyQuad("-5*x^2 + 22*x - 4",1);
    gsMyCircle<real_t> myCircle;

    int deg1 = 1;
    int deg2 = 2;
    int deg3 = 3;

    real_t a = 0; // starting knot
    real_t b = 1; // ending knot
    unsigned interior = 1; // number of interior knots
    unsigned multEnd1 = deg1+1; // multiplicity at the two end knots
    unsigned multEnd2 = deg2+1; // multiplicity at the two end knots
    unsigned multEnd3 = deg3+1; // multiplicity at the two end knots

    gsKnotVector<> kv1(a, b, interior, multEnd1);
    gsKnotVector<> kv2(a, b, interior, multEnd2);
    gsKnotVector<> kv3(a, b, interior, multEnd3);

    gsBSplineBasis<> bas1(kv1);
    gsBSplineBasis<> bas2(kv2);
    gsBSplineBasis<> bas3(kv3);


    int numRef = 10;
    bool passed = true;


// ---------  Convergence-rate test with gsMySinus (dim = 1)

    gsInfo<<"\nLocal interpolation-based error analysis (cubic):\n";
    passed &= errorAnalysis<real_t>(mySinus, bas3, 4, numRef);

    gsInfo<<"\nSchoenberg error analysis (linear):\n";
    passed &= errorAnalysis<real_t>(mySinus, bas1, 1, numRef);

    gsInfo<<"\nTaylor error analysis (linear):\n";
    passed &= errorAnalysis(mySinus, bas1, 2, numRef);

    gsInfo<<"\nEvaluation-based error analysis (linear, special case):\n";
    passed &= errorAnalysis(mySinus, bas1, 3, numRef);

    gsInfo<<"\nEvaluation-based error analysis (linear, general):\n";
    passed &= errorAnalysis(mySinus, bas1, 0, numRef);

    gsInfo<<"\nSchoenberg error analysis (quadratic):\n";
    passed &= errorAnalysis(mySinus, bas2, 1, numRef);

    gsInfo<<"\nTaylor error analysis (quadratic):\n";
    passed &= errorAnalysis(mySinus, bas2, 2, numRef);

    gsInfo<<"\nEvaluation-based error analysis (quadratic, special case):\n";
    passed &= errorAnalysis(mySinus, bas2, 3, numRef);

    gsInfo<<"\nEvaluation-based error analysis (quadratic, general):\n";
    passed &= errorAnalysis(mySinus, bas2, 0, numRef);


    gsInfo<<"\nSchoenberg error analysis (cubic):\n";
    passed &= errorAnalysis(mySinus, bas3, 1, numRef);

    gsInfo<<"\nTaylor error analysis (cubic):\n";
    passed &= errorAnalysis(mySinus, bas3, 2, numRef);

    gsInfo<<"\nEvaluation-based error analysis (cubic, special case):\n";
    passed &= errorAnalysis(mySinus, bas3, 3, numRef);

    gsInfo<<"\nEvaluation-based error analysis (cubic, general):\n";
    passed &= errorAnalysis(mySinus, bas3, 0, numRef);


// ---------  Convergence-rate test with gsMyCircle (dim = 2)

    gsInfo<<"\nSchoenberg error analysis (linear):\n";
    passed &= errorAnalysis<real_t>(myCircle, bas1, 1, numRef);

    gsInfo<<"\nTaylor error analysis (linear):\n";
    passed &= errorAnalysis(myCircle, bas1, 2, numRef);

    gsInfo<<"\nEvaluation-based error analysis (linear, special case):\n";
    passed &= errorAnalysis(myCircle, bas1, 3, numRef);

    gsInfo<<"\nEvaluation-based error analysis (linear, general):\n";
    passed &= errorAnalysis(myCircle, bas1, 0, numRef);


    gsInfo<<"\nSchoenberg error analysis (quadratic):\n";
    passed &= errorAnalysis(myCircle, bas2, 1, numRef);

    gsInfo<<"\nTaylor error analysis (quadratic):\n";
    passed &= errorAnalysis(myCircle, bas2, 2, numRef);

    gsInfo<<"\nEvaluation-based error analysis (quadratic, special case):\n";
    passed &= errorAnalysis(myCircle, bas2, 3, numRef);

    gsInfo<<"\nEvaluation-based error analysis (quadratic, general):\n";
    passed &= errorAnalysis(myCircle, bas2, 0, numRef);


    gsInfo<<"\nSchoenberg error analysis (cubic):\n";
    passed &= errorAnalysis(myCircle, bas3, 1, numRef);

    gsInfo<<"\nTaylor error analysis (cubic):\n";
    passed &= errorAnalysis(myCircle, bas3, 2, numRef);

    gsInfo<<"\nEvaluation-based error analysis (cubic, special case):\n";
    passed &= errorAnalysis(myCircle, bas3, 3, numRef);

    gsInfo<<"\nEvaluation-based error analysis (cubic, general):\n";
    passed &= errorAnalysis(myCircle, bas3, 0, numRef);


// --------- Taylor Polynomial-reconstruction test (dim = 1)

    gsInfo<<"\nTaylor  poly-reconstruction (linear) [deg=1]:\n";
    passed &= polyReconstructionTaylor(myPolyLin, bas1, 1);

    gsInfo<<"\nTaylor poly-reconstruction (quadratic) [deg=1]:\n";
    passed &= polyReconstructionTaylor(myPolyLin, bas2, 1);

    gsInfo<<"\nTaylor poly-reconstruction (linear) [deg=2]:\n";
    passed &= polyReconstructionTaylor(myPolyQuad, bas1, 2);

    gsInfo<<"\nTaylor poly-reconstruction (quadratic)  [deg=2]:\n";
    passed &= polyReconstructionTaylor(myPolyQuad, bas2, 2);

    return passed;
}


bool qi_2D()
{
    gsInfo<<"\n******** Running QI-2D ********\n";

    gsFunctionExpr<real_t> mySinus("sin(x)*cos(y)",2);
    int deg1 = 1;
    int deg2 = 2;
    int deg3 = 3;
    real_t a = 0; // starting knot
    real_t b = 1; // ending knot
    unsigned interior = 1; // number of interior knots
    unsigned multEnd1 = deg1+1; // multiplicity at the two end knots
    unsigned multEnd2 = deg2+1; // multiplicity at the two end knots
    unsigned multEnd3 = deg3+1; // multiplicity at the two end knots

    gsKnotVector<> kv1(a, b, interior, multEnd1);
    gsKnotVector<> kv2(a, b, interior, multEnd2);
    gsKnotVector<> kv3(a, b, interior, multEnd3);

    gsTensorBSplineBasis<2> bas1(kv1,kv1);
    gsTensorBSplineBasis<2> bas2(kv2,kv2);
    gsTensorBSplineBasis<2> bas3(kv3,kv3);

    int numRef = 5;
    bool passed = true;

// ---------  Convergence-rate test for trigonometric function

    gsInfo<<"\nLocal interpolation-based error analysis (linear):\n";
    passed &= errorAnalysis<real_t>(mySinus, bas1, 4, numRef);

    gsInfo<<"\nLocal interpolation-based error analysis (quadratic):\n";
    passed &= errorAnalysis<real_t>(mySinus, bas2, 4, numRef);

    gsInfo<<"\nLocal interpolation-based error analysis (cubic):\n";
    passed &= errorAnalysis<real_t>(mySinus, bas3, 4, numRef);

    gsInfo<<"\nSchoenberg error analysis (linear):\n";
    passed &= errorAnalysis<real_t>(mySinus, bas1, 1, numRef);

    gsInfo<<"\nSchoenberg error analysis (quadratic):\n";
    passed &= errorAnalysis(mySinus, bas2, 1, numRef);

    gsInfo<<"\nSchoenberg error analysis (cubic):\n";
    passed &= errorAnalysis(mySinus, bas3, 1, numRef);

    return passed;
}


bool qi_hs_2D()
{
    gsInfo<<"\n******** Running QI-Hierarchical-2D ********\n";

    gsFunctionExpr<real_t> mySinus("sin(x)*cos(y)",2);
    gsFunctionExpr<> myPolyLin("50*x- + 30*y + 28",2);
    gsFunctionExpr<> myPolyQuad("-5*x^2 + 3*x*y + y^2 + 22*x + y - 4",2);

    int deg1 = 1;
    int deg2 = 2;
    int deg3 = 3;

    real_t a = 0; // starting knot
    real_t b = 1; // ending knot
    unsigned interior = 1; // number of interior knots
    unsigned multEnd1 = deg1+1; // multiplicity at the two end knots
    unsigned multEnd2 = deg2+1; // multiplicity at the two end knots
    unsigned multEnd3 = deg3+1; // multiplicity at the two end knots

    gsKnotVector<> kv1(a, b, interior, multEnd1);
    gsKnotVector<> kv2(a, b, interior, multEnd2);
    gsKnotVector<> kv3(a, b, interior, multEnd3);

    gsTensorBSplineBasis<2> bas1(kv1,kv1);
    gsTensorBSplineBasis<2> bas2(kv2,kv2);
    gsTensorBSplineBasis<2> bas3(kv3,kv3);

    int numRef = 5;
    bool passed = true;

    gsTHBSplineBasis<2,real_t> thb1,thb2,thb3;
    thb1 = gsTHBSplineBasis<2,real_t>(bas1);
    thb2 = gsTHBSplineBasis<2,real_t>(bas2);
    thb3 = gsTHBSplineBasis<2,real_t>(bas3);

    gsMatrix<> refBoxes(2,2);
    refBoxes.col(0) << 0,0;
    refBoxes.col(1) << 0.25,0.25;
    thb1.refine( refBoxes );
    thb2.refine( refBoxes );
    thb3.refine( refBoxes );

    //gsWriteParaview(thb1,"basis");


// ---------  Convergence-rate test for trigonometric function

    gsInfo<<"\nLocal interpolation-based error analysis (linear):\n";
    passed &= errorAnalysis<real_t>(mySinus, thb1, 4, numRef);

    gsInfo<<"\nLocal interpolation-based error analysis (quadratic):\n";
    passed &= errorAnalysis<real_t>(mySinus, thb2, 4, numRef);

    gsInfo<<"\nLocal interpolation-based error analysis (cubic):\n";
    passed &= errorAnalysis<real_t>(mySinus, thb3, 4, numRef);


    gsMatrix<> coefs;
    gsQuasiInterpolate<real_t>::localIntpl(thb3, mySinus, coefs);
    gsTHBSpline<2> aa1(thb3,coefs);
    gsInfo<<"\nCheck projection (cubic): ";
    gsQuasiInterpolate<real_t>::localIntpl(thb3, aa1, coefs);
    gsTHBSpline<2> aa2(thb3,coefs);
    gsInfo<<"error="<< computeError(aa1,aa2,100) <<"\n";
                
    return passed;
}

bool qi_3D()
{
    gsInfo<<"\n******** Running QI-3D ********\n";

    gsFunctionExpr<real_t> mySinus("sin(x)*cos(y)*sin(z)",3);
    int deg1 = 1;
    int deg2 = 2;
    int deg3 = 3;
    real_t a = 0; // starting knot
    real_t b = 1; // ending knot
    unsigned interior = 1; // number of interior knots
    unsigned multEnd1 = deg1+1; // multiplicity at the two end knots
    unsigned multEnd2 = deg2+1; // multiplicity at the two end knots
    unsigned multEnd3 = deg3+1; // multiplicity at the two end knots

    gsKnotVector<> kv1(a, b, interior, multEnd1);
    gsKnotVector<> kv2(a, b, interior, multEnd2);
    gsKnotVector<> kv3(a, b, interior, multEnd3);

    gsTensorBSplineBasis<3> bas1(kv1,kv1,kv1);
    gsTensorBSplineBasis<3> bas2(kv2,kv2,kv2);
    gsTensorBSplineBasis<3> bas3(kv3,kv3,kv3);

    int numRef = 5;
    bool passed = true;


// ---------  Convergence-rate test for trigonometric function

    gsInfo<<"\nLocal interpolation-based error analysis (linear):\n";
    passed &= errorAnalysis<real_t>(mySinus, bas1, 4, numRef);

    gsInfo<<"\nLocal interpolation-based error analysis (quadratic):\n";
    passed &= errorAnalysis<real_t>(mySinus, bas2, 4, numRef);

    // gsInfo<<"\nLocal interpolation-based error analysis (cubic):\n";
    // passed &= errorAnalysis<real_t>(mySinus, bas3, 4, numRef);

    gsInfo<<"\nSchoenberg error analysis (linear):\n";
    passed &= errorAnalysis<real_t>(mySinus, bas1, 1, numRef);

    gsInfo<<"\nSchoenberg error analysis (quadratic):\n";
    passed &= errorAnalysis(mySinus, bas2, 1, numRef);

    // gsInfo<<"\nSchoenberg error analysis (cubic):\n";
    // passed &= errorAnalysis(mySinus, bas3, 1, numRef);

    return passed;
}


int main(int argc, char* argv[])
{
    bool passed = true;
    passed &=  qi_1D();
    passed &=  qi_2D();
    passed &=  qi_3D();

    passed &=  qi_hs_2D();
    return passed ? 0 : 1;
}
