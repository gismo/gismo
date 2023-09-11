/** @file polynomial_intersection.cpp

    @brief Demonstrates adaptive fitting of data samples

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): S. Imperatore , Y. Ji
*/

#include <gismo.h>

//#include <Eigen/Dense>
//#include <unsupported/Eigen/NonLinearOptimization>
//#include <unsupported/Eigen/NumericalDiff>

#include <lsqcpp.h>

using namespace gismo;

bool areBoundingBoxesIntersecting(const gsMatrix<>& box1, const gsMatrix<>& box2) {
  if (box1(0,0) > box2(0,1) || box2(0,0) > box1(0,1)) {
    return false;
  }
  if (box1(1,0) > box2(1,1) || box2(1,0) > box1(1,1)) {
    return false;
  }
  return true;
}

//// Generic functor
//// See http://gsEigen.tuxfamily.org/index.php?title=Functors
//// C++ version of a function pointer that stores meta-data about the function
//template<typename _Scalar, int NX = gsEigen::Dynamic, int NY = gsEigen::Dynamic>
//struct Functor
//{
//
//  // Information that tells the caller the numeric type (eg. double) and size (input / output dim)
//  typedef _Scalar Scalar;
//  enum { // Required by numerical differentiation module
//    InputsAtCompileTime = NX,
//    ValuesAtCompileTime = NY
//  };
//
//  // Tell the caller the matrix sizes associated with the input, output, and jacobian
//  typedef gsEigen::Matrix<Scalar,InputsAtCompileTime,1> InputType;
//  typedef gsEigen::Matrix<Scalar,ValuesAtCompileTime,1> ValueType;
//  typedef gsEigen::Matrix<Scalar,ValuesAtCompileTime,InputsAtCompileTime> JacobianType;
//
//  // Local copy of the number of inputs
//  int m_inputs, m_values;
//
//  // Two constructors:
//  Functor() : m_inputs(InputsAtCompileTime), m_values(ValuesAtCompileTime) {}
//  Functor(int inputs, int values) : m_inputs(inputs), m_values(values) {}
//
//  // Get methods for users to determine function input and output dimensions
//  int inputs() const { return m_inputs; }
//  int values() const { return m_values; }
//
//};
//
//// https://en.wikipedia.org/wiki/Test_functions_for_optimization
//// Booth Function
//// Implement f(x,y) = (x + 2*y -7)^2 + (2*x + y - 5)^2
//struct BoothFunctor : Functor<double>
//{
//  // Simple constructor
//  BoothFunctor(): Functor<double>(2,2) {}
//
//  // Implementation of the objective function
//  int operator()(const gsEigen::VectorXd &z, gsEigen::VectorXd &fvec) const {
//    double x = z(0);   double y = z(1);
//    /*
//     * Evaluate the Booth function.
//     * Important: LevenbergMarquardt is designed to work with objective functions that are a sum
//     * of squared terms. The algorithm takes this into account: do not do it yourself.
//     * In other words: objFun = sum(fvec(i)^2)
//     */
//    fvec(0) = x + 2*y - 7;
//    fvec(1) = 2*x + y - 5;
//    return 0;
//  }
//};

//template<typename Scalar>
typedef gsEigen::Matrix<real_t, gsEigen::Dynamic, 1> Vector;
typedef gsEigen::Matrix<real_t, gsEigen::Dynamic, gsEigen::Dynamic> Matrix;
typedef std::function<void(const Vector &, Vector &, Matrix &)> ErrorFunction;

ErrorFunction BoothFunctor(const Vector &z, Vector &fvec, Matrix &jac) {
    double x = z(0);   double y = z(1);
    fvec.resize(2);
    /*
     * Evaluate the Booth function.
     * Important: LevenbergMarquardt is designed to work with objective functions that are a sum
     * of squared terms. The algorithm takes this into account: do not do it yourself.
     * In other words: objFun = sum(fvec(i)^2)
     */
    fvec(0) = x + 2*y - 7;
    fvec(1) = 2*x + y - 5;
//  fvec(0) = x * x + y - 11;
//  fvec(1) = x + y * y - 7;
    return nullptr;
};

int main(int argc, char *argv[])
{
    // Options with default values

    std::string fn1 = "../filedata/cip/MBench_S1002_v3_curve_01.xml";
    std::string fn2 = "../filedata/cip/MBench_UIC60_v3_curve_01.xml";

    // Reading options from the command line
    gsCmdLine cmd("Intersection of two spline curves.");
    cmd.addString("c", "curve1", "Input curve 1", fn1);
    cmd.addString("s", "curve2", "Input curve 2", fn2);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }


    //gsGeometry<>::uPtr geo = filedata.getAnyFirst< gsGeometry<> >();
    gsFileData<> fd1(fn1);
    gsBSpline<>::uPtr geo = fd1.getAnyFirst< gsBSpline<> >();
    gsBSpline<> & curve1 = const_cast<gsBSpline<> &>(*geo);

    // index_t deg = curve1.degree();
    index_t deg = curve1.degree();
    index_t order = deg+1;
    gsInfo << "Bspline of degree " << deg << ", order "<< order << "\n";

    // gsKnotVector<> ku(0, 1, 5, deg+1);//start,end,interior knots, start/end multiplicites of knots
    // gsMatrix<> coefs(9, 2);
    //
    // coefs << 0, 0,
    //          0., 0.25,
    //          0.125, 0.5,
    //          0.25, 0.55,
    //          0.375, 0.5,
    //          0.5, 0.25,
    //          0.625, 0.25,
    //          0.8, 0.25,
    //          1, 0.5;
    //
    // gsMatrix<> coefs_plot;
    // coefs_plot = coefs.transpose();
    // gsWriteParaviewPoints(coefs_plot, "coefs_original");
    // gsBSpline<> curve1( ku, give(coefs));
    gsDebugVar(curve1.knots());


    gsKnotVector<> ku = curve1.knots();
    gsMatrix<> coefs = curve1.coefs();

    gsMatrix<> coefs_plot;
    coefs_plot = coefs.transpose();
    gsWriteParaviewPoints(coefs_plot, "coefs_original");









    gsInfo << "curve:\n" << curve1 << "\n";

//    gsWriteParaview( curve1, "c1", 1000);
    gsWriteParaview( curve1, "c1", 1000, false, true);


    // TODO: add bezier extraction, by proper knot insertion algorithm.
    auto multContainer = ku.multiplicities();
    gsInfo << multContainer[0] << "\n";
    //gsDebugVar(multContainer);
    // for(index_t i = 1; i < multContainer.size()-1; i++)
    // {
    //   c1.insertKnot(ku.knot(i), degree - multContainer[i]);
    // }
    //
    ku.increaseMultiplicity(deg-1, false);
    gsDebugVar(ku);

    gsMatrix<> X, uv;
    ku.greville_into(uv);
    gsDebugVar(uv); // problem with values outside the domain.
    uv = uv.cwiseMax(*ku.domainBegin()).cwiseMin(*ku.domainEnd()); // keep point within design bounds
    gsInfo << "Debug greville points with domain constraints.\n";
    gsInfo << "knot domain: [" << *ku.domainBegin() << ", " << *ku.domainEnd() << "]\n";
    gsInfo << uv.rows() << " x " << uv.cols() << "\n";
    gsDebugVar(uv);
    curve1.eval_into(uv, X);

    gsInfo << "Points:\n" << X.rows() << " x " << X.cols() << "\n";
    gsInfo << "Parameters:\n" << uv << "\n";
    gsWriteParaviewPoints(X, "points");


    gsBSplineBasis<> basis(ku);
    gsDebugVar(basis);

    gsFitting<> bezierExtraction(uv, X, basis);
    bezierExtraction.compute();

    gsMatrix<> newCoefs = bezierExtraction.result()->coefs();
//    gsDebugVar(basis);
//    gsDebugVar(newCoefs);

    gsBSpline<> bezierExtractionCurve( basis, newCoefs);
    gsWriteParaview(bezierExtractionCurve, "bezierC1", 1000);



    gsKnotVector<> kv(0, 1, 0, 2);//start,end,interior knots, start/end multiplicites of knots
    gsMatrix<> coefs2(2, 2);

    // coefs2 << 0.2, 0.2,
    //          0.6, 0.6;
    coefs2 << -10, -10,
             20, 20;

   gsBSpline<> curve2( kv, give(coefs2));

   gsInfo << "curve:\n" << curve2 << "\n";

//   gsWriteParaview( curve2, "line", 1000);
   gsWriteParaview( curve2, "line", 1000, false, true);

   // return 0;
   //find bounding box for bezier curves segments
   gsInfo << "----------------------------------------------------\n";
   gsInfo << deg <<"\n";
   gsKnotVector<> auxKnt(0, 1, 0, deg+1);//start,end,interior knots, start/end multiplicites of knots
   gsBSplineBasis<> basisAux(auxKnt);

   gsInfo << basisAux << "\n";
   gsMatrix<> coefsPatch(deg+1,2);
   gsInfo << "# Bezier patches =  " << basis.numElements() << "\n";
   gsInfo << newCoefs.rows() << " x " << newCoefs.cols() << "\n";

//   auto numElems = basis.numElements();
//   gsDebugVar(basis.numElements());

   gsMultiPatch<> mpCrv1, mpCrv2, markedCrv1;
   gsMatrix<> bbCrv1, bbCrv2;
   mpCrv2.addPatch(curve2);
   mpCrv2.boundingBox(bbCrv2);
   for(index_t i = 0; i < basis.numElements(); i++)
   {
     // gsInfo << newCoefs.middleRows(i*(deg), deg+1) << "\n\n";
     coefsPatch = newCoefs.middleRows(i*(deg), deg+1);
     gsBSpline<> bezierPatch(basisAux, coefsPatch);

     mpCrv1.clear();
     mpCrv1.addPatch(bezierPatch);
     mpCrv1.boundingBox(bbCrv1);
//     gsDebugVar(bbCrv1);

//     gsDebugVar(areBoundingBoxesIntersecting(bbCrv1, bbCrv2));
      if (areBoundingBoxesIntersecting(bbCrv1, bbCrv2)) {
        markedCrv1.addPatch(bezierPatch);
      }
   }

//  gsDebugVar(markedCrv1.nPatches());
//  gsDebugVar(bbCrv2);

  gsWriteParaview(markedCrv1, "c1", 1000, false, true);

   gsMultiPatch<> oneCrv = markedCrv1.patch(2);
  gsWriteParaview(oneCrv, "c1", 1000, false, true);

    // compute the intersection between oneCrv and curve2
//  gsDebugVar( oneCrv.coefs() );
//  gsDebugVar( curve2.coefs() );
//
//  gsDebugVar( oneCrv.basis(0) );
//  gsDebugVar( curve2.basis() );

  gsTensorBSplineBasis<2, real_t> basisBezierTensorProduct(auxKnt, kv);
//  gsDebugVar(basisBezierTensorProduct);
  gsMatrix<> coefsBezier(oneCrv.coefsSize()*mpCrv2.coefsSize(), 2);
  for (auto i=0; i!=oneCrv.coefsSize(); ++i) {
    for (auto j = 0; j != mpCrv2.coefsSize(); ++j) {
//      auto idxCurrPt = i * mpCrv2.coefsSize() + j;
      auto idxCurrPt = j*oneCrv.coefsSize() + i;
//      gsDebugVar(idxCurrPt);

      coefsBezier.row(idxCurrPt) = oneCrv.patch(0).coef(i) - mpCrv2.patch(0).coef(j);
    }
  }
//  gsDebugVar(coefsBezier);
  gsTensorBSpline<2,real_t> bezierTensorProduct(basisBezierTensorProduct,
                                                give(coefsBezier));

//  gsWriteParaview(bezierTensorProduct, "c1", 5000, false, false);

  // std::cout << "Testing the Booth function..." << std::endl;
  gsInfo << "Initial guess: \n";
//  gsEigen::VectorXd zInit(2); zInit << 0.5, 0.5;
  gsEigen::VectorXd zInit(2); zInit << 0.1, 0.1;
  gsMatrix<> xx = bezierTensorProduct.eval(zInit);
  gsWriteParaviewPoints(xx, "initialGuess");
  gsDebugVar(xx);

  ErrorFunction interSec = [&bezierTensorProduct](const Vector &z, Vector &fvec, Matrix &jac) {
//    double x = z(0);   double y = z(1);
    fvec.resize(2);
//
//    fvec(0) = x + 2*y - 7;
//    fvec(1) = 2*x + y - 5;
//  fvec(0) = x * x + y - 11;
//  fvec(1) = x + y * y - 7;

// TODO: fix this parameter transformation
//    gsVector<real_t> param = z;
//    param(0) = 1 / (1+ exp(-z(0)));
//    param(1) = 1 / (1+ exp(-z(1)));
////////////////////////////////////////////////////////////////////////////////

    fvec = bezierTensorProduct.eval(z);
    return nullptr;
  };

  // TODO: next step is to compute intersection using PP algorithm
//  gsMatrix<> ptCrv = curve1.eval(curve1.basis().knots(0).greville());
//  gsVector<> pt, res;
//
//  for (auto i = 0; i != ptCrv.cols(); ++i) {
//    pt = ptCrv.col(0);
//    real_t dist = curve2.closestPointTo(pt, res);
//    gsDebugVar(dist);
//  }

//  std::cout << "Testing the Booth function..." << std::endl;
//  gsEigen::VectorXd zInit(2); zInit << 1.87, 2.032;
//  std::cout << "zInit: " << zInit.transpose() << std::endl;
//  gsEigen::VectorXd zSoln(2); zSoln << 1.0, 3.0;
//  std::cout << "zSoln: " << zSoln.transpose() << std::endl;

//  BoothFunctor functor;
//  gsEigen::NumericalDiff<BoothFunctor> numDiff(functor);
//  lsq::LevenbergMarquardt<> lm(numDiff);
//  lm.parameters.maxfev = 1000;
//  lm.parameters.xtol = 1.0e-10;

//  lsq::LevenbergMarquardt lm();
//  lm().minimize(zInit);

  // Create LevenbergMarquardt object

  lsq::LevenbergMarquardt<real_t, ErrorFunction, lsq::NoCallback<real_t>,
      lsq::CentralDifferences<real_t>> lmOptimizer;
//  lsq::DoglegMethod<real_t, ErrorFunction, lsq::NoCallback<real_t>,
//                          lsq::CentralDifferences<real_t>> lmOptimizer;

//  lmOptimizer().setErrorFunction(BoothFunctor);
  // Configure the optimizer if necessary
//  lmOptimizer.setIncrease(2.0);
//  lmOptimizer.setDecrease(0.5);
//  lmOptimizer.setLambda(1.0);
//  lmOptimizer.setMaxItLM(100);

//  gsVector<> fval;
//  gsMatrix<> jac;
//  BoothFunctor(zInit, fval, jac);

//  lmOptimizer.setErrorFunction(BoothFunctor);
  lmOptimizer.setErrorFunction(interSec);
  lmOptimizer.setVerbosity(3);

  // Run the minimization
  auto res = lmOptimizer.minimize(zInit);

//  gsDebugVar( bezierTensorProduct.eval(res.xval) );

  auto param = res.xval;
  gsVector<> param1(1); param1 << param(0);
  gsVector<> param2(1); param2 << param(1);
  gsDebugVar(oneCrv.patch(0).eval(param1));
  gsDebugVar(mpCrv2.patch(0).eval(param2));

  // zInit
  gsWriteParaviewPoints(oneCrv.patch(0).eval(zInit.transpose()), "guessi_1");
  gsWriteParaviewPoints(mpCrv2.patch(0).eval(zInit.transpose()), "guessi_2");
  gsWriteParaviewPoints(oneCrv.patch(0).eval(param1), "guessf_1");
  gsWriteParaviewPoints(mpCrv2.patch(0).eval(param2), "guessf_2");

//  std::cout << "max fun eval: " << lm.parameters.maxfev << std::endl;
//  std::cout << "x tol: " << lm.parameters.xtol << std::endl;

  gsInfo << "iter count: " << res.iterations << "\n";
  gsInfo << "error: " << res.error << "\n";
  gsInfo << "return status: " << res.converged << "\n";
  gsInfo << "zSolver: " << res.xval.transpose() << "\n";
  gsInfo << "fVal: " << res.fval.transpose() << "\n";
  gsInfo << "-----------------------------------------------------------\n";

  return 0;
}
