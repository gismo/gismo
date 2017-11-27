/* @file newton_check.cpp

   Tests the newton iteration
*/

//#define TEST_INFO


#include "gismo_unittest.h" // Brings in G+Smo and the UnitTest++ framework

TEST(newton_check) // Declares test
{
    UNITTEST_TIME_CONSTRAINT(1000);// this will produce failure if test takes more than 1 sec

    gsGeometry<>::Ptr f = gsNurbsCreator<>::NurbsQuarterAnnulus();
    gsVector<> x = gsVector<>::vec(0.5, 0.5),
               y = gsVector<>::vec(0.1, 1.9);

    // For any function
    int iter = f->newtonRaphson(y, x, true);

    gsMatrix<> fx = f->eval(x);

    gsTestInfo << "Result:     " << x .transpose() << "\n";
    gsTestInfo << "Value:      " << fx.transpose() << "\n";
    real_t res = (fx - y).norm();
    gsTestInfo << "Res.norm:   " <<  res << "\n";
    gsTestInfo << "Iterations: " << iter << "\n";

    CHECK( iter >= 0   );
    CHECK( res <= 1e-5 );

    // For a gsGeometry
    gsTestInfo << "-- Invert point on geometry.\n";
    gsMatrix<> points(2,2);
    points.col(0) = gsVector<>::vec(0.2, 1.7);
    points.col(1) = gsVector<>::vec(2.0, 0.0);

    gsMatrix<> params;
    f->invertPoints(points, params);

    fx = f->eval(params);

    gsTestInfo << "Result:     " << params.asRowVector() << "\n";
    gsTestInfo << "Value:      " << fx    .asRowVector() << "\n";
    res = (fx - points).norm();
    gsTestInfo << "Res.norm:   " << res << "\n";

    CHECK( res <= 1e-5 );
}
