#include <gismo.h>

using namespace gismo;
//! [Include namespace]

int main(int argc, char *argv[])
{
    gsTensorNurbs<2> geom = * gsNurbsCreator<>::NurbsQuarterAnnulus(1,2);

    gsVector<> query(2), res(2);
    // (-3, 3) breaks both versions
    query << -3., 3.;
    res.setZero();
    real_t  distance = geom.closestPointTo(query, res);
    auto resXY = geom.eval( res );
    gsInfo << "Closest Point: " << resXY.transpose() << "\n";

    gsWriteParaview(geom, "NewtRaphTest");
    gsMatrix<> points(2,2);
    points.col(0) = query;
    points.col(1) = resXY;
    gsWriteParaviewPoints(points, "NewtRaphTestPoints");
    return 0;
}