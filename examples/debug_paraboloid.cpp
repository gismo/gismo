#include <gismo.h>
using namespace gismo;
int main() {
    gsKnotVector<> kv(0,1,0,3);
    gsTensorBSplineBasis<2> basis(kv,kv);
    gsMatrix<> coefs(9,3);
    coefs << 0,0,0.5, 0.5,0,1, 1,0,0.5,
             0,0.5,1, 0.5,0.5,1.5, 1,0.5,1,
             0,1,0.5, 0.5,1,1, 1,1,0.5;
    gsTensorBSpline<2> patch(basis, coefs);

    gsMatrix<> pt(2,1); pt << 0.25, 0.25;
    gsMapData<real_t> md;
    md.flags = NEED_GRAD_TRANSFORM | NEED_MEASURE;
    md.points = pt;
    patch.computeMap(md);

    gsInfo << "jacInvTr (raw):\n" << md.jacInvTr.transpose() << "\n";
    gsMatrix<> raw = md.jacInvTr.col(0); raw.resize(3,2);
    gsInfo << "Stored as 3x2 (= J*(J^TJ)^-1):\n" << raw << "\n";
    gsInfo << "Transposed 2x3 (= (J^TJ)^-1 J^T, left pseudo-inv):\n" << raw.transpose() << "\n";
    gsInfo << "Expected: [[2/3,-1/3,1/3],[-1/3,2/3,1/3]]\n";
    gsInfo << "i.e.:     [[" << 2./3 << "," << -1./3 << "," << 1./3 << "],["
           << -1./3 << "," << 2./3 << "," << 1./3 << "]]\n";
    gsInfo << "measure = " << md.measures << "\n";
    gsInfo << "expected meas = sqrt(3) = " << std::sqrt(3.) << "\n";
    return 0;
}
