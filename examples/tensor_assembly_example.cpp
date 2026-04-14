/** @file assembly_example.cpp

    @brief Tutorial on how to use the gsAssembler class

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris
*/

# include <gismo.h>

# include <gsTensor/gsTensorFunction.h>
# include <gsTensor/gsTensorAssembler.h>


using namespace gismo;

int main(int argc, char *argv[])
{
    short_t D = 2;
    index_t numRefine = 0,
        numElevate = 0,
        degree = 2,
        knots = 0;
    
    gsCmdLine cmd("Tutorial on assemblying a Poisson problem.");
    cmd.addInt("d", "dimension", "Dimension 2d/3d", D);
    cmd.addInt("p", "degree", "Degree of the initial tensor product basis", degree);
    cmd.addInt("k", "knots", "Number of knots of the initial tensor product basis", knots);
    cmd.addInt("u", "uniform", "Number of uniform refinement", numRefine);
    cmd.addInt("e", "elevate", "Number of degree elevation steps", numElevate);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }
    //! [Parse command line]
    
    gsTensorFunction<> RR(2,1); // creates a constant function = 1
    // gsDebugVar(RR);

    std::vector<gsKnotVector<> > knot_vector(D, gsKnotVector<>(0, 1, knots, degree + 1));
    gsBasis<>::uPtr basis = gsBSplineBasis<>::create(knot_vector);

    if (numElevate > 0)
        basis->degreeElevate(numElevate);
    for (int i = 0; i < numRefine; ++i)
        basis->uniformRefine();

    gsTensorAssembler<> ta;
    ta.compute(*basis, STIFFNESS, RR); // DIFFUSION ?
    // ta.compute(*basis, MASS, RR); // DIFFUSION ? bug for this line?

    gsFiberMatrix<real_t> mat = ta.kronecker().toFiberMatrix();

    gsInfo << "Matrix size: \n" << mat.rows() <<"\n";
    gsInfo << "Matrix : \n" << mat.toSparseMatrix().toDense() <<"\n";


    // comparison with standard assembly
    gsExprAssembler<> A(1, 1);
    std::string fn = "planar/unitsquare.xml";
    gsFileData<> fd(fn);
    gsInfo << "Loaded file "<< fd.lastPath() <<"\n";

    gsMultiPatch<real_t>::uPtr mp = gsReadFile<>(fn);
    gsInfo << " Got" << *mp << " \n";
    gsExprAssembler<>::geometryMap G = A.getMap(*mp);
    A.setIntegrationDomain(basis->domain());
    auto u = A.getSpace(*basis);
    auto v = A.getTestSpace(*basis, 1);
    auto w = A.getCoeff(RR);
    A.initMatrix();
    A.assemble( igrad(u, G) * igrad(u, G).tr() * meas(G) );

    gsMatrix<> mat_std = A.giveMatrix().toDense();
    gsInfo << "Standard assembly matrix: \n" << mat_std <<"\n";
    gsInfo << "Difference: \n" << (mat.toSparseMatrix().toDense() - mat_std).cwiseAbs().maxCoeff() <<"\n";


    // diffusion assembly test
    gsTensorFunction<> density_field(2,0); // creates a constant function = 1
    gsTensorAssembler<> ta_diffusion;
    ta_diffusion.compute(*basis, DIFFUSION, density_field); // DIFFUSION

    gsFiberMatrix<real_t> mat_diffusion = ta_diffusion.kronecker().toFiberMatrix();

    gsInfo << "Matrix size: \n" << mat_diffusion.rows() <<"\n";
    gsInfo << "Matrix : \n" << mat_diffusion.toSparseMatrix().toDense() <<"\n";

    return 0;
}

