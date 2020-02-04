/** @file biharmonic_twoPatch.cpp

    @brief A Biharmonic example for ONLY TWO-PATCHES

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): P. Weinmueller
*/

# include <gismo.h>
# include <gsAssembler/gsG1BiharmonicAssembler.h>

#include <gsG1Basis/gsG1Basis.h>
#include <gsG1Basis/gsG1System.h>

#include <gsG1Basis/gsNormL2.h>
#include <gsG1Basis/gsSeminormH1.h>
#include <gsG1Basis/gsSeminormH2.h>
#include <gsG1Basis/gsH1NormWithJump.h>

using namespace gismo;

int main(int argc, char *argv[])
{
    index_t numRefine = 4;
    index_t numDegree = 0;
    index_t regularity = 1;

    // For the spline space of the gluing data
    index_t p_tilde = 1;
    index_t r_tilde = 0;

    bool plot = false;
    index_t geometry = 0;
    bool direct = false;
    bool local = false;

    gsCmdLine cmd("Example for solving the biharmonic problem.");
    cmd.addInt("k", "refine", "Number of refinement steps", numRefine);
    cmd.addInt("p", "p_tilde", "Polynomial degree for tilde{p}", p_tilde);
    cmd.addInt("r", "r_tilde", "Regularity for tilde{r}", r_tilde);
    cmd.addSwitch( "plot", "Plot result in ParaView format", plot );
    cmd.addSwitch( "direct", "Construction of the G1 basis functions", direct );
    cmd.addSwitch( "local", "Gluing data local or global", local );
    cmd.addInt("g", "geometry", "Geometry", geometry);
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    // ======= Solution =========
    gsFunctionExpr<> source  ("256*pi*pi*pi*pi*(4*cos(4*pi*x)*cos(4*pi*y) - cos(4*pi*x) - cos(4*pi*y))",2);
    gsFunctionExpr<> laplace ("-16*pi*pi*(2*cos(4*pi*x)*cos(4*pi*y) - cos(4*pi*x) - cos(4*pi*y))",2);
    gsFunctionExpr<> solVal("(cos(4*pi*x) - 1) * (cos(4*pi*y) - 1)",2);
    gsFunctionExpr<>sol1der ("-4*pi*(cos(4*pi*y) - 1)*sin(4*pi*x)",
                              "-4*pi*(cos(4*pi*x) - 1)*sin(4*pi*y)",2);
    gsFunctionExpr<>sol2der ("-16*pi^2*(cos(4*pi*y) - 1)*cos(4*pi*x)",
                              "-16*pi^2*(cos(4*pi*x) - 1)*cos(4*pi*y)",
                              " 16*pi^2*sin(4*pi*x)*sin(4*pi*y)", 2);
    gsFunctionWithDerivatives<real_t> solution(solVal, sol1der, sol2der);

    // ======= Geometry =========
    std::string string_geo;
    bool oneComponent = false;
    switch(geometry)
    {
        case 0:
            string_geo = "planar/twoPatches/square_diagonal.xml";
            numDegree = 2; // 2 == degree 3
            break;
        case 1:
            string_geo = "planar/twoPatches/square_curved.xml";
            numDegree = 0; // 0 == degree 3
            break;
        case 2:
            string_geo = "planar/twoPatches/square_cuttedCorner.xml";
            numDegree = 0; // 0 == degree 3
            oneComponent = true;
            break;
        default:
            gsInfo << "No geometry is used! \n";
            break;
    }

    gsFileData<> fd(string_geo);
    gsInfo << "Loaded file "<< fd.lastPath() <<"\n";

    gsMultiPatch<> multiPatch;
    fd.getId(0, multiPatch); // id=0: Multipatch domain
    multiPatch.computeTopology();

// FOR LOOP:

//gsVector<index_t> level_pTilde(4);
//level_pTilde << 1, 2, 3, 7;
//for (index_t level_pt = 0; level_pt < level_pTilde.size(); level_pt++)
//{
//p_tilde = level_pTilde(level_pt);
//gsInfo << " ####### p_tilde " << p_tilde << " ######## \n";

//gsVector<index_t> level_refine(5);
//level_refine << 4, 9, 19, 39, 79;
//for (index_t level = 0; level < level_refine.size(); level++)
//{
    //numRefine = level_refine(level);
    //gsInfo << " ####### numRefine " << numRefine << " ######## \n";


    gsMultiBasis<> multiBasis(multiPatch);
    if (oneComponent)
        multiBasis.degreeElevate(1,0);
    multiBasis.degreeElevate(numDegree);
    //multiBasis.degreeElevate(1,1);

    index_t polynomDegree = multiBasis.minCwiseDegree();
    multiBasis.uniformRefine(numRefine, polynomDegree - regularity);
    // ==========================

/*    gsMultiPatch<> mp;
    mp = multiPatch;
    mp.uniformRefine(numRefine, polynomDegree - regularity);
    gsWriteParaview(mp,"geometry",5000,true);
*/
    // ==========================
    //gsInfo << "Basis at the interface is computing... \n";
    gsG1Basis<real_t> g1Basis(multiPatch,multiBasis,numRefine,regularity,
        p_tilde,r_tilde,direct,local,plot);
    g1Basis.assemble();
    g1Basis.solve();

    gsMultiPatch<real_t> basisG1_L, basisG1_R;
    g1Basis.constructSolution(basisG1_L,basisG1_R);

    if (plot)
        g1Basis.plotG1Basis(basisG1_L,basisG1_R);


    // ======= Boundary =========
    gsBoundaryConditions<> bcInfo, bcInfo2;
    for (gsMultiPatch<>::const_biterator bit = multiPatch.bBegin(); bit != multiPatch.bEnd(); ++bit)
    {
        bcInfo.addCondition( *bit, condition_type::dirichlet, &solVal ); // = 0
        bcInfo2.addCondition(*bit, condition_type::neumann, &laplace ); // = 0
    }
    // BiharmonicAssembler
    gsG1BiharmonicAssembler<real_t> g1BiharmonicAssembler(multiPatch, multiBasis, bcInfo, bcInfo2, source);

    g1BiharmonicAssembler.assemble();
    g1BiharmonicAssembler.computeDirichletDofsL2Proj(basisG1_L,basisG1_R,g1Basis.get_n_tilde(),
                                                        g1Basis.get_n_bar());

    //gsInfo << "Build system ... \n";
    gsG1System<real_t> g1System(basisG1_L, basisG1_R,
                                g1BiharmonicAssembler.get_g1dofs(),
                                g1BiharmonicAssembler.get_mapper(),
                                g1BiharmonicAssembler.get_mapper_boundary(),
                                g1BiharmonicAssembler.get_mapper_interface(),
                                g1BiharmonicAssembler.matrix().dim().first,
                                g1Basis.get_n_tilde(),
                                g1Basis.get_n_bar());
    g1System.assemble();


    // Solving system:
    gsMatrix<real_t> f = g1BiharmonicAssembler.rhs(); // with the second boundary condition
    gsSparseMatrix<real_t> K_sparse = g1BiharmonicAssembler.matrix();


    gsSparseMatrix<real_t> D_0_sparse = g1System.get_D_0_sparse();
    gsSparseMatrix<real_t> D_boundary_sparse = g1System.get_D_boundary_sparse();

    gsMatrix<real_t> g = g1System.get_g();

    //gsInfo << "Solving system... \n";
    gsSparseMatrix<real_t> A = D_0_sparse * K_sparse * D_0_sparse.transpose();
    gsVector<real_t> F = D_0_sparse * f - D_0_sparse * K_sparse * D_boundary_sparse.transpose() * g;

    //gsInfo << "System finished with " << A.nonZeros() << " non-zeros!\n";

    gsSparseSolver<real_t>::CGDiagonal solver;
    //gsSparseSolver<real_t>::LU solver;
    //solver.analyzePattern(BiharmonicAssembler.matrix() );
    //solver.factorize(BiharmonicAssembler.matrix());
    //gsInfo << "matrix: " << A.dim() << "\n";
    //gsInfo << "rhs: " << F.dim() << "\n";
    solver.compute(A);
    gsMatrix<> solVector= solver.solve(F);

    gsMultiPatch<> mpsol;
    g1BiharmonicAssembler.constructSolution(solVector.bottomRows(g1BiharmonicAssembler.matrix().dim().first),mpsol);
    gsField<> solField(multiPatch, mpsol);

    g1System.constructSolution_G1(solVector.topRows(basisG1_L.nPatches()),
                                  basisG1_L, basisG1_R);

    if (plot)
    {
        //gsInfo<<"Plotting in Paraview...\n";
        //gsWriteParaview(solField,"biharmonic_trafo",5000);

        const gsField<> exact( multiPatch, solution, false );
        gsWriteParaview<>( exact, "Biharmonic2d_exact", 5000);

        g1BiharmonicAssembler.writeParaview(solField,"biharmonic_trafo_g1",basisG1_L,basisG1_R,5000);

    }

    std::vector< gsMultiPatch<>> basisG1;
    basisG1.resize(2);
    basisG1.at(0) = basisG1_L;
    basisG1.at(1) = basisG1_R;

    gsNormL2<real_t> error(solField,solVal,basisG1);
    error.compute();

    //gsInfo << "L2 error : "<< error.value() << "\n";

    //gsInfo << "L2 error: H1 Seminorm error: H2 Seminorm error: H2 Norm error: errorJump error: \n";
    gsInfo << error.value();

    gsSeminormH1<real_t> errorH1(solField,solVal,basisG1);
    errorH1.compute();

    gsInfo << " " << errorH1.value();

    gsSeminormH2<real_t> errorH2(solField,solVal,basisG1);
    errorH2.compute();

    gsInfo << " " << errorH2.value();

    gsInfo << " " << math::sqrt(errorH2.value()*errorH2.value() +
        errorH1.value()*errorH1.value() + error.value()*error.value());

    boundaryInterface  iFace = *multiPatch.iBegin();
    gsH1NormWithJump<real_t> errorJump(solField,basisG1,iFace);
    errorJump.compute();
    gsInfo << " " << errorJump.value() << "\n";


//} // level
//} // level_pt
    return  0;

}