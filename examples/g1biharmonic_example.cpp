/** @file biharmonic_multiPatch.cpp

    @brief A Biharmonic example for ONLY TWO-PATCHES

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): P. Weinmueller
*/
# include <gismo.h>
# include <gsG1Basis/gsG1AuxiliaryEdgeMultiplePatches.h>
# include <gsG1Basis/gsG1BasisVertex.h>
# include <gsAssembler/gsG1BiharmonicAssembler.h>

using namespace gismo;

int main(int argc, char *argv[])
{
    // Geometry data
    index_t geometry = 0; // Which geometry

    index_t numRefine = 4;
    index_t numDegree = 0;
    index_t regularity = 1;

    // For the spline space of the gluing data
    index_t p_tilde = 1;
    index_t r_tilde = 0;

    index_t threads = 1;

    bool plot = false;
    bool direct = false;
    bool local = false;
    bool loop = false;
    bool local_g1 = false;

    gsCmdLine cmd("Example for solving the biharmonic problem.");
    cmd.addInt("k", "refine", "Number of refinement steps", numRefine);
    cmd.addInt("p", "p_tilde", "Polynomial degree for tilde{p}", p_tilde);
    cmd.addInt("r", "r_tilde", "Regularity for tilde{r}", r_tilde);
    cmd.addSwitch( "plot", "Plot result in ParaView format", plot );
    cmd.addSwitch( "direct", "Construction of the G1 basis functions", direct );
    cmd.addSwitch( "loop", "If you want to solve several levels", loop );
    cmd.addSwitch( "local_g1", "If you want to solve several levels", local_g1 );
    cmd.addInt("g", "geometry", "Geometry", geometry);
    cmd.addInt("t", "threads", "Threads", threads);
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
    switch(geometry)
    {
        case 0:
            string_geo = "planar/multiPatches/4_square_diagonal.xml";
            numDegree = 2; // 2 == degree 3
            break;
        case 1:
            string_geo = "planar/multiPatches/6_square_diagonal.xml";
            numDegree = 2; // 2 == degree 3
            break;
        case 2:
            string_geo = "planar/multiPatches/square_curved.xml";
            numDegree = 0; // 0 == degree 3
            break;
        case 3:
            string_geo = "planar/twoPatches/square_curved.xml";
            numDegree = 0; // 0 == degree 3
            break;
        case 4:
            string_geo = "planar/twoPatches/square_curved_deg_5.xml";
            numDegree = 0; // 0 == degree 3
            break;
        case 5:
            string_geo = "planar/twoPatches/square_curved_deg_7.xml";
            numDegree = 0; // 0 == degree 3
            break;
        case 6:
            string_geo = "planar/twoPatches/square_non_conform.xml";
            numDegree = 2; // 0 == degree 3
            break;
        case 7:
            string_geo = "planar/twoPatches/square_bent.xml";
            numDegree = 0; // 0 == degree 3
            break;
        case 8:
            string_geo = "planar/twoPatches/square_complex_bent.xml";
            numDegree = 0; // 0 == degree 3
            break;
        case 9:
            string_geo = "planar/twoPatches/square_diagonal.xml";
            numDegree = 2; // 2 == degree 3
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


    gsWriteParaview(multiPatch,"geometry",5000,true);

    gsOptionList optionList;
    optionList.addInt("p_tilde","Grad",p_tilde);
    optionList.addInt("r_tilde","Reg",r_tilde);
    optionList.addInt("regularity","Regularity of the initial geometry",regularity);
    optionList.addSwitch("local","Local projection for gluing data",local);
    optionList.addSwitch("direct","Local projection for gluing data",direct);
    optionList.addSwitch("plot","Plot in Paraview",plot);
    optionList.addInt("refine","Refinement",numRefine);
    optionList.addInt("degree","Degree",numDegree);

    //multiPatch.patch(1).degreeElevate(1,0);
    multiPatch.degreeElevate(optionList.getInt("degree"));

    multiPatch.uniformRefine_withSameRegularity(optionList.getInt("refine"), optionList.getInt("regularity"));
    gsMultiBasis<> mb(multiPatch);

    gsInfo << "basis : " << mb.basis(0) << "\n";
    gsInfo << "basis : " << mb.basis(1) << "\n";

/*
    // Interface loop
    for (const boundaryInterface &  item : multiPatch.interfaces() )
    {
        gsG1AuxiliaryEdgeMultiplePatches a(multiPatch, item.first().patch, item.second().patch);
        a.computeG1InterfaceBasis(optionList);
    }

    for (gsMultiPatch<>::const_biterator bit = multiPatch.bBegin(); bit != multiPatch.bEnd(); ++bit)
    {
        gsInfo << "Patch: " << bit->patch << "\n";
        gsInfo << "m_index: " << bit->m_index << "\n";
        gsG1AuxiliaryEdgeMultiplePatches a(multiPatch, bit->patch);
        a.computeG1BoundaryBasis(optionList, bit->m_index);
    }
*/

    gsBSplineBasis<> bsp_temp = dynamic_cast<gsBSplineBasis<> & >(mb.basis(0).component(0)); // TODO for all components
    real_t p = bsp_temp.maxDegree(); // p_max
    real_t h_geo = bsp_temp.knots().at(p + 2); // h_max (?)
    real_t val = 2; // TODO Change

    gsMatrix<> zero;
    zero.setZero(2,1);
    real_t sigma = 0.0;
    for (index_t i = 0; i < val; i++)
        sigma += multiPatch.patch(i).deriv(zero).lpNorm<Eigen::Infinity>();
    sigma *= h_geo/(val*p);
    sigma = 1 / sigma;

    gsInfo << "SIGMA : " << sigma << "\n";
    optionList.addReal("sigma","Sigma for the Vertex",sigma);
    std::vector<bool> isBoundary_0, isBoundary_1;
    isBoundary_0.push_back(false); // interface at u
    isBoundary_0.push_back(false); // boundary at v
    isBoundary_1.push_back(false); // boundary at u
    isBoundary_1.push_back(false); // interface at v

    gsG1BasisVertex<real_t> g1BasisVertex_0(multiPatch.patch(0),mb.basis(0),isBoundary_0,optionList);
    gsG1BasisVertex<real_t> g1BasisVertex_1(multiPatch.patch(1),mb.basis(1),isBoundary_1,optionList);
    gsMultiPatch<> g1Basis_0, g1Basis_1;
    g1BasisVertex_1.constructSolution(g1Basis_1);
    g1BasisVertex_1.plotG1BasisBoundary(g1Basis_1,multiPatch.patch(1),"BasisVertex");
    g1BasisVertex_0.constructSolution(g1Basis_0);
    g1BasisVertex_0.plotG1BasisBoundary(g1Basis_0,multiPatch.patch(0),"BasisVertex0");

// NEW NEW NEW NEW NEW NEW NEW NEW NEW

    gsBoundaryConditions<> bcInfo, bcInfo2;
    for (gsMultiPatch<>::const_biterator bit = multiPatch.bBegin(); bit != multiPatch.bEnd(); ++bit)
    {
        bcInfo.addCondition( *bit, condition_type::dirichlet, &solVal ); // = 0
        bcInfo2.addCondition(*bit, condition_type::neumann, &laplace ); // = 0
    }

    // BiharmonicAssembler
    gsG1BiharmonicAssembler<real_t> g1BiharmonicAssembler(multiPatch, mb, bcInfo, bcInfo2, source);
    g1BiharmonicAssembler.assemble();

    // TODO g1BiharmonicAssembler.computeDirichletDofsL2Proj(basisG1, n_tilde, n_bar );

} // main