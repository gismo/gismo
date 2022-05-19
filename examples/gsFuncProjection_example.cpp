/** @file biharmonic2_example.cpp

    @brief Tutorial on how to use expression assembler and the (approx.) C1 basis function
                to solve the Biharmonic equation

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): P. Weinmueller
*/

//! [Include namespace]
#include <gismo.h>

#include <gsUnstructuredSplines/src/gsApproxC1Spline.h>
#include <gsUnstructuredSplines/src/gsDPatch.h>
#include <gsUnstructuredSplines/src/gsAlmostC1.h>
#include <gsUnstructuredSplines/src/gsC1SurfSpline.h>

#include <gsUtils/gsQuasiInterpolate.h>
#include <gsUtils/gsFuncProjection.h>

using namespace gismo;
//! [Include namespace]

/**
 * Smoothing method:
 * - m 0 == Approx C1 method
 * - m 1 == D-Patch method
 * - m 2 == Almost C1 method
 * - m 3 == Nitsche's method
 */
enum MethodFlags
{
    APPROXC1       = 0, // Approx C1 Method
    DPATCH         = 1, // D-Patch
    ALMOSTC1       = 2, // Almost C1
    NITSCHE        = 3, // Nitsche
    SURFASG1       = 5, // AS-G1
    // Add more [...]
};

void setMapperForBiharmonic(gsBoundaryConditions<> & bc, gsMappedBasis<2,real_t> & bb2, gsDofMapper & mapper)
{
    mapper.setIdentity(bb2.nPatches(), bb2.size(), 1);

    gsMatrix<index_t> bnd;
    for (typename gsBoundaryConditions<real_t>::const_iterator
                 it = bc.begin("Dirichlet"); it != bc.end("Dirichlet"); ++it)
    {
        bnd = bb2.basis(it->ps.patch).boundary(it->ps.side());
        mapper.markBoundary(it->ps.patch, bnd, 0);
    }

    for (typename gsBoundaryConditions<real_t>::const_iterator
             it = bc.begin("Neumann"); it != bc.end("Neumann"); ++it)
    {
        bnd = bb2.basis(it->ps.patch).boundaryOffset(it->ps.side(),1);
        mapper.markBoundary(it->ps.patch, bnd, 0);
    }
    mapper.finalize();
}

void setMapperForBiharmonic(gsBoundaryConditions<> & bc, gsMultiBasis<> & dbasis, gsDofMapper & mapper)
{
    mapper.init(dbasis);

    for (gsBoxTopology::const_iiterator it = dbasis.topology().iBegin();
         it != dbasis.topology().iEnd(); ++it) // C^0 at the interface
    {
        dbasis.matchInterface(*it, mapper);
    }

    gsMatrix<index_t> bnd;
    for (typename gsBoundaryConditions<real_t>::const_iterator
                 it = bc.begin("Dirichlet"); it != bc.end("Dirichlet"); ++it)
    {
        bnd = dbasis.basis(it->ps.patch).boundary(it->ps.side());
        mapper.markBoundary(it->ps.patch, bnd, 0);
    }

    for (typename gsBoundaryConditions<real_t>::const_iterator
                 it = bc.begin("Neumann"); it != bc.end("Neumann"); ++it)
    {
        bnd = dbasis.basis(it->ps.patch).boundaryOffset(it->ps.side(),1);
        mapper.markBoundary(it->ps.patch, bnd, 0);
    }
    mapper.finalize();
}

void gsDirichletNeumannValuesL2Projection(gsMultiPatch<> & mp, gsMultiBasis<> & dbasis, gsBoundaryConditions<> & bc,
                                           gsMappedBasis<2,real_t> & bb2, const expr::gsFeSpace<real_t> & u)
{
    const gsDofMapper & mapper = u.mapper();

    gsMatrix<index_t> bnd = mapper.findFree(mapper.numPatches()-1);
    gsDofMapper mapperBdy;
    mapperBdy.setIdentity(bb2.nPatches(), bb2.size(), 1);  // bb2.nPatches() == 1
    mapperBdy.markBoundary(0, bnd, 0);
    mapperBdy.finalize();

    gsExprAssembler<real_t> A(1,1);
    A.setIntegrationElements(dbasis);

    auto G = A.getMap(mp);
    auto uu = A.getSpace(bb2);
    auto g_bdy = A.getBdrFunction(G);

    uu.setupMapper(mapperBdy);
    gsMatrix<real_t> & fixedDofs_A = const_cast<expr::gsFeSpace<real_t>&>(uu).fixedPart();
    fixedDofs_A.setZero( uu.mapper().boundarySize(), 1 );

    real_t lambda = 1e-5;

    A.initSystem();
    A.assembleBdr(bc.get("Dirichlet"), uu * uu.tr() * meas(G));
    A.assembleBdr(bc.get("Dirichlet"), uu * g_bdy * meas(G));
    A.assembleBdr(bc.get("Neumann"),
                  lambda * (igrad(uu, G) * nv(G).normalized()) * (igrad(uu, G) * nv(G).normalized()).tr() * meas(G));
    A.assembleBdr(bc.get("Neumann"),
                  lambda *  (igrad(uu, G) * nv(G).normalized()) * (g_bdy.tr() * nv(G).normalized()) * meas(G));

    gsSparseSolver<real_t>::SimplicialLDLT solver;
    solver.compute( A.matrix() );
    gsMatrix<real_t> & fixedDofs = const_cast<expr::gsFeSpace<real_t>& >(u).fixedPart();
    fixedDofs = solver.solve(A.rhs());
}

void gsDirichletNeumannValuesL2Projection(gsMultiPatch<> & mp, gsMultiBasis<> & dbasis, gsBoundaryConditions<> & bc, const expr::gsFeSpace<real_t> & u)
{
    gsDofMapper mapper = u.mapper();
    gsDofMapper mapperBdy(dbasis, u.dim());
    for (gsBoxTopology::const_iiterator it = dbasis.topology().iBegin();
         it != dbasis.topology().iEnd(); ++it) // C^0 at the interface
    {
        dbasis.matchInterface(*it, mapperBdy);
    }
    for (size_t np = 0; np < mp.nPatches(); np++)
    {
        gsMatrix<index_t> bnd = mapper.findFree(np);
        mapperBdy.markBoundary(np, bnd, 0);
    }
    mapperBdy.finalize();

    gsExprAssembler<real_t> A(1,1);
    A.setIntegrationElements(dbasis);

    auto G = A.getMap(mp);
    auto uu = A.getSpace(dbasis);
    auto g_bdy = A.getBdrFunction(G);

    uu.setupMapper(mapperBdy);
    gsMatrix<real_t> & fixedDofs_A = const_cast<expr::gsFeSpace<real_t>&>(uu).fixedPart();
    fixedDofs_A.setZero( uu.mapper().boundarySize(), 1 );

    real_t lambda = 1e-5;

    A.initSystem();
    A.assembleBdr(bc.get("Dirichlet"), uu * uu.tr() * meas(G));
    A.assembleBdr(bc.get("Dirichlet"), uu * g_bdy * meas(G));
    A.assembleBdr(bc.get("Neumann"),
                  lambda * (igrad(uu, G) * nv(G).normalized()) * (igrad(uu, G) * nv(G).normalized()).tr() * meas(G));
    A.assembleBdr(bc.get("Neumann"),
                  lambda *  (igrad(uu, G) * nv(G).normalized()) * (g_bdy.tr() * nv(G).normalized()) * meas(G));

    gsSparseSolver<real_t>::SimplicialLDLT solver;
    solver.compute( A.matrix() );
    gsMatrix<real_t> & fixedDofs = const_cast<expr::gsFeSpace<real_t>& >(u).fixedPart();
    gsMatrix<real_t> fixedDofs_temp = solver.solve(A.rhs());

    // Reordering the dofs of the boundary
    fixedDofs.setZero(mapper.boundarySize(),1);
    index_t sz = 0;
    for (size_t np = 0; np < mp.nPatches(); np++)
    {
        gsMatrix<index_t> bnd = mapperBdy.findFree(np);
        bnd.array() += sz;
        for (index_t i = 0; i < bnd.rows(); i++)
        {
            index_t ii = mapperBdy.asVector()(bnd(i,0));
            fixedDofs(mapper.global_to_bindex(mapper.asVector()(bnd(i,0))),0) = fixedDofs_temp(ii,0);
        }
        sz += mapperBdy.patchSize(np,0);
    }
}

void computeStabilityParameter(gsMultiPatch<> mp, gsMultiBasis<> dbasis, gsMatrix<real_t> & mu_interfaces)
{
    mu_interfaces.setZero();

    index_t i = 0;
    for ( typename gsMultiPatch<real_t>::const_iiterator it = mp.iBegin(); it != mp.iEnd(); ++it, ++i)
    {
        gsMultiPatch<> mp_temp;
        mp_temp.addPatch(mp.patch(it->first().patch));
        mp_temp.addPatch(mp.patch(it->second().patch));
        mp_temp.computeTopology();

        gsMultiBasis<> dbasis_temp;
        dbasis_temp.addBasis(dbasis.basis(it->first().patch).clone().release());
        dbasis_temp.addBasis(dbasis.basis(it->second().patch).clone().release());

        gsBoundaryConditions<> bc;

//        patchSide pS1 = mp_temp.interfaces()[0].first();
//        patchSide pS2 = mp_temp.interfaces()[0].second();
//
//
//        index_t side = pS1.index() < 3 ? (pS1.index() == 1 ? 2 : 1) : (pS1.index() == 3 ? 4 : 3);
//        bc.addCondition(patchSide(pS1.patchIndex(), side), condition_type::dirichlet, 0);
//
//        side = pS2.index() < 3 ? (pS2.index() == 1 ? 2 : 1) : (pS2.index() == 3 ? 4 : 3);
//        bc.addCondition(patchSide(pS2.patchIndex(), side), condition_type::dirichlet, 0);

        // Make the Eigenvalue problem to a homogeneous one
        for (gsMultiPatch<>::const_biterator bit = mp_temp.bBegin(); bit != mp_temp.bEnd(); ++bit)
            bc.addCondition(*bit, condition_type::dirichlet, 0);

        gsExprAssembler<real_t> A2(1, 1), B2(1, 1);

        // Elements used for numerical integration
        A2.setIntegrationElements(dbasis_temp);
        B2.setIntegrationElements(dbasis_temp);

        // Set the geometry map
        auto GA = A2.getMap(mp_temp);
        auto GB = B2.getMap(mp_temp);

        // Set the discretization space
        auto uA = A2.getSpace(dbasis_temp);
        auto uB = B2.getSpace(dbasis_temp);

        uA.setup(bc, dirichlet::homogeneous, 0);
        uB.setup(bc, dirichlet::homogeneous,0);
        //uA.setup(0);
        //uB.setup(0);

        A2.initSystem();
        B2.initSystem();

        real_t c = 0.25;
        A2.assembleIfc(mp_temp.interfaces(),
                       c * ilapl(uA.left(), GA.left()) * ilapl(uA.left(), GA.left()).tr() * nv(GA.left()).norm(),
                       c * ilapl(uA.left(), GA.left()) * ilapl(uA.right(), GA.right()).tr() * nv(GA.left()).norm(),
                       c * ilapl(uA.right(), GA.right()) * ilapl(uA.left(), GA.left()).tr() * nv(GA.left()).norm(),
                       c * ilapl(uA.right(), GA.right()) * ilapl(uA.right(), GA.right()).tr() * nv(GA.left()).norm());

        B2.assemble(ilapl(uB, GB) * ilapl(uB, GB).tr() * meas(GB));

        // TODO INSTABLE && SLOW
        Eigen::MatrixXd AA = A2.matrix().toDense().cast<double>();
        Eigen::MatrixXd BB = B2.matrix().toDense().cast<double>();
        Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> ges(AA, BB);

        real_t m_h      = dbasis_temp.basis(0).getMinCellLength(); //*dbasis.basis(0).getMinCellLength();
        mu_interfaces(i,0) = 16.0 * m_h * ges.eigenvalues().array().maxCoeff();
/*
        gsSparseSolver<>::SimplicialLDLT sol;
        sol.compute(B2.matrix());
        gsSparseMatrix<> R = sol.matrixU();
        gsSparseMatrix<> RT = sol.matrixL();
        gsMatrix<> AAA = RT.toDense().inverse() * AA * R.toDense().inverse();

        gsConjugateGradient<> cg(AAA);

        cg.setCalcEigenvalues(true);
        cg.setTolerance(1e-15);
        cg.setMaxIterations(100000);

        gsMatrix<> rhs, result;
        rhs.setRandom( AAA.rows(), 1 );
        result.setRandom( AAA.rows(), 1 );

        cg.solve(rhs,result);

        gsInfo << "Tol: " << cg.error() << "\n";
        gsInfo << "Max it: " << cg.iterations() << "\n";

        gsMatrix<real_t> eigenvalues;
        cg.getEigenvalues(eigenvalues);

        gsInfo << "Cond Number: " << eigenvalues.bottomRows(1)(0,0) << "\n";
*/
    }
}


int main(int argc, char *argv[])
{
    //! [Parse command line]
    bool plot  = false;
    bool mesh  = false;


    index_t numRefine  = 3;
    index_t degree = 3;
    index_t smoothness = 2;

    index_t gluingDataDegree = -1;
    index_t gluingDataSmoothness = -1;

    std::string output;
    std::string fn;

    gsCmdLine cmd("Tutorial on solving a Biharmonic problem with different spaces.");
    // Flags related to the input/geometry
    cmd.addString( "f", "file", "Input geometry file from path (with .xml)", fn );

    // Flags related to the discrete settings
    cmd.addInt( "p", "degree", "Set the polynomial degree of the basis.", degree );
    cmd.addInt( "s", "smoothness", "Set the smoothness of the basis.",  smoothness );
    cmd.addInt( "r", "numRefine", "Number of refinement-loops.",  numRefine );

    cmd.addSwitch("plot", "Create a ParaView visualization file with the solution", plot);
    cmd.addSwitch("mesh", "Plot the mesh", mesh);

    cmd.addString("o", "output", "Output in xml (for python)", output);
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }
    //! [Parse command line]

    //! [Initialize data]
    gsMultiPatch<real_t> mp;

    //! [Read Argument inputs]
    //! [Read geometry]
    std::string string_geo;
    string_geo = fn;

    gsInfo << "Filedata: " << string_geo << "\n";
    gsReadFile<>(string_geo, mp);
    mp.clearTopology();
    mp.computeTopology();
    gsMultiBasis<real_t> dbasis(mp, true);//true: poly-splines (not NURBS)

    // Elevate and p-refine the basis to order p + numElevate
    // where p is the highest degree in the bases
    dbasis.setDegree( degree); // preserve smoothness
    //dbasis.degreeElevate(degree- mp.patch(0).degree(0));

    if (degree-mp.patch(0).degree(0) >0)
        mp.degreeElevate(degree-mp.patch(0).degree(0));
    else if (degree-mp.patch(0).degree(0) <0)
    {
        // Perform
        gsMultiPatch<> result;
        gsFuncProjection<real_t>::L2(dbasis,mp,result);
        gsWriteParaview(result,"mp_corrected",1000,true);
        gsWrite(result,"mp_corrected");
        mp = result;
        gsInfo<<"\nProjected to lower-order basis of degree "<<degree<<"\n\n";
    }

    // h-refine each basis
    for (int r =0; r < numRefine; ++r)
        dbasis.uniformRefine(1, degree-smoothness);

    //! [Export visualization in ParaView]
    if (plot)
    {
        gsInfo<<"Plotting in Paraview...\n";
        gsWriteParaview( mp, "geom",1000,true);
    }
    else
        gsInfo << "Done. No output created, re-run with --plot to get a ParaView "
                  "file containing the solution.\n";
    //! [Export visualization in ParaView]

    // //! [Export data to xml]
    // if (!output.empty())
    // {
    //     index_t cols = method == MethodFlags::NITSCHE ? 7+penalty.cols() : 7;
    //     gsMatrix<real_t> error_collection(l2err.rows(), cols);
    //     error_collection.col(0) = meshsize;
    //     error_collection.col(1) = dofs;
    //     error_collection.col(2) = l2err;
    //     error_collection.col(3) = h1err;
    //     error_collection.col(4) = h2err;
    //     error_collection.col(5) = IFaceErr;
    //     error_collection.col(6) = cond_num;
    //     if (method == MethodFlags::NITSCHE)
    //         error_collection.block(0,7,penalty.rows(),penalty.cols()) = penalty;

    //     gsFileData<real_t> xml_out;
    //     xml_out << error_collection;
    //     xml_out.addString("Meshsize, dofs, l2err, h1err, h2err, iFaceErr, cond_num, (penalty)","Label");
    //     xml_out.addString(std::to_string(degree),"Degree");
    //     xml_out.addString(std::to_string(smoothness),"Regularity");
    //     xml_out.addString(std::to_string(numRefine),"NumRefine");
    //     xml_out.addString(std::to_string(method),"Method");
    //     // Add solution
    //     // [...]
    //     xml_out.save(output);
    //     gsInfo << "XML saved to " + output << "\n";
    // }
    // //! [Export data to xml]

    // if (!write.empty())
    // {
    //     std::ofstream file(write.c_str());
    //     file<<"Meshsize, dofs, l2err, h1err, h2err, iFaceErr"<<"\n";
    //     for (index_t k=0; k<meshsize.size(); ++k)
    //     {
    //         file<<meshsize[k]<<","<<dofs[k]<<","<<l2err[k]<<","<<h1err[k]<<","<<h2err[k]<<","<<IFaceErr[k]<<"\n";
    //     }
    //     file.close();
    // }


    return EXIT_SUCCESS;
}// end main
