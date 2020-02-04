/** @file biharmonic_multiPatch.cpp

    @brief A Biharmonic example for ONLY TWO-PATCHES

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): P. Weinmueller
*/
# include <gismo.h>
# include <omp.h>

#include <gsG1Basis/gsNormL2.h>
#include <gsG1Basis/gsSeminormH1.h>
#include <gsG1Basis/gsSeminormH2.h>
#include <gsG1Basis/gsH1NormWithJump.h>

#include <gsG1Basis/gsG1Basis_mp.h>
#include <gsG1Basis/gsG1BasisLocal_mp.h>
#include <gsAssembler/gsG1BiharmonicAssembler.h>
#include <gsG1Basis/gsG1System_mp.h>

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

    std::vector<std::vector<patchCorner>> allcornerLists;
    for (index_t n = 0; n < multiPatch.nPatches(); ++n)
    {
        for(index_t j=1;j<=4;++j)
        {
            std::vector<patchCorner> cornerLists;
            patchCorner start(n, j);
            multiPatch.getCornerList(start, cornerLists);
            bool alreadyReached = false;
            for(size_t k = 0;k<allcornerLists.size();++k)
                for(size_t l = 0;l<allcornerLists[k].size();++l)
                    if(allcornerLists[k][l].patch==n && allcornerLists[k][l].m_index==j)
                        alreadyReached = true;
            if (cornerLists.size() > 1 && !alreadyReached)
                allcornerLists.push_back(cornerLists);
        }
    }
    for (std::vector<std::vector<patchCorner>>::iterator it = allcornerLists.begin(); it!=allcornerLists.end(); ++it)
    {
        gsInfo << "Corner in " << it->at(0).m_index << " in Patch " << it->at(0).patch << "\n";
        for (std::vector<patchCorner>::iterator it_corner = it->begin(); it_corner!=it->end(); ++it_corner)
        {
            gsInfo << "Corner : " << it_corner->patch << " : " << it_corner->m_index << "\n";
        }
    }


    // FOR LOOP:

std::vector<index_t> level_pTilde;
if (loop)
{
    level_pTilde.push_back(1);
    level_pTilde.push_back(2);
    level_pTilde.push_back(3);
    level_pTilde.push_back(7);
}
else
{
    level_pTilde.push_back(p_tilde);
}

for (index_t level_pt = 0; level_pt < level_pTilde.size(); level_pt++)
{
p_tilde = level_pTilde.at(level_pt);
gsInfo << " ####### p_tilde " << p_tilde << " ######## \n";

std::vector<index_t> level_refine;
if (loop)
{
    level_refine.push_back(1);
    level_refine.push_back(4);
    level_refine.push_back(9);
    level_refine.push_back(19);
    level_refine.push_back(39);
}
else
{
    level_refine.push_back(numRefine);
}
for (index_t level = 0; level < level_refine.size(); level++)
{
    numRefine = level_refine.at(level);
    //gsInfo << " ####### numRefine " << numRefine << " ######## \n";

    gsMultiPatch<> multiPatch;
    fd.getId(0, multiPatch); // id=0: Multipatch domain
    multiPatch.computeTopology();

    // REFINE GEOMETRY
    multiPatch.degreeElevate(numDegree);

    index_t polynomDegree = multiPatch.patch(0).basis().minDegree(); // assume same degree for all patches

    // Add one knot for left patch for non-conforming case
    real_t temp_knot = 0.4;
    gsTensorBSpline<2,real_t> & temp_basis = dynamic_cast<gsTensorBSpline<2,real_t> &>(multiPatch.patch(0));
    temp_basis.insertKnot(temp_knot,1,polynomDegree - regularity);
    // Add one knot for right patch for non-conforming case
    temp_knot = 0.6;
    gsTensorBSpline<2,real_t> & temp_basis2 = dynamic_cast<gsTensorBSpline<2,real_t> &>(multiPatch.patch(1));
    temp_basis2.insertKnot(temp_knot,1,polynomDegree - regularity);

    // refinement
    multiPatch.uniformRefine(numRefine, polynomDegree - regularity);
    multiPatch.patch(0).uniformRefine(1,polynomDegree - regularity); // Left patch one more refine level

    gsMultiBasis<> multiBasis(multiPatch);

    gsInfo << "MultiBasis: " << multiBasis << "\n";
    gsInfo << "MultiBasis: " << multiBasis.basis(1) << "\n";

    if (plot)
    {
        gsMatrix<real_t> coefs;
        gsMultiPatch<real_t> mp_patches_L, mp_patches_R;
        coefs.setZero(multiBasis.basis(0).size(),1);
        for (index_t i = 0; i < coefs.size(); i++)
        {
            coefs.at(i) = 1;
            mp_patches_L.addPatch(multiBasis.basis(0).makeGeometry(coefs));
            coefs.setZero(multiBasis.basis(0).size(),1);
        }
        coefs.setZero(multiBasis.basis(1).size(),1);
        for (index_t i = 0; i < coefs.size(); i++)
        {
            coefs.at(i) = 1;
            mp_patches_R.addPatch(multiBasis.basis(1).makeGeometry(coefs));
            coefs.setZero(multiBasis.basis(1).size(),1);
        }
        std::string baseName = "basisfunction";
        const std::string baseName1(baseName + "_L");
        gsParaviewCollection collection1(baseName1);
        const std::string baseName2(baseName + "_R");
        gsParaviewCollection collection2(baseName2);

        std::string fileName, fileName2;
        for (unsigned i = 0; i < mp_patches_L.nPatches(); i++)
        {

            fileName = baseName1 + "_" + util::to_string(i);

            gsField<> temp_field_L(multiPatch.patch(0),mp_patches_L.patch(i));
            gsWriteParaview(temp_field_L,fileName,5000);
            collection1.addTimestep(fileName,i,"0.vts");

        }
        for (unsigned i = 0; i < mp_patches_R.nPatches(); i++)
        {

            fileName2 = baseName2 + "_" + util::to_string(i);

            gsField<> temp_field_R(multiPatch.patch(1),mp_patches_R.patch(i));
            gsWriteParaview(temp_field_R,fileName2,5000);
            collection2.addTimestep(fileName2,i,"0.vts");

        }
        collection1.save();
        collection2.save();

        //gsField<> field(multiPatch.patch(0),mp_patches);
        //gsField<> field2(multiPatch.patch(1),multiBasis.basis(1));

    }

    gsWriteParaview(multiPatch,"geometry",5000,true);

    std::multimap<index_t, std::map<index_t, std::map<index_t, gsMultiPatch<real_t>>>> basisG1; // first index of patch
                                                                             // second index of side
                                                                             // Third index of interface + basis

    std::vector<index_t> n_tilde, n_bar;


    omp_set_num_threads(threads);
    omp_set_nested(1);

#pragma omp parallel for
    // Compute for each interface the gluing data and the basis functions
    for (unsigned i = 0 ; i < multiPatch.interfaces().size(); i++)  // Iterate over the interfaces
    {
        gsMultiPatch<real_t> bG1_L, bG1_R;
        const boundaryInterface & iFace = multiPatch.interfaces()[i];// assume one single interface

        index_t n_tilde_temp, n_bar_temp;
        if (local_g1)
        {
            gsG1BasisLocal_mp<real_t> g1BasisLocal_mp(multiPatch, multiBasis, iFace, p_tilde, r_tilde, direct, local, plot);
            g1BasisLocal_mp.assemble();
            g1BasisLocal_mp.solve();


            g1BasisLocal_mp.constructSolution(bG1_L,bG1_R);
            n_tilde_temp = g1BasisLocal_mp.get_n_tilde();
            n_bar_temp = g1BasisLocal_mp.get_n_bar();

            if (plot)
                g1BasisLocal_mp.plotG1Basis(bG1_L,bG1_R,"gdBasis_" + std::to_string(i));
        }
        else
        {
            gsG1Basis_mp<real_t> g1Basis_mp(multiPatch, multiBasis, iFace, p_tilde, r_tilde, direct, local, plot);
            g1Basis_mp.assemble();
            g1Basis_mp.solve();


            g1Basis_mp.constructSolution(bG1_L,bG1_R);
            n_tilde_temp = g1Basis_mp.get_n_tilde();
            n_bar_temp = g1Basis_mp.get_n_bar();

            if (plot)
                g1Basis_mp.plotG1Basis(bG1_L,bG1_R,"gdBasis_" + std::to_string(i));
        }

        std::map<index_t, gsMultiPatch<real_t>> l,r;
        l.insert(std::pair< index_t, gsMultiPatch<real_t>>(i, bG1_L));
        r.insert(std::pair< index_t, gsMultiPatch<real_t>>(i, bG1_R));

        std::map<index_t, std::map<index_t, gsMultiPatch<real_t>>> left, right;
        left.insert ( std::pair<index_t, std::map<index_t, gsMultiPatch<real_t>>>(iFace.second().index(), l) );
        right.insert( std::pair<index_t, std::map<index_t, gsMultiPatch<real_t>>>(iFace.first().index(), r));

        basisG1.insert( std::pair<index_t, std::map<index_t, std::map<index_t, gsMultiPatch<real_t>>>>(iFace.second().patch, left));
        basisG1.insert( std::pair<index_t, std::map<index_t, std::map<index_t, gsMultiPatch<real_t>>>>(iFace.first().patch, right));
        n_tilde.push_back(n_tilde_temp);
        n_bar.push_back(n_bar_temp);
    }


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
    g1BiharmonicAssembler.computeDirichletDofsL2Proj(basisG1, n_tilde, n_bar );

    gsG1System_mp<real_t> g1SystemMp(basisG1,
                                     g1BiharmonicAssembler.get_g1dofs(),
                                     g1BiharmonicAssembler.get_mapper(),
                                     g1BiharmonicAssembler.get_mapper_boundary(),
                                     g1BiharmonicAssembler.get_mapper_interface(),
                                     g1BiharmonicAssembler.matrix().dim().first,
                                     n_tilde,
                                     n_bar,
                                     multiBasis);

    g1SystemMp.assemble();

    // Solving system:
    gsMatrix<real_t> f = g1BiharmonicAssembler.rhs(); // with the second boundary condition
    gsSparseMatrix<real_t> K_sparse = g1BiharmonicAssembler.matrix();


    gsSparseMatrix<real_t> D_0_sparse = g1SystemMp.get_D_0_sparse();
    gsSparseMatrix<real_t> D_boundary_sparse = g1SystemMp.get_D_boundary_sparse();

    gsMatrix<real_t> g = g1SystemMp.get_g();

    //gsInfo << "Solving system... \n";
    gsSparseMatrix<real_t> A = D_0_sparse * K_sparse * D_0_sparse.transpose();
    gsVector<real_t> F = D_0_sparse * f - D_0_sparse * K_sparse * D_boundary_sparse.transpose() * g;

    //gsInfo << "System finished with " << A.dim() << " non-zeros!\n";

    gsSparseSolver<real_t>::CGDiagonal solver;
    //gsSparseSolver<real_t>::LU solver;
    //solver.analyzePattern(BiharmonicAssembler.matrix() );
    //solver.factorize(BiharmonicAssembler.matrix());
    //gsInfo << "matrix: " << K_sparse.dim() << "\n";
    //gsInfo << "rhs: " << F << "\n";
    solver.compute(A);
    gsMatrix<> solVector= solver.solve(F);
    //gsInfo << "rhs: " << F << "\n";
    //gsInfo << "Solving finished! \n";

    gsMultiPatch<> mpsol;
    g1BiharmonicAssembler.constructSolution(solVector.bottomRows(g1BiharmonicAssembler.matrix().dim().first),mpsol);
    gsField<> solField(multiPatch, mpsol);

    g1SystemMp.constructSolution_G1(solVector.topRows(g1SystemMp.get_dim_g1()),basisG1);

    if (plot)
    {
        //gsInfo<<"Plotting in Paraview...\n";
        //gsWriteParaview(solField,"biharmonic_trafo",5000);

        const gsField<> exact( multiPatch, solution, false );
        gsWriteParaview<>( exact, "Biharmonic2d_exact", 5000);

        g1BiharmonicAssembler.writeParaview(solField,"biharmonic_trafo_mp_g1",basisG1,5000);

    }
    // Error analysis

    //gsInfo << "L2 error: H1 Seminorm error: H2 Seminorm error: H2 Norm error: errorJump error: \n";

    omp_set_nested(1);
    omp_set_num_threads(1);
    //gsInfo << "Starting the error computing! \n";
#pragma omp parallel for
    for (index_t e = 0; e < 4; ++e)
    {
        if (e == 0)
        {

            // ERROR als visitor schreiben usw... dann kann man das
            // Parallelisieren
            gsNormL2<real_t> error(solField, solVal, basisG1);
            error.compute();
            gsInfo << error.value();
        }
        else if (e == 1)
        {
            gsSeminormH1<real_t> errorH1(solField, solVal, basisG1);
            errorH1.compute();
            gsInfo << " " << errorH1.value();
        }
        else if (e == 2)
        {
            gsSeminormH2<real_t> errorH2(solField, solVal, basisG1);
            errorH2.compute();
            gsInfo << " " << errorH2.value();
        }
        else if (e == 3)
        {

        }

    }
    for (index_t i = 0; i < multiPatch.interfaces().size(); ++i)
    {
        const boundaryInterface & iFace = multiPatch.interfaces()[i];
        gsH1NormWithJump<real_t> errorJump(solField, basisG1, iFace);
        errorJump.compute();
        gsInfo << " " << errorJump.value() << "\n";
    }
    //gsInfo << " " << math::sqrt(errorH2.value()*errorH2.value() +
    //errorH1.value()*errorH1.value() + error.value()*error.value());

} // level
} // level_pt

} // main