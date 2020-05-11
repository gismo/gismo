/** @file approxProjection_example.cpp

    @brief Test the rate for the approximated L2 porjection

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): P. Weinmueller
*/
# include <omp.h>

# include <gismo.h>

# include <gsG1Basis/gsApproxProjectionAssembler.h>
# include <gsG1Basis/gsApproxNormL2.h>

using namespace gismo;

int main(int argc, char *argv[])
{
    // Geometry data
    index_t loop = 1; // Number of refinement steps
    index_t geometry = 1; // Which geometry

    index_t numRefine = 4; // Initial refinement
    index_t regularity = 1; // Regularity

    // For the spline space of the gluing data
    index_t p_tilde = 1;
    index_t r_tilde = 0;

    // Which basis function should computed
    index_t basisID = 0;
    index_t degree_target = 0;

    index_t threads = 4; // For parallel computing

    index_t gluingData = 0; // global

    bool plot = false;
    bool localApprox = false;

    gsCmdLine cmd("Example for solving the biharmonic problem.");
    cmd.addInt("k", "refine", "Number of refinement steps", numRefine);
    cmd.addInt("p", "p_tilde", "Polynomial degree for tilde{p}", p_tilde);
    cmd.addInt("r", "r_tilde", "Regularity for tilde{r}", r_tilde);
    cmd.addInt("g", "geometry", "Geometry", geometry);
    cmd.addInt( "l", "loop", "The number of refinement steps", loop);
    cmd.addInt("b", "basisID", "Which basis functions should computed", basisID);
    cmd.addInt( "d", "degree_target", "degree_target", degree_target);
    cmd.addInt( "a", "gluingData", "gluingData", gluingData);
    cmd.addInt( "t", "threads", "threads", threads);
    cmd.addSwitch( "plot", "Plot result in ParaView format", plot );
    cmd.addSwitch( "localApprox", "localApprox", localApprox );
    try { cmd.getValues(argc,argv); } catch (int rv) {  }

    gsOptionList optionList;
    optionList.addInt("loop","Loop", loop);
    optionList.addInt("gluingData","gluingData",gluingData);


    optionList.addInt("numRefine","Number of refinement", numRefine);

    optionList.addInt("p_tilde","Grad",p_tilde);
    optionList.addInt("r_tilde","Reg",r_tilde);
    optionList.addInt("regularity","Regularity of the initial geometry",regularity);
    optionList.addInt("basisID","basisID", basisID);

    optionList.addSwitch("plot","Plot in Paraview",plot);
    optionList.addSwitch("localApprox", "Local approx projection", localApprox);

    // ======= Geometry =========
    std::string string_geo;
    index_t numDegree = 0;
    switch(geometry)
    {
        case 0:
            string_geo = "planar/twoPatches/square_diagonal.xml";
            numDegree = 2;
            break;
        case 1:
            string_geo = "planar/twoPatches/square_curved.xml";
            numDegree = 0;
            break;
        case 2:
            string_geo = "planar/twoPatches/2patch_curved.xml";
            numDegree = 0;
            break;
        case 3:
            string_geo = "planar/twoPatches/2patch_C1curved.xml";
            numDegree = 0;
            break;


        default:
            gsInfo << "No geometry is used! \n";
            break;
    }
    gsFileData<> fd(string_geo);
    gsInfo << "Loaded file "<< fd.lastPath() <<"\n";

    gsMultiPatch<> multiPatch_init;
    fd.getId(0, multiPatch_init); // id=0: Multipatch domain
    multiPatch_init.computeTopology();

    //multiPatch.patch(1).degreeElevate(1,0);
    multiPatch_init.degreeElevate(numDegree + degree_target);
    //optionList.setInt("regularity",2);

    gsVector<index_t> num_knots(optionList.getInt("loop"));
    num_knots[0] = optionList.getInt("numRefine");
    for (index_t i = 1; i < optionList.getInt("loop"); i++)
        num_knots[i] = num_knots[i-1]*2 + 1;

    gsVector<real_t> l2Error_vec(optionList.getInt("loop") + 1);
    gsVector<real_t> h1SemiError_vec(optionList.getInt("loop") + 1);
    l2Error_vec.setZero();
    h1SemiError_vec.setZero();

    for (index_t refinement_level = 0; refinement_level < optionList.getInt("loop"); refinement_level++)
    {
        gsMultiPatch<> multiPatch(multiPatch_init);
        multiPatch.uniformRefine_withSameRegularity(num_knots[refinement_level], optionList.getInt("regularity"));

        gsInfo << "###### Level: " << refinement_level << " with " << num_knots[refinement_level]
               << " inner knots ###### " << "\n";

        gsMultiBasis<> mb(multiPatch);

        // PREPARE THE SPACE FOR PROJECTION
        index_t m_uv = 1; // v
        gsBSplineBasis<> basis_target = dynamic_cast<gsBSplineBasis<> &>(mb.basis(0).component(m_uv)); // 0 -> v, 1 -> u

        // Compute the assembler
        gsApproxProjectionAssembler<real_t> approxProjectionAssembler(basis_target, multiPatch, optionList);

        gsInfo << "Solving...\n";
        gsSparseSolver<real_t>::CGDiagonal solver;
        gsMatrix<> sol;
        solver.compute(approxProjectionAssembler.matrix());
        sol = solver.solve(approxProjectionAssembler.rhs());
        gsInfo << "Finish solving...\n";

        gsMultiPatch<> result_mp;
        approxProjectionAssembler.constructSolution(sol, result_mp);

        gsMatrix<> points(1, 1);
        points << 0;
        gsInfo << "Eval at zero: " << result_mp.patch(0).eval(points) << "\n";

        if (plot)
        {
            gsInfo << "Plot paraview\n";
            gsWriteParaview(result_mp, "result_approxProjection", 5000);

        }

#ifdef _OPENMP
        omp_set_num_threads(threads);
        omp_set_nested(1);
#endif
        gsApproxNormL2<real_t> approxNormL2(multiPatch, result_mp);
        approxNormL2.compute(optionList);
        l2Error_vec[refinement_level] = approxNormL2.value();
        //h1SemiError_vec[refinement_level] = errorSemiH1.value();
    }


    if (optionList.getInt("loop") > 1)
    {
        gsInfo << "=====================================================================\n";
        gsMatrix<> rate(optionList.getInt("loop") + 1,1);
        rate.setZero();

        printf("%-5d & %-14.6e & %-5.2f \\\\ \n", num_knots[0],
               l2Error_vec[0], rate(0,0));
        for (index_t i = 1; i < optionList.getInt("loop"); i++)
        {
            rate(i,0) = log2(l2Error_vec[i-1] / l2Error_vec[i]);
            printf("%-5d & %-14.6e & %-5.2f  \\\\ \n", num_knots[i],
                   l2Error_vec[i], rate(i,0));
        }

    }
    else
    {
        gsInfo << "=====================================================================\n";
        gsInfo << "L2 Error: " << l2Error_vec[0] << "\n";
        gsInfo << "H1 Semi-error: " << h1SemiError_vec[0] << "\n";
        gsInfo << "=====================================================================\n";

    }

}

