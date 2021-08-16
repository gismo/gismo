/** @file gsCompositeBasis_test.h

    @brief File testing the gsCompositeBasis class.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): F. Buchegger
*/

#include <gismo.h>

#include <gsMSplines/gsCompositeIncrSmoothnessBasis.h>
#include <gsMSplines/gsCompositeIncrSmoothnessGeom.h>

#include <gsMSplines/gsDPatch.h>

#include <gsKLShell/gsThinShellAssembler.h>
#include <gsKLShell/gsMaterialMatrixLinear.h>


using namespace gismo;


void writeToCSVfile(std::string name, gsMatrix<> matrix)
{
    std::ofstream file(name.c_str());
    for(int  i = 0; i < matrix.rows(); i++){
        for(int j = 0; j < matrix.cols(); j++){
           std::string str = std::to_string(matrix(i,j));
           if(j+1 == matrix.cols()){
               file<<str;
           }else{
               file<<str<<',';
           }
        }
        file<<'\n';
    }
  }


int main(int argc, char *argv[])
{
    bool plot       = false;
    bool last       = false;
    index_t numRefine  = 2;
    index_t numElevate = 1;
    index_t geometry = 1;
    index_t smoothing = 0;
    std::string input;

    std::string fn("planar/multipatch_triangle.xml");

    gsCmdLine cmd("Composite basis tests.");
    cmd.addPlainString("filename", "G+Smo input geometry file.", input);
    cmd.addInt( "e", "degreeElevation",
                "Number of degree elevation steps to perform before solving (0: equalize degree in all directions)", numElevate );
    cmd.addInt( "r", "uniformRefine", "Number of Uniform h-refinement steps to perform before solving",  numRefine );
    cmd.addInt( "g", "geometry", "which geometry",  geometry );
    cmd.addInt( "s", "smooth", "Smoothing method to use",  smoothing );
    cmd.addSwitch("plot", "plot",plot);
    cmd.addSwitch("last", "last case only",last);
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    gsMultiPatch<> mp;

    gsBoundaryConditions<> bc;

    if (input.empty())
    {
        if (geometry==0)
        {
            mp.addPatch( gsNurbsCreator<>::BSplineSquare(1, 1, 1) ); // patch 0
            // mp.embed(3);
            mp.computeTopology();

            bc.addCondition(0,boundary::north, condition_type::dirichlet, 0 );
            bc.addCondition(0,boundary::south, condition_type::dirichlet, 0 );
            bc.addCondition(0,boundary::west, condition_type::dirichlet, 0 );
            bc.addCondition(0,boundary::east, condition_type::dirichlet, 0 );
        }
        else if (geometry==1)
        {
            mp.addPatch( gsNurbsCreator<>::BSplineSquare(1) ); // patch 0
            mp.addPatch( gsNurbsCreator<>::BSplineSquare(1) ); // patch 1
            // mp.embed(3);

            gsMatrix<> ones;
            ones = gsMatrix<>::Ones(mp.patch(1).coefs().rows(),1); // patch 1

            // translate second patch to the right
            mp.patch(1).coefs().col(0) += ones; // patch 1

            mp.addInterface(0,boundary::side::east,1,boundary::side::west);

            mp.addAutoBoundaries();

            bc.addCondition(0,boundary::north, condition_type::dirichlet, 0 );
            bc.addCondition(0,boundary::south, condition_type::dirichlet, 0 );
            bc.addCondition(0,boundary::west, condition_type::dirichlet, 0 );

            bc.addCondition(1,boundary::north, condition_type::dirichlet, 0 );
            bc.addCondition(1,boundary::south, condition_type::dirichlet, 0 );
            bc.addCondition(1,boundary::east, condition_type::dirichlet, 0 );
        }
        else if (geometry==2)
        {
            mp.addPatch( gsNurbsCreator<>::BSplineSquare(1) ); // patch 0
            mp.addPatch( gsNurbsCreator<>::BSplineSquare(1) ); // patch 1
            mp.addPatch( gsNurbsCreator<>::BSplineSquare(1) ); // patch 2
            // mp.embed(3);

            gsMatrix<> ones;
            ones = gsMatrix<>::Ones(mp.patch(1).coefs().rows(),1); // patch 1

            // translate second patch up
            mp.patch(1).coefs().col(1) += ones; // patch 1

            // translate third patch up and left
            ones = gsMatrix<>::Ones(mp.patch(2).coefs().rows(),1); // patch 2
            mp.patch(2).coefs().col(0) -= ones; // patch 2
            mp.patch(2).coefs().col(1) += ones; // patch 2


            mp.addInterface(0,boundary::side::north,1,boundary::side::south);
            mp.addInterface(1,boundary::side::west,2,boundary::side::east);

            mp.addAutoBoundaries();

            // // Needed
            // mp.uniformRefine(1);

            bc.addCondition(0,boundary::south, condition_type::dirichlet, 0 );
            bc.addCondition(0,boundary::west, condition_type::dirichlet, 0 );
            bc.addCondition(0,boundary::east, condition_type::dirichlet, 0 );

            bc.addCondition(1,boundary::north, condition_type::dirichlet, 0 );
            bc.addCondition(1,boundary::east, condition_type::dirichlet, 0 );

            bc.addCondition(2,boundary::north, condition_type::dirichlet, 0 );
            bc.addCondition(2,boundary::west, condition_type::dirichlet, 0 );
            bc.addCondition(2,boundary::south, condition_type::dirichlet, 0 );

        }
        else if (geometry==3)
        {
            mp.addPatch( gsNurbsCreator<>::BSplineSquare(1) ); // patch 0
            mp.addPatch( gsNurbsCreator<>::BSplineSquare(1) ); // patch 1
            mp.addPatch( gsNurbsCreator<>::BSplineSquare(1) ); // patch 2
            mp.addPatch( gsNurbsCreator<>::BSplineSquare(1) ); // patch 3
            // mp.embed(3);

            gsMatrix<> ones;
            ones = gsMatrix<>::Ones(mp.patch(1).coefs().rows(),1); // patch 1

            // translate second patch to the right
            mp.patch(1).coefs().col(0) += ones; // patch 1

            // translate fourth patch to the right
            ones = gsMatrix<>::Ones(mp.patch(2).coefs().rows(),1); // patch 2
            mp.patch(2).coefs().col(1) += ones; // patch 2

            // translate fourth patch to the right
            ones = gsMatrix<>::Ones(mp.patch(3).coefs().rows(),1); // patch 3
            mp.patch(3).coefs().col(0) += ones; // patch 3
            mp.patch(3).coefs().col(1) += ones; // patch 3

            mp.addInterface(0,boundary::side::north,2,boundary::side::south);
            mp.addInterface(0,boundary::side::east,1,boundary::side::west);
            mp.addInterface(2,boundary::side::east,3,boundary::side::west);
            mp.addInterface(1,boundary::side::north,3,boundary::side::south);

            mp.addAutoBoundaries();

            bc.addCondition(0,boundary::west, condition_type::dirichlet, 0 );
            bc.addCondition(0,boundary::south, condition_type::dirichlet, 0 );

            bc.addCondition(1,boundary::south, condition_type::dirichlet, 0 );
            bc.addCondition(1,boundary::east, condition_type::dirichlet, 0 );

            bc.addCondition(2,boundary::north, condition_type::dirichlet, 0 );
            bc.addCondition(2,boundary::west, condition_type::dirichlet, 0 );

            bc.addCondition(3,boundary::north, condition_type::dirichlet, 0 );
            bc.addCondition(3,boundary::east, condition_type::dirichlet, 0 );
        }
        else if (geometry==4)
        {
            gsReadFile<>("planar/hexagon_3p.xml", mp);
            // mp.embed(3);
            mp.clearTopology();
            mp.computeTopology();
        mp.degreeElevate(numElevate);

            bc.addCondition(0,boundary::north, condition_type::dirichlet, 0 );
            bc.addCondition(0,boundary::east, condition_type::dirichlet, 0 );

            bc.addCondition(1,boundary::north, condition_type::dirichlet, 0 );
            bc.addCondition(1,boundary::west, condition_type::dirichlet, 0 );

            bc.addCondition(2,boundary::south, condition_type::dirichlet, 0 );
            bc.addCondition(2,boundary::west, condition_type::dirichlet, 0 );
        }
        else if (geometry==5)
        {
            gsReadFile<>("planar/hexagon_5p.xml", mp);
            // mp.embed(3);
            mp.clearTopology();
            mp.computeTopology();
        mp.degreeElevate(numElevate);

            bc.addCondition(0,boundary::south, condition_type::dirichlet, 0 );
            bc.addCondition(0,boundary::east, condition_type::dirichlet, 0 );

            bc.addCondition(1,boundary::south, condition_type::dirichlet, 0 );
            bc.addCondition(1,boundary::west, condition_type::dirichlet, 0 );

            bc.addCondition(2,boundary::north, condition_type::dirichlet, 0 );
            bc.addCondition(2,boundary::east, condition_type::dirichlet, 0 );

            bc.addCondition(3,boundary::north, condition_type::dirichlet, 0 );
            bc.addCondition(3,boundary::west, condition_type::dirichlet, 0 );

        }
        else if (geometry==6)
        {
            gsReadFile<>("planar/triangle2d_u3p.xml", mp);
            // mp.embed(3);
            mp.clearTopology();
            mp.computeTopology();

            bc.addCondition(0,boundary::south, condition_type::dirichlet, 0 );
            bc.addCondition(0,boundary::west, condition_type::dirichlet, 0 );

            bc.addCondition(1,boundary::south, condition_type::dirichlet, 0 );
            bc.addCondition(1,boundary::east, condition_type::dirichlet, 0 );

            bc.addCondition(2,boundary::south, condition_type::dirichlet, 0 );
            bc.addCondition(2,boundary::east, condition_type::dirichlet, 0 );
        }
        else
            GISMO_ERROR("Geometry with index "<<geometry<<" unknown.");
    }
    else
    {
        gsReadFile<>(input, mp);
        mp.clearTopology();
        mp.computeTopology();
    }

    mp.embed(3);

    gsMultiPatch<> geom = mp;

    // p-refine
    if (numElevate!=0)
        mp.degreeElevate(numElevate);

    if (last)
    {
        // h-refine
        for (int r =0; r < numRefine; ++r)
            mp.uniformRefine();

        numRefine = 0;
    }

    gsWriteParaview(mp,"mp",1000,true,false);

    // NOTE: Points go wrong!
    gsPointLoads<real_t> pLoads = gsPointLoads<real_t>();
    gsVector<> point(2); point<< 1, 0.5 ;
    gsVector<> load (3); load << 0.0, 0.0, 0.0 ;
    pLoads.addLoad(point, load, 0 );

    real_t thickness = 1.0;
    real_t E_modulus = 1.0;
    real_t PoissonRatio = 0.0;
    gsFunctionExpr<> force("0","0","1",3);
    gsFunctionExpr<> t(std::to_string(thickness),3);
    gsFunctionExpr<> E(std::to_string(E_modulus),3);
    gsFunctionExpr<> nu(std::to_string(PoissonRatio),3);


    std::vector<gsFunction<>*> parameters(2);
    parameters[0] = &E;
    parameters[1] = &nu;

    gsMaterialMatrixLinear<3,real_t> materialMatrix(mp,t,parameters);

    gsThinShellAssemblerBase<real_t>* assembler;

    //! [Solver loop]
    gsVector<> l2err(numRefine+1), h1err(numRefine+1), linferr(numRefine+1),
        b2err(numRefine+1), b1err(numRefine+1), binferr(numRefine+1);

    gsSparseSolver<>::CGDiagonal solver;

    gsVector<> solVector;

    gsMappedBasis<2,real_t> bb2;

    gsSparseMatrix<> global2local;
    gsMatrix<> coefs;
    gsStopwatch time;

    gsMultiBasis<> dbasis(mp);

    for( index_t r = 0; r<=numRefine; ++r)
    {

        gsInfo<<"--------------------------------------------------------------\n";
        time.restart();
        if (smoothing==0)
        {
            gsCompositeIncrSmoothnessGeom<2,real_t> cgeom(mp,3);
            gsCompositeIncrSmoothnessBasis<2,real_t> & basis = dynamic_cast<gsCompositeIncrSmoothnessBasis<2,real_t>&>(cgeom.getCompBasis());

            global2local = basis.getMapper().asMatrix();
            geom = cgeom.exportToPatches();
        }
        else if (smoothing==1)
        {
            gsDPatch<2,real_t> dpatch(mp);
            dpatch.matrix_into(global2local);
            global2local = global2local.transpose();
            geom = dpatch.exportToPatches();
            dbasis = dpatch.localBasis();
        }
        else
            GISMO_ERROR("Option "<<smoothing<<" for smoothing does not exist");

        gsInfo<<"\tAssembly of mapping:\t"<<time.stop()<<"\t[s]\n";

        bb2.init(dbasis,global2local);
        // gsMappedSpline<2,real_t> mspline(bb2,coefs);
        // geom = mspline.exportToPatches();


        assembler = new gsThinShellAssembler<3, real_t, true>(mp,dbasis,bc,force,&materialMatrix);
        assembler->setSpaceBasis(bb2);
        assembler->setPointLoads(pLoads);

        time.restart();
        // Initialize the system

        assembler->assemble();
        gsSparseMatrix<> matrix = assembler->matrix();
        gsVector<> vector = assembler->rhs();

        gsInfo<<"\tSystem assembly:\t"<<time.stop()<<"\t[s]\n";

        time.restart();

        solver.compute( matrix );
        solVector = solver.solve(vector);

        gsInfo<<"\tSolving system:\t\t"<<time.stop()<<"\t[s]\n";
        time.restart();

        // l2err[r]= math::sqrt( ev.integral( (f - s).sqNorm()*meas(G) ) / ev.integral(f.sqNorm()*meas(G)) );

        // h1err[r]= l2err[r] + math::sqrt(ev.integral( ( igrad(f) - grad(s)*jac(G).inv() ).sqNorm()*meas(G) )/ev.integral( igrad(f).sqNorm()*meas(G) ) );

        // linferr[r] = ev.max( f-s ) / ev.max(f);

        gsInfo<<"\tError computations:\t"<<time.stop()<<"\t[s]\n"; // This takes longer for the D-patch, probably because there are a lot of points being evaluated, all containing the linear combinations of the MSplines

        // TO DO: Refine the mb
        mp.uniformRefine();
        dbasis = gsMultiBasis<>(mp);

        gsWriteParaview(dbasis.basis(0),"basis");
    }
    //! [Solver loop]

    // gsInfo<< "\n\nCG it.: "<<CGiter.transpose()<<"\n";

    // //! [Error and convergence rates]
    // gsInfo<<"\n* Error\n";
    // gsInfo<< "H1    "<<std::scientific<<h1err.transpose()<<"\n";
    // gsInfo<< "L2    "<<std::scientific<<std::setprecision(3)<<l2err.transpose()<<"\n";
    // gsInfo<< "Linf  "<<std::scientific<<linferr.transpose()<<"\n";

    // gsInfo<<"\n* EoC\n";
    // gsInfo<< "H1c   0 "<< std::fixed<<std::setprecision(2)
    //       <<( h1err.head(numRefine).array() /
    //           h1err.tail(numRefine).array() ).log().transpose() / std::log(2.0) <<"\n";
    // gsInfo<< "L2c   0 " << std::fixed<<std::setprecision(2)
    //       << ( l2err.head(numRefine).array() /
    //            l2err.tail(numRefine).array() ).log().transpose() / std::log(2.0) <<"\n";

    // gsInfo<<   "Linfc 0 "<< std::fixed<<std::setprecision(2)
    //       <<( linferr.head(numRefine).array() /
    //           linferr.tail(numRefine).array() ).log().transpose() / std::log(2.0) <<"\n";

    //! [Error and convergence rates]

    //! [Export visualization in ParaView]
    if (plot)
    {
        assembler->plotSolution("solution", solVector);

        gsMultiPatch<> deformation = assembler->constructDisplacement(solVector);

        gsField<> solField(geom, deformation);
        gsInfo<<"Plotting in Paraview...\n";
        gsWriteParaview<>( solField, "Deformation", 1000, true);

        deformation = assembler->constructSolution(solVector);
        gsWriteParaview<>( deformation, "deformed_geom", 1000, true);


        gsWriteParaview( geom, "geom",100,true);
    }
    //! [Export visualization in ParaView]

    return EXIT_SUCCESS;


}
