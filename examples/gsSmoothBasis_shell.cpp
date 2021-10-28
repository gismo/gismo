/** @file gsCompositeBasis_test.h

    @brief File testing the gsCompositeBasis class.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): F. Buchegger
*/

#include <gismo.h>

// #include <gsMSplines/gsCompositeIncrSmoothnessBasis.h>
// #include <gsMSplines/gsCompositeIncrSmoothnessGeom.h>

#include <gsMSplines/gsDPatch.h>

#include <gsUnstructuredSplines/gsApproxC1Spline.h>
#include <gsUnstructuredSplines/gsC1SurfSpline.h>

#include <gsKLShell/gsThinShellAssembler.h>
#include <gsKLShell/gsMaterialMatrixLinear.h>
#include <gsKLShell/gsFunctionSum.h>

#include <gsUtils/gsQuasiInterpolate.h>


#include <gsAssembler/gsExprAssembler.h>



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
    bool info       = false;
    bool writeMatrix= false;
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
    cmd.addSwitch("writeMat", "Write projection matrix",writeMatrix);
    cmd.addSwitch( "info", "Print information", info );

    // to do:
    // smoothing method add nitsche @Pascal


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

            for (index_t d = 0; d!=3; d++)
            {
                bc.addCondition(0,boundary::north, condition_type::dirichlet, 0, 0, false, d );
                bc.addCondition(0,boundary::south, condition_type::dirichlet, 0, 0, false, d);
                bc.addCondition(0,boundary::west, condition_type::dirichlet, 0, 0, false, d);
                bc.addCondition(0,boundary::east, condition_type::dirichlet, 0, 0, false, d);
            }
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

            // gsReadFile<>("domain2d/twoSquares.xml", mp);

            mp.addAutoBoundaries();
            for (index_t d = 0; d!=3; d++)
            {
                bc.addCondition(0,boundary::north, condition_type::dirichlet, 0, 0, false, d);
                bc.addCondition(0,boundary::south, condition_type::dirichlet, 0, 0, false, d);
                bc.addCondition(0,boundary::west, condition_type::dirichlet, 0, 0, false, d);

                bc.addCondition(1,boundary::north, condition_type::dirichlet, 0, 0, false, d);
                bc.addCondition(1,boundary::south, condition_type::dirichlet, 0, 0, false, d);
                bc.addCondition(1,boundary::east, condition_type::dirichlet, 0, 0, false, d);
            }
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
            for (index_t d = 0; d!=3; d++)
            {
                bc.addCondition(0,boundary::south, condition_type::dirichlet, 0, 0, false, d);
                bc.addCondition(0,boundary::west, condition_type::dirichlet, 0, 0, false, d);
                bc.addCondition(0,boundary::east, condition_type::dirichlet, 0, 0, false, d);

                bc.addCondition(1,boundary::north, condition_type::dirichlet, 0, 0, false, d);
                bc.addCondition(1,boundary::east, condition_type::dirichlet, 0, 0, false, d);

                bc.addCondition(2,boundary::north, condition_type::dirichlet, 0, 0, false, d);
                bc.addCondition(2,boundary::west, condition_type::dirichlet, 0, 0, false, d);
                bc.addCondition(2,boundary::south, condition_type::dirichlet, 0, 0, false, d);
            }

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

            for (index_t d = 0; d!=3; d++)
            {
                bc.addCondition(0,boundary::west, condition_type::dirichlet, 0, 0, false, d);
                bc.addCondition(0,boundary::south, condition_type::dirichlet, 0, 0, false, d);

                bc.addCondition(1,boundary::south, condition_type::dirichlet, 0, 0, false, d);
                bc.addCondition(1,boundary::east, condition_type::dirichlet, 0, 0, false, d);

                bc.addCondition(2,boundary::north, condition_type::dirichlet, 0, 0, false, d);
                bc.addCondition(2,boundary::west, condition_type::dirichlet, 0, 0, false, d);

                bc.addCondition(3,boundary::north, condition_type::dirichlet, 0, 0, false, d);
                bc.addCondition(3,boundary::east, condition_type::dirichlet, 0, 0, false, d);
            }
        }
        else if (geometry==4)
        {
            gsReadFile<>("planar/hexagon_3p.xml", mp);
            // mp.embed(3);
            mp.clearTopology();
            mp.computeTopology();
            mp.degreeElevate(numElevate);

            for (index_t d = 0; d!=3; d++)
            {
                bc.addCondition(0,boundary::north, condition_type::dirichlet, 0, 0, false, d);
                bc.addCondition(0,boundary::east, condition_type::dirichlet, 0, 0, false, d);

                bc.addCondition(1,boundary::north, condition_type::dirichlet, 0, 0, false, d);
                bc.addCondition(1,boundary::west, condition_type::dirichlet, 0, 0, false, d);

                bc.addCondition(2,boundary::south, condition_type::dirichlet, 0, 0, false, d);
                bc.addCondition(2,boundary::west, condition_type::dirichlet, 0, 0, false, d);
            }
        }
        else if (geometry==5)
        {
            gsReadFile<>("planar/hexagon_5p.xml", mp);
            // mp.embed(3);
            mp.clearTopology();
            mp.computeTopology();
            mp.degreeElevate(numElevate);

            for (index_t d = 0; d!=3; d++)
            {
                bc.addCondition(0,boundary::south, condition_type::dirichlet, 0, 0, false, d);
                bc.addCondition(0,boundary::east, condition_type::dirichlet, 0, 0, false, d);

                bc.addCondition(1,boundary::south, condition_type::dirichlet, 0, 0, false, d);
                bc.addCondition(1,boundary::west, condition_type::dirichlet, 0, 0, false, d);

                bc.addCondition(2,boundary::north, condition_type::dirichlet, 0, 0, false, d);
                bc.addCondition(2,boundary::east, condition_type::dirichlet, 0, 0, false, d);

                bc.addCondition(3,boundary::north, condition_type::dirichlet, 0, 0, false, d);
                bc.addCondition(3,boundary::west, condition_type::dirichlet, 0, 0, false, d);
            }

        }
        else if (geometry==6)
        {
            gsReadFile<>("planar/triangle2d_u3p.xml", mp);
            // mp.embed(3);
            mp.clearTopology();
            mp.computeTopology();

            for (index_t d = 0; d!=3; d++)
            {
                bc.addCondition(0,boundary::south, condition_type::dirichlet, 0, 0, false, d);
                bc.addCondition(0,boundary::west, condition_type::dirichlet, 0, 0, false, d);

                bc.addCondition(1,boundary::south, condition_type::dirichlet, 0, 0, false, d);
                bc.addCondition(1,boundary::east, condition_type::dirichlet, 0, 0, false, d);

                bc.addCondition(2,boundary::south, condition_type::dirichlet, 0, 0, false, d);
                bc.addCondition(2,boundary::east, condition_type::dirichlet, 0, 0, false, d);
            }
        }
        else if (geometry==7)
        {
            real_t L = 1;
            index_t N = 5;
            index_t M = 5;
            mp = gsNurbsCreator<>::BSplineSquareGrid(N,M,L,0,0);
            mp.computeTopology();

            gsMatrix<> bbox, pbbox;
            mp.boundingBox(bbox);
            // std::vector<index_t> boundaryPatches(mp.nPatches());
            for (size_t p = 0; p!= mp.nPatches(); p++)
            {
                gsMultiPatch<> mp_tmp(mp.patch(p));
                mp_tmp.boundingBox(pbbox);
                if ( (bbox(0,0) - pbbox(0,0)) ==0 )
                {
                    gsInfo<<"Patch "<<p<<" is a boundary patch on west side!\n";
                    for (index_t d = 0; d!=3; d++)
                    {
                        bc.addCondition(p,boundary::west, condition_type::dirichlet, 0, 0, false, d);
                    }
                }
                if ( (bbox(0,1) - pbbox(0,1)) ==0 )
                {
                    gsInfo<<"Patch "<<p<<" is a boundary patch on east side!\n";
                    for (index_t d = 0; d!=3; d++)
                    {
                        bc.addCondition(p,boundary::east, condition_type::dirichlet, 0, 0, false, d);
                    }
                }
                if ( (bbox(1,0) - pbbox(1,0)) ==0 )
                {
                    gsInfo<<"Patch "<<p<<" is a boundary patch on south side!\n";
                    for (index_t d = 0; d!=3; d++)
                    {
                        bc.addCondition(p,boundary::south, condition_type::dirichlet, 0, 0, false, d);
                    }
                }
                if ( (bbox(1,1) - pbbox(1,1)) ==0 )
                {
                    gsInfo<<"Patch "<<p<<" is a boundary patch on north side!\n";
                    for (index_t d = 0; d!=3; d++)
                    {
                        bc.addCondition(p,boundary::north, condition_type::dirichlet, 0, 0, false, d);
                    }
                }
                // if ( (bbox(0,0) - pbbox(0,0)) * (bbox(0,1) - pbbox(0,1)) * (bbox(0,1) - pbbox(0,1)) * (bbox(1,1) - pbbox(1,1)) != 0 )
                //     boundaryPatches[p] = 0;
            }
        }
        else if (geometry==8)
        {
            index_t N = 5;
            index_t M = 5;
            gsReadFile<>("planar/diamondTile.xml", mp);

            std::vector<gsMultiPatch<>> container(2);
            gsMultiPatch<> mp_copy = mp;
            mp.clear();
            mp = mp_copy;
            gsNurbsCreator<>::mirror2D(mp,1);
            container[0] = mp_copy;
            container[1] = mp;
            mp = gsNurbsCreator<>::makeGrid(container,N,M);

            gsMatrix<> bbox, pbbox;
            mp.boundingBox(bbox);
            // std::vector<index_t> boundaryPatches(mp.nPatches());
            for (size_t p = 0; p!= mp.nPatches(); p++)
            {
                gsMultiPatch<> mp_tmp(mp.patch(p));
                mp_tmp.boundingBox(pbbox);
                if ( (bbox(0,0) - pbbox(0,0)) ==0 )
                {
                    gsInfo<<"Patch "<<p<<" is a boundary patch on west side!\n";
                    for (index_t d = 0; d!=3; d++)
                    {
                        bc.addCondition(p,boundary::west, condition_type::dirichlet, 0, 0, false, d);
                    }
                }
                if ( (bbox(0,1) - pbbox(0,1)) ==0 )
                {
                    gsInfo<<"Patch "<<p<<" is a boundary patch on east side!\n";
                    for (index_t d = 0; d!=3; d++)
                    {
                        bc.addCondition(p,boundary::east, condition_type::dirichlet, 0, 0, false, d);
                    }
                }
                if ( (bbox(1,0) - pbbox(1,0)) ==0 )
                {
                    gsInfo<<"Patch "<<p<<" is a boundary patch on south side!\n";
                    for (index_t d = 0; d!=3; d++)
                    {
                        bc.addCondition(p,boundary::south, condition_type::dirichlet, 0, 0, false, d);
                    }
                }
                if ( (bbox(1,1) - pbbox(1,1)) ==0 )
                {
                    gsInfo<<"Patch "<<p<<" is a boundary patch on north side!\n";
                    for (index_t d = 0; d!=3; d++)
                    {
                        bc.addCondition(p,boundary::north, condition_type::dirichlet, 0, 0, false, d);
                    }
                }
                // if ( (bbox(0,0) - pbbox(0,0)) * (bbox(0,1) - pbbox(0,1)) * (bbox(0,1) - pbbox(0,1)) * (bbox(1,1) - pbbox(1,1)) != 0 )
                //     boundaryPatches[p] = 0;
            }
        }
        else if(geometry > 999)
        {
            input = "planar/geometries/g" + std::to_string(geometry) + ".xml";
            gsReadFile<>(input, mp);
            mp.clearTopology();
            mp.computeTopology();
            for (gsMultiPatch<>::const_biterator bit = mp.bBegin(); bit != mp.bEnd(); ++bit)
                for (index_t d = 0; d!=3; d++)
                    bc.addCondition(bit->patch, bit->side(), condition_type::dirichlet, 0, 0, false, d);

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
        // TODO
        gsMultiBasis<real_t> multiBasis(mp);
        // h-refine
        for (int r =0; r < numRefine; ++r)
            mp.uniformRefine(1,multiBasis.maxCwiseDegree()-1);

        numRefine = 0;
    }

    gsWriteParaview(mp,"mp",1000,true,false);
    for (size_t p = 0; p!=mp.nPatches(); ++p)
        gsDebugVar(mp.patch(p));

    real_t thickness = 1.0;
    real_t E_modulus = 1.0;
    real_t PoissonRatio = 0.0;
    gsFunctionExpr<> force("0","0","1",3);
    gsFunctionExpr<> t(std::to_string(thickness),3);
    gsFunctionExpr<> E(std::to_string(E_modulus),3);
    gsFunctionExpr<> nu(std::to_string(PoissonRatio),3);

    gsPointLoads<real_t> pLoads = gsPointLoads<real_t>();
    if (geometry==1)
    {
        gsVector<> point(2); point<< 1, 0.5 ;
        gsVector<> load (3); load << 0.0, 0.0, 1.0 ;
        pLoads.addLoad(point, load, 0 );
        force = gsFunctionExpr<>("0","0","0",3);
    }


    std::vector<gsFunction<>*> parameters(2);
    parameters[0] = &E;
    parameters[1] = &nu;

    gsMaterialMatrixLinear<3,real_t> materialMatrix(mp,t,parameters);

    gsThinShellAssemblerBase<real_t> * assembler;

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
    gsWriteParaview(mp.basis(0),"basis",1000);

    for( index_t r = 0; r<=numRefine; ++r)
    {

        gsInfo<<"--------------------------------------------------------------\n";
        time.restart();
        if (smoothing==0)
        {
            // gsCompositeIncrSmoothnessGeom<2,real_t> cgeom(mp,3);
            // gsCompositeIncrSmoothnessBasis<2,real_t> & basis = dynamic_cast<gsCompositeIncrSmoothnessBasis<2,real_t>&>(cgeom.getCompBasis());

            // global2local = basis.getMapper().asMatrix();
            // geom = cgeom.exportToPatches();
        }
        else if (smoothing==1)
        {
            gsDPatch<2,real_t> dpatch(mp);
            dpatch.matrix_into(global2local);
            global2local = global2local.transpose();
            geom = dpatch.exportToPatches();
            dbasis = dpatch.localBasis();
        }
        else if (smoothing==2) // Pascal
        {
            mp.embed(2);
            gsApproxC1Spline<2,real_t> approxC1(mp,dbasis);
            approxC1.options().setSwitch("info",info);
            approxC1.options().setSwitch("plot",plot);
            // approxC1.options().setInt("gluingDataDegree",)
            // approxC1.options().setInt("gluingDataRegularity",)

            gsDebugVar(approxC1.options());

            approxC1.init();
            approxC1.compute();
            mp.embed(3);

            global2local = approxC1.getSystem();
            global2local = global2local.transpose();
            global2local.pruned(1,1e-10);
            geom = mp;
            approxC1.getMultiBasis(dbasis);
        }
        else if (smoothing==3) // Andrea
        {
            gsC1SurfSpline<2,real_t> smoothC1(mp,dbasis);
            smoothC1.init();
            smoothC1.compute();

            global2local = smoothC1.getSystem();
            global2local = global2local.transpose();
            smoothC1.getMultiBasis(dbasis);
        }
        else
            GISMO_ERROR("Option "<<smoothing<<" for smoothing does not exist");

        gsInfo<<"\tAssembly of mapping:\t"<<time.stop()<<"\t[s]\n";

        if (writeMatrix)
        {
            gsWrite(global2local,"mat");
            //gsWrite(geom,"geom");
            //gsWrite(dbasis,"dbasis");
        }

        bb2.init(dbasis,global2local);
        // gsMappedSpline<2,real_t> mspline(bb2,coefs);
        // geom = mspline.exportToPatches();

        //mem-leak here
        assembler = new gsThinShellAssembler<3, real_t, true>(geom,dbasis,bc,force,&materialMatrix);
        if (smoothing==1)
            assembler->options().setInt("Continuity",-1);
        else if (smoothing==2)
            assembler->options().setInt("Continuity",-1);
        assembler->setSpaceBasis(bb2);
        assembler->setPointLoads(pLoads);
        // gsOptionList options = assembler->options();
        // options.setInt("Continuity",1);
        // assembler->setOptions(options);

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

        l2err[r]= 0;
        h1err[r]= 0;
        linferr[r] = 0;


        gsInfo<<"\tError computations:\t"<<time.stop()<<"\t[s]\n"; // This takes longer for the D-patch, probably because there are a lot of points being evaluated, all containing the linear combinations of the MSplines

        // TO DO: Refine the mb
        mp.uniformRefine();
        dbasis = gsMultiBasis<>(mp);
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
        /// Make a gsMappedSpline to represent the solution
        // 1. Get all the coefficients (including the ones from the eliminated BCs.)
        gsMatrix<real_t> solFull = assembler->fullSolutionVector(solVector);

        // 2. Reshape all the coefficients to a Nx3 matrix
        GISMO_ASSERT(solFull.rows() % 3==0,"Rows of the solution vector does not match the number of control points");
        solFull.resize(solFull.rows()/3,3);

        // 3. Make the mapped spline
        gsMappedSpline<2,real_t> mspline(bb2,solFull);

        gsFunctionSum<real_t> def(&mp,&mspline);

        // 4. Plot the mapped spline on the original geometry
        gsField<> solField(mp, mspline,true);
        gsInfo<<"Plotting in Paraview...\n";
        gsWriteParaview<>( solField, "Deformation", 1000, false);

        // 5. Plot the mapped spline on the deformed geometry
        gsField<> defField(mp, def,true);
        gsInfo<<"Plotting in Paraview...\n";
        gsWriteParaview<>( defField, "mp_def", 1000, true);

        // QUASI INTERPOLATION
        // /*

            // gsMultiPatch<> mpatches = mbasis.exportToPatches(tmp);
            gsMultiPatch<> mp2;
            for (size_t p = 0; p!=mp.nPatches(); p++)
            {
                gsMatrix<> coefs;
                gsQuasiInterpolate<real_t>::localIntpl(mp.basis(p), mspline.piece(p), coefs);
                mp2.addPatch(mp.basis(p).makeGeometry( give(coefs) ));
            }

            gsField<> solfield(mp,mp2,true);
            gsWriteParaview(solfield,"solfield");

        // */


        /*

        // L2 PROJECTION


        gsMultiBasis<> mb(mp);
        gsBoundaryConditions<> bc_empty;

        typedef gsExprAssembler<>::geometryMap geometryMap;
        typedef gsExprAssembler<>::space       space;
        typedef gsExprAssembler<>::solution    solution;
        gsExprAssembler<> L2Projector(1,1);
        geometryMap G   = L2Projector.getMap(mp);

        L2Projector.setIntegrationElements(mb);
        space v = L2Projector.getSpace(bb2, 1);//m-splines
        space u = L2Projector.getTestSpace(v,mb);//TP splines
        //solution sol = L2Projector.getSolution(v,solFull);
        auto sol = L2Projector.getCoeff(mspline);

        u.setup(bc_empty,dirichlet::homogeneous);
        v.setup(bc_empty,dirichlet::homogeneous);


        gsExprEvaluator<> ev(L2Projector);
        gsMatrix<> pt(2,1);
        pt.setConstant(0.25);
        ev.writeParaview(sol,G,"solution");

        L2Projector.initSystem(3);
        L2Projector.assemble(u * v.tr(), u * sol.tr() );
        gsMatrix<> result = L2Projector.matrix().toDense().
            colPivHouseholderQr().solve(L2Projector.rhs());

        gsMultiPatch<> mp_res;
        gsMatrix<> coefs;
        index_t offset = 0;
        index_t blocksize = 0;
        for (index_t p = 0; p != mp.nPatches(); p++)
        {
            blocksize = mp.patch(p).coefs().rows();
            gsDebugVar(blocksize);
            gsDebugVar(offset);
            gsDebugVar(result.rows());
            gsDebugVar(mb.basis(p).size());
            coefs = result.block(offset,0,blocksize,result.cols());
            mp_res.addPatch(mb.basis(p).makeGeometry(give(coefs)));
            offset += blocksize;
        }

        gsField<> solfield_L2(mp,mp_res,true);
        gsWriteParaview(solfield_L2,"solfield_L2");

        // gsSparseSolver<>::QR solver( L2Projector.matrix() );
        // gsMatrix<> result = solver.solve(L2Projector.rhs().col(k));

        */

        //

        // Interpolate at anchors
        /*
        // gsMultiBasis<> mb(mp);
        gsMultiPatch<> mp2;
        for (size_t p = 0; p!=mp.nPatches(); ++p)
        {
            gsMatrix<> anchors;
            mb.basis(p).anchors_into(anchors);
            gsMatrix<> result;
            bb2.eval_into(p,anchors,result);
            mp2.addPatch(mb.basis(p).interpolateAtAnchors(result));
        }

        gsField<> solField(mp,mp2);

        gsWriteParaview(solField,"beer");
        */
        // gsDebugVar(solFull.rows());
        // gsDebugVar(mbasis.size());

        // gsMatrix<> u(2,1);
        // u.setConstant(0.25);

        // gsMatrix<> B ;
        // gsMatrix<index_t> actives;

        // mspline.basis().eval_into(u,B);
        // mbasis.active_into(u,actives);


        // gsField<> solField(mp.patch(0),mspline,true);

        // gsWriteParaview(solField,"mspline");



        return 0;

        // solFull.resize(solFull.rows()/3,3);
        // gsDebugVar(solFull);


        index_t compSize = solFull.rows() / 3;
        for (size_t d = 0; d!=3; d++)
        {
            gsDebugVar(global2local.rows());
            gsDebugVar(global2local.cols());
            gsDebugVar(coefs.rows());
            gsDebugVar(coefs.cols());

            gsMatrix<> coefs = solFull.block(d*compSize,0,compSize,1);
            global2local = coefs.asDiagonal() * global2local;

            gsMappedBasis<2,real_t> mbasis(dbasis,global2local);
            gsMappedSpline<2,real_t> mspline(mbasis,coefs);
            gsMultiPatch<> test = mspline.exportToPatches();
            gsDebugVar(test);
        }

        gsDebugVar(global2local.rows());
        gsDebugVar(global2local.cols());
        global2local = solFull.asDiagonal() * global2local;
        gsDebugVar(global2local);

        // gsMappedSpline<2,real_t> mspline(dbasis,global2local);

        // gsMultiPatch<> test = mspline.exportToPatches();
        // gsWriteParaview("test",test,1000,true);




        // gsMultiBasis<> multiBasis(mp);
        // for (size_t np = 0; np < mp.nPatches(); np++)
        // {
        //     gsBasis<> & basis = multiBasis.basis(np);
        //     gsMatrix<> points2D = basis.anchors();
        //     // Update dbasis with solution
        //     typename gsGeometry<>::uPtr patch = basis.interpolateAtAnchors(dbasis.basis(np).eval(points2D));
        // }

        // assembler->plotSolution("solution", solVector);

        // gsMultiPatch<> deformation = assembler->constructDisplacement(solVector);
        // // gsMultiPatch<> mp_def = assembler->constructSolution(solVector);

        // gsField<> solField(geom, deformation);
        // gsInfo<<"Plotting in Paraview...\n";
        // gsWriteParaview<>( solField, "Deformation", 1000, true);

        // deformation = assembler->constructSolution(solVector);
        // gsWriteParaview<>( deformation, "deformed_geom", 1000, true);


        // gsWriteParaview( geom, "geom",100,true);
    }
    //! [Export visualization in ParaView]

    return EXIT_SUCCESS;


}
