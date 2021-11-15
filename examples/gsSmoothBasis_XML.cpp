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
// #include <gsUnstructuredSplines/gsC1Spline.h>

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

    std::string fn1,fn2,fn3,fn4;
    fn1 = "pde/1p_square_geom.xml";
    fn2 = "pde/1p_square_bvp.xml";
    fn3 = "options/solver_options.xml";

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

    /*
        to do:
        - remove hard-coded IDs from XML reader
     */

    gsFileData<> fd;
    gsInfo<<"Reading geometry from "<<fn1<<"...";
    gsReadFile<>(fn1, mp);
    gsInfo<<"Finished\n";

    fd.read(fn2);
    // index_t num = 0;
    // gsInfo<<"Reading BCs from "<<fn2<<"...";
    // num = fd.template count<gsBoundaryConditions<>>();
    // GISMO_ENSURE(num==1,"Number of boundary condition objects in XML should be 1, but is "<<num);
    // fd.template getFirst<gsBoundaryConditions<>>(bc); // Multipatch domain
    // gsInfo<<"Finished\n";

    gsFunctionExpr<> hom("0",3);
    for (index_t d = 0; d!=3; d++)
    {
        bc.addCondition(0,boundary::north, condition_type::dirichlet, &hom, 0, false, d );
        bc.addCondition(0,boundary::south, condition_type::dirichlet, &hom, 0, false, d);
        bc.addCondition(0,boundary::west, condition_type::dirichlet, &hom, 0, false, d);
        bc.addCondition(0,boundary::east, condition_type::dirichlet, &hom, 0, false, d);
    }
    bc.setGeoMap(mp);

    gsDebugVar(bc);

    // Loads
    gsFunctionExpr<> force, pressure;
    gsInfo<<"Reading force function from "<<fn2<<" (ID=21) ...";
    fd.getId(21, force); // id=1: source function
    gsInfo<<"Finished\n";
    // fd.getId(22, pressure); // id=1: source function ------- TO DO!
    // gsInfo<<"Pressure function "<< force << "\n";

    // Loads
    gsPointLoads<real_t> pLoads = gsPointLoads<real_t>();
    gsMatrix<> points,loads;
    gsInfo<<"Reading point load locations from "<<fn2<<" (ID=30) ...";
    fd.getId(30,points);
    gsInfo<<"Finished\n";
    gsInfo<<"Reading point loads from "<<fn2<<" (ID=31) ...";
    fd.getId(31,loads);
    gsInfo<<"Finished\n";
    for (index_t k =0; k!=points.cols(); k++)
        pLoads.addLoad(points.col(k), loads.col(k), 0 ); // in parametric domain!

    // Material properties
    gsFunctionExpr<> t,E,nu,rho;
    gsInfo<<"Reading thickness from "<<fn2<<" (ID=10) ...";
    fd.getId(10,t);
    gsInfo<<"Finished\n";

    gsInfo<<"Reading Young's Modulus from "<<fn2<<" (ID=11) ...";
    fd.getId(11,E);
    gsInfo<<"Finished\n";

    gsInfo<<"Reading Poisson ratio from "<<fn2<<" (ID=12) ...";
    fd.getId(12,nu);
    gsInfo<<"Finished\n";

    gsInfo<<"Reading density from "<<fn2<<" (ID=13) ...";
    fd.getId(13,rho);
    gsInfo<<"Finished\n";

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
    for (size_t p = 0; p!=mp.nPatches(); ++p)
        gsDebugVar(mp.patch(p));

    std::vector<gsFunction<>*> parameters(2);
    parameters[0] = &E;
    parameters[1] = &nu;

    gsMaterialMatrixLinear<3,real_t> materialMatrix(mp,t,parameters,rho);

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
        if (r < numRefine)
        {
            mp.uniformRefine();
            dbasis = gsMultiBasis<>(mp);
        }
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
        gsField<> solField(geom, mspline,true);
        gsInfo<<"Plotting in Paraview...\n";
        gsWriteParaview<>( solField, "Deformation", 1000, false);

        // 5. Plot the mapped spline on the deformed geometry
        gsField<> defField(geom, def,true);
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
