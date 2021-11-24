/** @file gsSmoothBasis_test.h

    @brief File testing the various gsUnstructuredSpline classes .

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): F. Buchegger
*/

#include <gismo.h>

#include <gsUnstructuredSplines/gsMPBESBasis.h>
#include <gsUnstructuredSplines/gsMPBESGeom.h>

#include <gsUnstructuredSplines/gsDPatch.h>

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

            // Needed
            mp.uniformRefine(1);
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

        }
        else if (geometry==4)
        {
            gsReadFile<>("planar/hexagon_3p.xml", mp);
            // mp.embed(3);
            mp.clearTopology();
            mp.computeTopology();
        }
        else if (geometry==5)
        {
            gsReadFile<>("planar/hexagon_5p.xml", mp);
            // mp.embed(3);
            mp.clearTopology();
            mp.computeTopology();

            // Needed because of the multiple EVs on one patch
            mp.uniformRefine();
        }
        else if (geometry==6)
        {
            gsReadFile<>("planar/triangle2d_u3p.xml", mp);
            // mp.embed(3);
            mp.clearTopology();
            mp.computeTopology();

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

    gsMultiBasis<> mb(mp);
    gsExprAssembler<> A(1,1);
    gsExprEvaluator<> ev(A);
    A.setIntegrationElements(mb);

    //gsInfo<<"Active options:\n"<< A.options() <<"\n";
    typedef gsExprAssembler<>::geometryMap geometryMap;
    typedef gsExprAssembler<>::space       space;
    typedef gsExprAssembler<>::solution    solution;

    // Set the discretization space
    gsMappedBasis<2,real_t> bb2;
    space u = A.getSpace(bb2);

    // Get the map and function
    geometryMap G = A.getMap(geom);
    gsFunctionExpr<> ff;
    ff = gsFunctionExpr<>("cos(pi*x/8) * sin(pi*y/8)", 2);
    auto f = A.getCoeff(ff, G);

    // Solution vector and solution variable
    gsMatrix<> sVector;
    solution s = A.getSolution(u, sVector);

    gsSparseSolver<>::CGDiagonal solver;

    //! [Solver loop]
    gsVector<> l2err(numRefine+1), h1err(numRefine+1), linferr(numRefine+1),
        b2err(numRefine+1), b1err(numRefine+1), binferr(numRefine+1);
    gsVector<index_t>  CGiter(numRefine+1);

    gsSparseMatrix<> global2local;
    gsMatrix<> coefs;
    gsStopwatch time;
    for( index_t r = 0; r<=numRefine; ++r)
    {
        gsInfo<<"--------------------------------------------------------------\n";
        time.restart();
        if (smoothing==0)
        {
            gsMPBESSpline<2,real_t> cgeom(mp,3);
            gsMappedBasis<2,real_t> basis = cgeom.getMappedBasis();

            global2local = basis.getMapper().asMatrix();
            geom = cgeom.exportToPatches();
            auto container = basis.getBasesCopy();
            mb = gsMultiBasis<>(container,mp.topology());
        }
        else if (smoothing==1)
        {
            gsDPatch<2,real_t> dpatch(mp);
            dpatch.matrix_into(global2local);
            global2local = global2local.transpose();
            geom = dpatch.exportToPatches();
            mb = dpatch.localBasis();
        }
        else
            GISMO_ERROR("Option "<<smoothing<<" for smoothing does not exist");

        gsInfo<<"\tAssembly of mapping:\t"<<time.stop()<<"\t[s]\n";

        bb2.init(mb,global2local);
        // gsMappedSpline<2,real_t> mspline(bb2,coefs);
        // geom = mspline.exportToPatches();


        u.setup(bc, dirichlet::interpolation, -1);

        time.restart();
        // Initialize the system
        A.initSystem();

        gsInfo<<"\tDegrees of freedom:\t"<< A.numDofs() <<"\n";

        // Compute the system matrix and right-hand side
        A.assemble(  u * u.tr(), u * f);

        // Poisson
//        A.assemble( igrad(u,G) * igrad(u,G).tr() * meas(G), - u * ilapl(f) * meas(G) );

        // Biharmonic
        //A.assemble(  ilapl(u,G) * ilapl(u,G).tr() * meas(G), - u * ilapl(f) * meas(G) );

        gsInfo<<"\tSystem assembly:\t"<<time.stop()<<"\t[s]\n";

        time.restart();

        solver.compute( A.matrix() );
        sVector = solver.solve(A.rhs());
        CGiter[r] = solver.iterations();

        gsInfo<<"\tSolving system:\t\t"<<time.stop()<<"\t[s]\n";
        time.restart();

        l2err[r]= math::sqrt( ev.integral( (f - s).sqNorm()*meas(G) ) / ev.integral(f.sqNorm()*meas(G)) );

        h1err[r]= l2err[r] + math::sqrt(ev.integral( ( igrad(f) - grad(s)*jac(G).inv() ).sqNorm()*meas(G) )/ev.integral( igrad(f).sqNorm()*meas(G) ) );

        linferr[r] = ev.max( f-s ) / ev.max(f);

        gsInfo<<"\tError computations:\t"<<time.stop()<<"\t[s]\n"; // This takes longer for the D-patch, probably because there are a lot of points being evaluated, all containing the linear combinations of the MSplines

        // TO DO: Refine the mb
        mp.uniformRefine();
        mb = gsMultiBasis<>(mp);
    }
    //! [Solver loop]

    gsDofMapper mapper = u.mapper();
    gsMatrix<> m_ddofs = u.fixedPart();
    const index_t dim = u.dim();
    gsMatrix<> cc(bb2.size(),dim);
    cc.setZero();
    for ( size_t p =0; p!=mapper.numPatches(); ++p) // Deform the geometry
    {
        for (index_t c = 0; c!=dim; c++) // for all components
        {
            gsDebugVar(p);
            gsDebugVar(c);
            gsDebugVar(bb2.size(p));
            gsDebugVar(bb2.localSize(p));
            gsDebugVar(mapper.size(p));
            gsDebugVar(mapper.patchSize(p,c));
            // loop over all basis functions (even the eliminated ones)
            for (index_t i = 0; i < mapper.patchSize(p,c); ++i)
            {
                const int ii = mapper.index(i, p, c);
                if ( mapper.is_free_index(ii) ) // DoF value is in the solVector
                    cc(i,c) = sVector.at(ii);
                else // eliminated DoF: fill with Dirichlet data
                {
                    cc(i,c) =  m_ddofs.at( mapper.global_to_bindex(ii) );
                }
            }

        }
        gsDebugVar(cc);
        gsMatrix<> tmp;

    }


    gsInfo<< "\n\nCG it.: "<<CGiter.transpose()<<"\n";

    //! [Error and convergence rates]
    gsInfo<<"\n* Error\n";
    gsInfo<< "H1    "<<std::scientific<<h1err.transpose()<<"\n";
    gsInfo<< "L2    "<<std::scientific<<std::setprecision(3)<<l2err.transpose()<<"\n";
    gsInfo<< "Linf  "<<std::scientific<<linferr.transpose()<<"\n";

    gsInfo<<"\n* EoC\n";
    gsInfo<< "H1c   0 "<< std::fixed<<std::setprecision(2)
          <<( h1err.head(numRefine).array() /
              h1err.tail(numRefine).array() ).log().transpose() / std::log(2.0) <<"\n";
    gsInfo<< "L2c   0 " << std::fixed<<std::setprecision(2)
          << ( l2err.head(numRefine).array() /
               l2err.tail(numRefine).array() ).log().transpose() / std::log(2.0) <<"\n";

    gsInfo<<   "Linfc 0 "<< std::fixed<<std::setprecision(2)
          <<( linferr.head(numRefine).array() /
              linferr.tail(numRefine).array() ).log().transpose() / std::log(2.0) <<"\n";

    //! [Error and convergence rates]

    //! [Export visualization in ParaView]
    if (plot)
    {
        gsInfo<<"Plotting in Paraview...\n";
        ev.options().setSwitch("plot.elements", false);
        ev.options().setInt   ("plot.npts"    , 500);
        ev.writeParaview( s, G, "solution");
        //ev.writeParaview( grad(s), G, "solution_grad");
        ev.writeParaview( f, G, "solution_ex");
        //ev.writeParaview( grad(f), G, "solution_ex_grad");
        ev.writeParaview( (f-s), G, "error_pointwise");
        gsWriteParaview( geom, "geom",100,true);
    }
    //! [Export visualization in ParaView]

    return EXIT_SUCCESS;


}
