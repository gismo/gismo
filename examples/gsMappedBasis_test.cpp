#/** @file gsMappedBasis_test.cpp

    @brief Example using the gsMappedBasis class

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris
*/

//! [Include namespace]
#include <gismo.h>

using namespace gismo;
//! [Include namespace]

#include <cfenv> //for float exceptions

// Fills in ones on the diagonal of zero-columns of \a mat
void completeWithIdentity(gsSparseMatrix<> & mat)
{
    for (index_t k = 0; k!=mat.outerSize(); ++k)
        if ( 0 == mat.innerVector(k).nonZeros() )
            mat.coeffRef(k,k) = 1;
}

int main(int argc, char *argv[])
{
    // Enable floating point exceptions
    feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);

    //! [Parse command line]
    bool plot = false;
    int numRefine  = 3;
    int qRule  = 1;

    gsFileManager::addSearchPaths(gsFileManager::getHomePath() + "/dbox/Florian_ThreePatches/star");
    std::string fn("star3_");

    gsCmdLine cmd("Example using mapped spline bases.");
    cmd.addInt( "r", "uniformRefine", "Number of Uniform h-refinement steps to perform before solving",  numRefine );
    cmd.addString( "f", "file", "Input XML file prefix", fn );
    cmd.addSwitch("plot", "Create a ParaView visualization file with the solution", plot);
    cmd.addInt( "q", "quadrature", "1 or 2", qRule );

    cmd.getValues(argc,argv);
    //! [Parse command line]

    //! [Problem setup]
    int init_knots = 0;
    gsFileData<> fd;
    gsMultiPatch<> mp;
    gsMultiBasis<> mb;
    gsSparseMatrix<> cf;

    gsExprAssembler<> A(1,1);
    gsExprEvaluator<> ev(A);
    //ev.options().setInt("quB", 2);// increase norm quadrature accuracy
    //ev.options().addInt("quRule","",2); // 1:Gauss-Legendre, 2:Gauss-Lobatto
    //A.options().addInt("quRule","", qRule); // 1:Gauss-Legendre, 2:Gauss-Lobatto
    //A.options().addInt("quB","", 1);
    A.setIntegrationElements(mb);

    //gsInfo<<"Active options:\n"<< A.options() <<"\n";
    typedef gsExprAssembler<>::geometryMap geometryMap;
    typedef gsExprAssembler<>::variable    variable;
    typedef gsExprAssembler<>::space       space;
    typedef gsExprAssembler<>::solution    solution;

    // Set the discretization space
    short_t d = 2; //mb.domainDim();
    gsMappedBasis<2,real_t> bb2;
    gsMappedBasis<3,real_t> bb3;
    space u = (2==d ? A.getSpace(bb2) : A.getSpace(bb3) );

    // Get the map and function
    geometryMap G = A.getMap(mp);
    gsFunctionExpr<> ff;
    ff = gsFunctionExpr<>("cos(pi*x/8) * sin(pi*y/8)", 2);
    variable f = A.getCoeff(ff, G);

    // Solution vector and solution variable
    gsMatrix<> sVector;
    solution s = A.getSolution(u, sVector);

    gsSparseSolver<>::CGDiagonal solver;

    //! [Solver loop]
    gsVector<> l2err(numRefine+1), h1err(numRefine+1), linferr(numRefine+1),
        b2err(numRefine+1), b1err(numRefine+1), binferr(numRefine+1);
    gsVector<index_t>  CGiter(numRefine+1);
    gsInfo<< "(dot1=assembled, dot2=solved, dot3=got_error)\n"
        "\nDoFs: ";

    for( index_t r = 0; r<=numRefine; ++r)
    {
        fd.read(fn + util::to_string(init_knots++) + ".xml.gz");
        fd.getFirst(mp);
        fd.getFirst(mb);
        fd.getFirst(cf);
        completeWithIdentity(cf);
        if (2==d)
            bb2.init(mb,cf);
        if (3==d)
            bb3.init(mb,cf);
            

        //u.setup(bc, dirichlet::interpolation, 0);
        
        // Initialize the system
        A.initSystem();

        gsInfo<< A.numDofs() <<std::flush;

        // Compute the system matrix and right-hand side
        A.assemble(  u * u.tr(), u * f );

        // Poisson
//        A.assemble( igrad(u,G) * igrad(u,G).tr() * meas(G), - u * ilapl(f) * meas(G) );

        // Biharmonic
        //A.assemble(  ilapl(u,G) * ilapl(u,G).tr() * meas(G), - u * ilapl(f) * meas(G) );

        gsInfo<< "." <<std::flush;// Assemblying done

        solver.compute( A.matrix() );
        sVector = solver.solve(A.rhs());
        CGiter[r] = solver.iterations();

        gsInfo<< "." <<std::flush; // Linear solving done

        l2err[r]= math::sqrt( ev.integral( (f - s).sqNorm()*meas(G) ) / ev.integral(f.sqNorm()*meas(G)) );

        h1err[r]= l2err[r] + math::sqrt(ev.integral( ( igrad(f) - grad(s)*jac(G).inv() ).sqNorm()*meas(G) )/ev.integral( igrad(f).sqNorm()*meas(G) ) );
        
        linferr[r] = ev.max( f-s ) / ev.max(f);
        
        gsInfo<< ". " <<std::flush; // Error computations done
    }
    //! [Solver loop]

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
        ev.options().setInt   ("plot.npts"    , 100);
        ev.writeParaview( s, G, "solution");
        //ev.writeParaview( grad(s), G, "solution_grad");
        ev.writeParaview( f, G, "solution_ex");
        //ev.writeParaview( grad(f), G, "solution_ex_grad");
        ev.writeParaview( (f-s), G, "error_pointwise");
        gsFileManager::open("error_pointwise.pvd");
    }
    //! [Export visualization in ParaView]

    return EXIT_SUCCESS;

}// end main
