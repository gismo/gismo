/** @file mortaring_example.cpp

    @brief Provides an example for mortar.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): S. Takacs
*/

#include <ctime>

#include <gismo.h>
#include <gsAssembler/gsJumpAssembler.h>

using namespace gismo;


int main(int argc, char *argv[])
{
    /************** Define command line options *************/

    std::string geometry("domain2d/annulus_dg2.xml");
    index_t refinements = 3;
    index_t degree = 2;
    bool nitsche = false;
    std::string boundaryConditions("d");
    std::string out;
    bool plot = false;

    gsCmdLine cmd("Solves a PDE with an isogeometric discretization using a multigrid solver.");
    cmd.addString("g", "Geometry",              "Geometry file", geometry);
    cmd.addInt   ("r", "Refinements",           "Number of uniform h-refinement steps to perform before solving", refinements);
    cmd.addInt   ("p", "Degree",                "Degree of the B-spline discretization space", degree);
    cmd.addSwitch("",  "Nitsche",               "Use Nitsche method for Dirichlet boundary conditions", nitsche);
    cmd.addString("b", "BoundaryConditions",    "Boundary conditions", boundaryConditions);
    cmd.addString("",  "out",                   "Write solution and used options to file", out);
    cmd.addSwitch(     "plot",                  "Plot the result with Paraview", plot);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    gsInfo << "Run mortaring_example with options:\n" << cmd << "\n";

    /******************* Define geometry ********************/

    gsInfo << "Define geometry... " << std::flush;

    //! [Define Geometry]
    gsMultiPatch<>::uPtr mpPtr = gsReadFile<>(geometry);
    //! [Define Geometry]
    if (!mpPtr)
    {
        gsInfo << "No geometry found in file " << geometry << ".\n";
        return EXIT_FAILURE;
    }
    gsMultiPatch<>& mp = *mpPtr;

    //! [Define Geometry2]
    //for (index_t i=0; i<splitPatches; ++i)
    //{
    //    gsInfo << "split patches uniformly... " << std::flush;
    //    mp = mp.uniformSplit();
    //}
    //! [Define Geometry2]

    gsInfo << "done.\n";

    /************** Define boundary conditions **************/

    gsInfo << "Define boundary conditions... " << std::flush;

    //! [Define Source]
    // Right-hand-side
    gsFunctionExpr<> f( "2*sin(x)*cos(y)", mp.geoDim() );

    // Dirichlet function
    gsFunctionExpr<> gD( "sin(x)*cos(y)", mp.geoDim() );

    // Neumann
    gsConstantFunction<> gN( 1.0, mp.geoDim() );

    gsBoundaryConditions<> bc;
    //! [Define Source]
    {
        const index_t len = boundaryConditions.length();
        index_t i = 0;
        for (gsMultiPatch<>::const_biterator it = mp.bBegin(); it < mp.bEnd(); ++it)
        {
            char b_local;
            if ( len == 1 )
                b_local = boundaryConditions[0];
            else if ( i < len )
                b_local = boundaryConditions[i];
            else
            {
                gsInfo << "\nNot enough boundary conditions given.\n";
                return EXIT_FAILURE;
            }

            if ( b_local == 'd' )
                bc.addCondition( *it, condition_type::dirichlet, &gD );
            else if ( b_local == 'n' )
                bc.addCondition( *it, condition_type::neumann, &gN );
            else
            {
                gsInfo << "\nInvalid boundary condition given; only 'd' (Dirichlet) and 'n' (Neumann) are supported.\n";
                return EXIT_FAILURE;
            }

            ++i;
        }
        if ( len > i )
            gsInfo << "\nToo many boundary conditions have been specified. Ignoring the remaining ones.\n";
        gsInfo << "done. "<<i<<" boundary conditions set.\n";
    }


    /************ Setup bases and adjust degree *************/

    //! [Define Basis]
    gsMultiBasis<> mb(mp);
    //! [Define Basis]

    gsInfo << "Setup bases and adjust degree... " << std::flush;

    //! [Set degree and refine]
    for ( size_t i = 0; i < mb.nBases(); ++ i )
        mb[i].setDegreePreservingMultiplicity(degree);

    for ( index_t i = 0; i < refinements; ++i )
        mb.uniformRefine();
    //! [Set degree and refine]

    gsInfo << "done.\n";

    /********* Setup assembler and assemble matrix **********/

    gsInfo << "Setup assembler and assemble matrix... " << std::flush;

    //! [Assemble]
    gsPoissonAssembler<> assembler(
        mp,
        mb,
        bc,
        f,
        nitsche ? dirichlet::nitsche : dirichlet::elimination,
        iFace::none // iFace::dg
    );
    assembler.assemble();
    //! [Assemble]

    gsSparseMatrix<> mat = assembler.matrix();

    std::vector<gsSparseMatrix<>> constraints;
    index_t nConstraints = 0;
    for ( typename gsMultiPatch<>::const_iiterator it = mp.iBegin(); it != mp.iEnd(); ++it )
    {
        gsJumpAssembler<real_t> ifAssembler(
            mp,
            mb,
            bc,
            f,
            nitsche ? dirichlet::nitsche : dirichlet::elimination,
            it
        );
        ifAssembler.assemble();
        //mat += ifAssembler.matrix(); //// TODO
        
        // TODO: assume all patches to have same number of dofs!
        size_t dofsPerPatch = mat.rows()/8;
        size_t k1 = it->first().patch;
        gsSparseMatrix<> tmp = ifAssembler.matrix().middleRows(k1*dofsPerPatch,dofsPerPatch);
        std::vector<index_t> rowsToBeChosen;
        for (index_t i=0; i<tmp.rows(); ++i)
            for (index_t j=0; j<tmp.cols(); ++j)
                if (tmp(i,j)!=0)
                {
                    rowsToBeChosen.push_back(i);
                    break;
                }
        gsSparseMatrix<> chser(rowsToBeChosen.size(), tmp.rows());
        for (size_t i=0; i<rowsToBeChosen.size(); ++i)
            chser(i, rowsToBeChosen[i]) = 1;
        
        //gsInfo << "chser " << chser.rows() << "x" << chser.cols() << "\n";
        //gsInfo << "tmp " << tmp.rows() << "x" << tmp.cols() << "\n";
        //gsInfo << "ifAssembler.matrix()  " << ifAssembler.matrix() .rows() << "x" << ifAssembler.matrix() .cols() << "\n";
        
        gsSparseMatrix<> loc = chser * tmp.middleCols(k1*dofsPerPatch,dofsPerPatch) * chser.transpose();
        gsInfo << "loc=" << loc << "\n\n";
        gsMatrix<> locinv;
        makeSparseCholeskySolver(loc)->toMatrix(locinv);
        gsSparseMatrix<> locinv2(locinv.rows(), locinv.cols());
        for (index_t i=0; i<locinv.rows(); ++i)
            for (index_t j=0; j<locinv.cols(); ++j)
                locinv2(i,j) = locinv(i,j);
        
        constraints.push_back(locinv2 * chser * tmp);
        nConstraints += rowsToBeChosen.size();
        
    }
        
    gsSparseMatrix<> allConstraints(nConstraints,assembler.matrix().cols());
    {
        index_t r = 0;
        for (size_t i=0; i<constraints.size(); ++i)
            for (index_t j=0; j<constraints[i].rows(); ++j)
            {
                for (index_t c=0; c<constraints[i].cols(); ++c)
                    if (constraints[i](j,c)*constraints[i](j,c)>1e-6)
                        allConstraints(r,c) = constraints[i](j,c);
                ++r;
            }
    }
    
    gsInfo << "Constraints\n" << allConstraints << "\n\n";
    
    gsSparseMatrix<> blockMatrix(assembler.matrix().rows()+nConstraints,assembler.matrix().cols()+nConstraints);
    for (index_t i=0; i<assembler.matrix().rows(); ++i)
        for (index_t j=0; j<assembler.matrix().cols(); ++j)
            if (assembler.matrix()(i,j)!=0)
                blockMatrix(i,j) = assembler.matrix()(i,j);
    for (index_t i=0; i<nConstraints; ++i)
        for (index_t j=0; j<assembler.matrix().cols(); ++j)
            if (allConstraints(i,j)!=0)
            {
                blockMatrix(assembler.matrix().rows()+i,j) = allConstraints(i,j);
                blockMatrix(j,assembler.matrix().cols()+i) = allConstraints(i,j);
            }
    
    
    gsMatrix<> rhs;
    rhs.setZero(blockMatrix.rows(), 1);
    rhs.topRows(assembler.matrix().rows()) = assembler.rhs();
    
    gsInfo << "done: "<<assembler.matrix().rows()<<" dofs.\n";

    /**************** Setup solver and solve ****************/

    gsInfo << "Setup solver and solve... " << std::flush;

    //! [Solve]
    gsMatrix<> xx;
    makeSparseLUSolver(blockMatrix)->apply(rhs, xx);
    gsMatrix<> x = xx.topRows(assembler.matrix().rows());
    //! [Solve]

    gsInfo << "done.\n\n";

    /******************** Print end Exit ********************/

    if (!out.empty())
    {
        gsFileData<> fd;
        std::time_t time = std::time(NULL);
        fd.add(cmd);
        fd.add(x);
        fd.addComment(std::string("mortaring_example   Timestamp:")+std::ctime(&time));
        fd.save(out);
        gsInfo << "Write solution to file " << out << "\n";
    }

    if (plot)
    {
        // Construct the solution as a scalar field
        gsMultiPatch<> mpsol;
        assembler.constructSolution(x, mpsol);
        gsField<> sol( assembler.patches(), mpsol );

        // Write solution to paraview files
        gsInfo << "Write Paraview data to file mortaring_result.pvd\n";
        gsWriteParaview<>(sol, "mortaring_result", 1000);
        //gsFileManager::open("mortaring_result.pvd");
    }
    if (!plot&&out.empty())
    {
        gsInfo << "Done. No output created, re-run with --plot to get a ParaView "
                  "file containing the solution or --out to write solution to xml file.\n";
    }
    return EXIT_SUCCESS;
}


