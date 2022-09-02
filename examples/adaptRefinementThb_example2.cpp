/** @file adaptRefinementThb_example.cpp

    @brief Tutorial on how to use G+Smo to solve the Poisson equation,
    see the \ref PoissonTutorial

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): S. Kleiss
*/

//! [Include namespace]
# include <gismo.h>
# include <gsAssembler/gsAdaptiveRefUtils.h>
# include <gsAssembler/gsAdaptiveMeshing.h>

using namespace gismo;
//! [Include namespace]

template<typename T>
class gsElementErrorPlotter : public gsFunction<T>
{
public:
    gsElementErrorPlotter(const gsBasis<T>& mp, const std::vector<T>& errors ) : m_mp(mp),m_errors(errors)
    {

    }

    virtual void eval_into(const gsMatrix<T>& u, gsMatrix<T>& res) const
    {
        // Initialize domain element iterator -- using unknown 0
        res.setZero(1,u.cols());
        for(index_t i=0; i<u.cols();++i)
        {
            int iter =0;
            // Start iteration over elements

            typename gsBasis<T>::domainIter domIt = m_mp.makeDomainIterator();
            for (; domIt->good(); domIt->next() )
            {
                 bool flag = true;
                const gsVector<T>& low = domIt->lowerCorner();
                const gsVector<T>& upp = domIt->upperCorner();


                for(int d=0; d<domainDim();++d )
                {
                    if(low(d)> u(d,i) || u(d,i) > upp(d))
                    {
                        flag = false;
                        break;
                    }
                }
                if(flag)
                {
                     res(0,i) = m_errors.at(iter);
                     break;
                }
                iter++;
            }
        }
    }

    short_t domainDim() const { return m_mp.dim();}

private:
    const gsBasis<T>& m_mp;
    const std::vector<T>& m_errors;
};


int main(int argc, char *argv[])
{
    //! [Parse command line]
    bool plot = false;
    bool coarsening = false;
    bool admissible = false;

    index_t Ext = 0;
    index_t CExt = 0;

    index_t numElevate = 0;
    index_t numRefine = 0;

    index_t type = 2;

    // Number of refinement loops to be done
    int numRefinementLoops = 4;

    gsCmdLine cmd("Tutorial on solving a Poisson problem.");
    cmd.addSwitch("plot", "Create a ParaView visualization file with the solution", plot);
    cmd.addSwitch("coarsen", "Also coarsen", coarsening);
    cmd.addSwitch("admissible", "Admissible refinement", admissible);
    cmd.addInt("r", "numLoop", "number of refinement loops", numRefinementLoops);
    cmd.addInt("R", "uniformRefine", "Number of Uniform h-refinement loops",  numRefine );
    cmd.addInt("e", "degreeElevation","Number of degree elevation steps to perform before solving (0: equalize degree in all directions)", numElevate );
    cmd.addInt("E", "Extension", "Extension", Ext);
    cmd.addInt("C", "CExtension", "CExtension", CExt);
    cmd.addInt("T", "BasisType", "BasisType; 0 = BSpline, 1 = THB, 2 = HB", type);
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }
    //! [Parse command line]

    // --------------- specify exact solution and right-hand-side ---------------

    //! [Function data]
    // Define exact solution (will be used for specifying Dirichlet boundary conditions
    gsFunctionExpr<> g("if( y>0, ( (x^2+y^2)^(1.0/3.0) )*sin( (2*atan2(y,x) - pi)/3.0 ), ( (x^2+y^2)^(1.0/3.0) )*sin( (2*atan2(y,x)+3*pi)/3.0 ) )", 2);
    // Define source function
    gsFunctionExpr<> f("0",2);
    //! [Function data]

    // Print out source function and solution
    gsInfo<<"Source function " << f << "\n";
    gsInfo<<"Exact solution "  << g << "\n\n";


    // --------------- read geometry from file ---------------

    //! [GetGeometryData]
    // Read xml and create gsMultiPatch
    std::string fileSrc( "planar/lshape2d_3patches_thb.xml" );
    gsMultiPatch<real_t> patches;
    gsReadFile<real_t>( fileSrc, patches);
    //! [GetGeometryData]
    gsInfo << "The domain is a "<< patches <<"\n";

    //! [computeTopology]
    // Get all interfaces and boundaries:
    patches.computeTopology();
    //! [computeTopology]

    //! [GetGeometryDataTens]
    std::string fileSrcTens( "planar/lshape2d_3patches_tens.xml" );
    gsMultiPatch<real_t> patchesTens;
    gsReadFile<real_t>( fileSrcTens, patchesTens);
    patchesTens.computeTopology();
    //! [GetGeometryDataTens]


    // --------------- add bonudary conditions ---------------
    //! [Boundary conditions]
    gsBoundaryConditions<> bcInfo;

    // For simplicity, set Dirichlet boundary conditions
    // given by exact solution g on all boundaries:
    for ( gsMultiPatch<>::const_biterator
            bit = patches.bBegin(); bit != patches.bEnd(); ++bit)
    {
       bcInfo.addCondition( *bit, condition_type::dirichlet, &g );
    }
    //! [Boundary conditions]


    // --------------- set up basis ---------------

    //! [GetBasisFromTHB]
    // Copy basis from the geometry
    gsMultiBasis<> bases( patches );
    //! [GetBasisFromTHB]

    // Elevate and p-refine the basis to order p + numElevate
    // where p is the highest degree in the bases
    bases.setDegree( bases.maxCwiseDegree() + numElevate);

    // h-refine each basis
    for (int r =0; r < numRefine; ++r)
        bases.uniformRefine();

    //! [GetBasisFromTens]
    // Copy tensor basis
    gsMultiBasis<real_t> basesTens( patchesTens );

    // Create a "basisContainer"
    std::vector< gsBasis<real_t>* > basisContainer;

    // fill the "basisContainer" with patch-wise...
    for ( size_t i = 0; i < basesTens.nBases(); i++)
    {
        if (type==0)
            basisContainer.push_back(&basesTens.basis(i));
        else if (type==1)
            basisContainer.push_back(new gsTHBSplineBasis<2,real_t>( basesTens.basis(i) ));
        else if (type==2)
            basisContainer.push_back(new gsHBSplineBasis<2,real_t>( basesTens.basis(i) ));
        else
            GISMO_ERROR("Basis type unknown");
    }

    // finally, create the gsMultiBasis containing gsTHBSpline ...
    gsMultiBasis<real_t> basesFromTens( basisContainer, patches );
    //! [GetBasisFromTens]


    //! [initialRefinements]
    // Number of initial uniform refinement steps:
    int numInitUniformRefine  = 2;

    gsMultiPatch<> patches0 = patches;

    for (int i = 0; i < numInitUniformRefine; ++i)
    {
     bases.uniformRefine();
     patches.uniformRefine();
    }
    //! [initialRefinements]


    // --------------- set up adaptive refinement loop ---------------

    //! [adaptRefSettings]


    // Specify cell-marking strategy...
    MarkingStrategy adaptRefCrit = PUCA;
    //MarkingStrategy adaptRefCrit = GARU;
    //MarkingStrategy adaptRefCrit = errorFraction;

    // ... and parameter.
    const real_t adaptRefParam = 0.9;
    //! [adaptRefSettings]

    // --------------- adaptive refinement loop ---------------

    gsParaviewCollection collection("adaptRef");

    gsAdaptiveMeshing<real_t> mesher(bases);
    mesher.options().setInt("CoarsenRule",2);
    mesher.options().setInt("RefineRule",2);
    mesher.options().setReal("CoarsenParam",0.1);
    mesher.options().setReal("RefineParam",0.1);
    mesher.options().setInt("CoarsenExtension",CExt);
    mesher.options().setInt("RefineExtension",Ext);
    mesher.options().setSwitch("Admissible",admissible);
    mesher.options().setInt("Admissibility",0);
    // mesher.options().setInt("MaxLevel",3);
    mesher.getOptions();

    //! [beginRefLoop]
    gsField<> solField;
    for( int refLoop = 0; refLoop <= numRefinementLoops; refLoop++)
    {
        // bases = gsMultiBasis<>(patches);
        gsDebug<<"-------------------Loop = "<<refLoop<<"\n";
        //! [beginRefLoop]
        // gsWriteParaview(bases.basis(0),"basis0",200);
        // gsWriteParaview(bases.basis(1),"basis1",200);
        // gsWriteParaview(bases.basis(2),"basis2",200);
        // gsDebugVar(bases.basis(0).numElements());
        // gsDebugVar(bases.basis(1).numElements());
        // gsDebugVar(bases.basis(2).numElements());
        // --------------- solving ---------------

        //! [solverPart]
        // Construct assembler
        gsPoissonAssembler<real_t> PoissonAssembler(patches,bases,bcInfo,f);
        PoissonAssembler.options().setInt("DirichletValues", dirichlet::l2Projection);

        // Generate system matrix and load vector
        PoissonAssembler.assemble();

        // Initialize the conjugate gradient solver
        gsSparseSolver<>::CGDiagonal solver( PoissonAssembler.matrix() );

        // Solve the linear system
        gsMatrix<> solVector = solver.solve( PoissonAssembler.rhs() );

        // Construct the isogeometric solution
        gsMultiPatch<> sol;
        PoissonAssembler.constructSolution(solVector, sol);
        // Associate the solution to the patches (isogeometric field)
        gsField<> solField(patches, sol);

        gsInfo<<"NumDofs = "<<solVector.size()<<"\n";

        gsWriteParaview<>(solField, "adaptRef" + std::to_string(2*refLoop), 1000, true);
        for (size_t p = 0; p!=patches.nPatches(); ++p)
        {
            collection.addTimestep("adaptRef"+ std::to_string(2*refLoop) + std::to_string(p),2*refLoop,".vts");
            collection.addTimestep("adaptRef"+ std::to_string(2*refLoop) + std::to_string(p) + "_mesh",2*refLoop,".vtp");
        }

        //! [solverPart]

        // --------------- error estimation/computation ---------------

        //! [errorComputation]
        // Compute the error in the H1-seminorm ( = energy norm in this example )
        // using the known exact solution.
        gsExprEvaluator<> ev;
        gsDebugVar(PoissonAssembler.multiBasis().totalElements());
        ev.setIntegrationElements(PoissonAssembler.multiBasis());
        gsExprEvaluator<>::geometryMap Gm = ev.getMap(patches);
        auto is = ev.getVariable(sol);
        auto ms = ev.getVariable(g, Gm);

        // Get the element-wise norms.
        ev.integralElWise( ( igrad(is,Gm) - igrad(ms)).sqNorm()*meas(Gm) );
        std::vector<real_t> eltErrs  = ev.elementwise();
        //! [errorComputation]

        index_t offset = 0;
        for (index_t p = 0; p!=patches.nPatches(); p++)
        {
            auto first = eltErrs.begin() + offset;
            offset += bases.basis(p).numElements();
            auto last = eltErrs.begin() + offset;
            std::vector<real_t> tmpErrors(first,last);
            gsElementErrorPlotter<real_t> err_eh(bases.basis(p),tmpErrors);
            const gsField<> elemError_eh( patches.patch(p), err_eh, true );
            gsWriteParaview<>( elemError_eh, "error_elem_ref" + std::to_string(p), 1000, false);
        }

        // --------------- adaptive refinement ---------------
        for (index_t k = 0; k!=eltErrs.size(); k++)
            gsDebug<<eltErrs[k]<<"\n";

        gsHBoxContainer<2,real_t> refined, unrefined;
        if (coarsening)
        {
            mesher.markRef_into(eltErrs,refined);
            mesher.markCrs_into(eltErrs,refined,unrefined);
            mesher.refine(refined);
            mesher.unrefine(unrefined);
        }
        else
        {
            mesher.markRef_into(eltErrs,refined);
            gsDebugVar(refined);
            mesher.refine(refined);
        }

        for (size_t p = 0; p!=patches.nPatches(); ++p)
        {
            gsHBoxContainer<2,real_t> patchContainer = refined.patch(p);
            gsWriteParaview(patchContainer,"Refined_patch_" + std::to_string(p));
        }

        gsWriteParaview<>(solField, "adaptRef" + std::to_string(2*refLoop+1), 1000, true);
        for (size_t p = 0; p!=patches.nPatches(); ++p)
        {
            collection.addTimestep("adaptRef"+ std::to_string(2*refLoop+1) + std::to_string(p),2*refLoop+1,".vts");
            collection.addTimestep("adaptRef"+ std::to_string(2*refLoop+1) + std::to_string(p) + "_mesh",2*refLoop+1,".vtp");
        }

        //! [adaptRefinementPart]

        //! [repairInterfaces]
        // Call repair interfaces to make sure that the new meshes
        // match along patch interfaces.
        // patches.repairInterfaces();
        gsDebugVar(bases.totalElements());
        bases.repairInterfaces( patches.interfaces() );
        gsWriteParaview(bases.basis(0),"basis0",200);
        gsWriteParaview(bases.basis(1),"basis1",200);
        gsWriteParaview(bases.basis(2),"basis2",200);
        gsDebugVar(bases.totalElements());
        // bases = gsMultiBasis<>(patches);
        //! [repairInterfaces]

        // bases.repairInterfaces( patches.interfaces() );


        gsWriteParaview(patches,"mp",1000,true);

        //! [Export to Paraview]
        // Export the final solution
        // if( plot && refLoop == numRefinementLoops )
        // {
           // Write the computed solution to paraview files

        // }
        //! [Export to Paraview]
        mesher.rebuild();
    }

    collection.save();

    //! [Plot in Paraview]
    if( plot )
    {
       // Run paraview
       // gsFileManager::open("adaptRef.pvd");
    }
    //! [Plot in Paraview]
    else
    {
       gsInfo<<"Done. No output created, re-run with --plot to get a ParaView "
               "file containing Plotting image data.\n";
    }

    return EXIT_SUCCESS;

}// end main
