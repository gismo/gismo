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

    index_t type = 2;

    // Number of refinement loops to be done
    int numRefinementLoops = 4;

    gsCmdLine cmd("Tutorial on solving a Poisson problem.");
    cmd.addSwitch("plot", "Create a ParaView visualization file with the solution", plot);
    cmd.addSwitch("coarsen", "Also coarsen", coarsening);
    cmd.addSwitch("admissible", "Admissible refinement", admissible);
    cmd.addInt("r", "numLoop", "number of refinement loops", numRefinementLoops);
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


    // patches.degreeElevate(2);

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

    gsAdaptiveMeshing<real_t> patchMesher(patches);
    patchMesher.options().setInt("CoarsenRule",2);
    patchMesher.options().setInt("RefineRule",2);
    patchMesher.options().setReal("CoarsenParam",0.9);
    patchMesher.options().setReal("RefineParam",0.9);
    patchMesher.options().setInt("CoarsenExtension",CExt);
    patchMesher.options().setInt("RefineExtension",Ext);
    patchMesher.options().setSwitch("Admissible",admissible);
    patchMesher.options().setInt("Admissibility",0);
    // patchMesher.options().setInt("MaxLevel",3);
    patchMesher.getOptions();

    gsAdaptiveMeshing<real_t> mesher(bases);
    mesher.options().setInt("CoarsenRule",2);
    mesher.options().setInt("RefineRule",2);
    mesher.options().setReal("CoarsenParam",0.9);
    mesher.options().setReal("RefineParam",0.9);
    mesher.options().setInt("CoarsenExtension",CExt);
    mesher.options().setInt("RefineExtension",Ext);
    mesher.options().setSwitch("Admissible",admissible);
    mesher.options().setInt("Admissibility",0);
    mesher.options().setInt("Verbose",1);
    // mesher.options().setInt("MaxLevel",3);
    mesher.getOptions();

    //! [beginRefLoop]
    gsField<> solField;
    for( int refLoop = 0; refLoop <= numRefinementLoops; refLoop++)
    {
        // bases = gsMultiBasis<>(patches);
        gsDebug<<"-------------------Loop = "<<refLoop<<"\n";
        //! [beginRefLoop]

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

        //! [adaptRefinementPart]
        // Mark elements for refinement, based on the computed local errors and
        // the refinement-criterion and -parameter.
        std::vector<bool> elMarked( eltErrs.size() );

        gsRefineMarkedElements( bases, elMarked, 1 );

        std::vector<bool> elCMarked( eltErrs.size() );
        for (index_t k=0; k!=eltErrs.size(); k++)
            eltErrs[k] = -eltErrs[k];

        gsMarkElementsForRef( eltErrs, adaptRefCrit, adaptRefParam, elCMarked);

        for (index_t k=0; k!=elMarked.size(); k++)
            gsInfo<<eltErrs[k]<<"\t"<<elCMarked[k]<<"\n";

        // Get the element-wise norms.
        ev.integralElWise( ( igrad(is,Gm) - igrad(ms)).sqNorm()*meas(Gm) );
        std::vector<real_t> eltCErrs  = ev.elementwise();

        //! [errorComputation]

        for (index_t k=0; k!=eltCErrs.size(); k++)
            gsInfo<<eltCErrs[k]<<"\n";

        offset = 0;
        for (index_t p = 0; p!=patches.nPatches(); p++)
        {
            auto first = eltCErrs.begin() + offset;
            offset += bases.basis(p).numElements();
            auto last = eltCErrs.begin() + offset;
            std::vector<real_t> tmpErrors(first,last);
            gsElementErrorPlotter<real_t> err_eh(bases.basis(p),tmpErrors);
            const gsField<> elemError_eh( patches.patch(p), err_eh, true );
            gsWriteParaview<>( elemError_eh, "error_elem_crs" + std::to_string(p), 1000, false);
        }

        for (index_t k = 0; k!=eltErrs.size(); k++)
            gsDebug<<eltErrs[k]<<"\n";

        mesher.mark(eltCErrs);
        mesher.refine();

        gsMarkElementsForRef( eltCErrs, adaptRefCrit, adaptRefParam, elCMarked);

        for (index_t k=0; k!=elMarked.size(); k++)
            gsInfo<<eltCErrs[k]<<"\t"<<elCMarked[k]<<"\n";

        // _toContainer
        std::vector<gsHBoxContainer<2,real_t>> container(patches.nPieces());
        index_t c = 0;
        gsBasis<> * basis = nullptr;

        gsMultiPatch<> * mp;
        gsMultiBasis<> * mb;
        typename gsBasis<>::domainIter domIt;
        gsHDomainIterator<real_t,2> * domHIt = nullptr;
        for (index_t patchInd=0; patchInd < patches.nPieces(); ++patchInd)
        {
            // Initialize domain element iterator
            if ( (mp = dynamic_cast<gsMultiPatch<>*>(&patches)) ) basis = &(mp->basis(patchInd));
            if ( (mb = dynamic_cast<gsMultiBasis<>*>(&patches)) ) basis = &(mb->basis(patchInd));
            GISMO_ASSERT(basis!=nullptr,"Object is not gsMultiBasis or gsMultiPatch");
            // for all elements in patch pn
            domIt  = basis->makeDomainIterator();
            domHIt = dynamic_cast<gsHDomainIterator<real_t,2> *>(domIt.get());
            GISMO_ENSURE(domHIt!=nullptr,"Domain should be 2 dimensional for flattening");

            for (; domHIt->good(); domHIt->next() )
            {
                if (elCMarked[c])
                    container.at(patchInd).add(gsHBox<2,real_t>(domHIt,patchInd));
                c++;
            }
        }
        // ! _toContainer

        return 0;


        if (coarsening)
        {
            mesher.mark(eltErrs);
            mesher.adapt();
            // patchMesher.mark(eltErrs);
            // patchMesher.adapt();
        }
        else
        {
            mesher.mark(eltErrs);
            mesher.refine();
            // patchMesher.mark(eltErrs);
            // patchMesher.refine();
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
        bases.repairInterfaces( patches.interfaces() );
        // bases = gsMultiBasis<>(patches);
        // patches.repairInterfaces();
        //! [repairInterfaces]

        // bases.repairInterfaces( patches.interfaces() );


        gsWriteParaview(patches,"mp",1000,true);
        gsWriteParaview(bases.basis(0),"basis",1000,true);

        //! [Export to Paraview]
        // Export the final solution
        // if( plot && refLoop == numRefinementLoops )
        // {
           // Write the computed solution to paraview files

        // }
        //! [Export to Paraview]

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
