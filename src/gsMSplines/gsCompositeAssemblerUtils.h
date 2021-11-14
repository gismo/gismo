/** @file gsCompositeAssemblerUtils.h

    @brief Provides utils for assemblers working with gsMappedBasis or Geometry.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): F. Buchegger
*/

#pragma once

#include <gsMSplines/gsMappedBasis.h>
#include <gsMSplines/gsCompositeHBasis.h>
#include <gsMSplines/gsCompositeGeom.h>
#include <gsMSplines/gsMappedSingleBasis.h>

#include <gsRecipeAssembler/gsRecipeAssemblerQMOpt.h>
#include <gsRecipeAssembler/gsRecipeAssemblerFitting.h>
#include <gsRecipeAssembler/gsRecipeAssemblerPoisson.h>
#include <gsRecipeAssembler/gsRecipeAssemblerBiharmonicSogn.h>

#include <gsAssembler/gsNormL2.h>
#include <gsAssembler/gsSeminormH1.h>
#include <gsAssembler/gsSeminormH2.h>

#include <gsCore/gsFieldCreator.h>
#include <gsIO/gsWriteParaview.h>
#include <gsUtils/gsStopwatch.h>
#include <gsSolver/gsSolverUtils.h>

#include <gsUtils/gsExportMatrix.h>

namespace gismo
{

//============================================== NORMS ===============================================//

real_t compositeL2Error(const gsField<real_t>& computedSol,
                        const gsFunction<real_t>& exactSol,
                        bool exactSolIsParam=false)
{
    gsNormL2<real_t> l2DistNorm(computedSol,exactSol,exactSolIsParam);
    //l2DistNorm.compute();
    gsNormL2<real_t> l2Norm(computedSol);
    //l2Norm.compute();
    return l2DistNorm.compute()/l2Norm.compute();
}

real_t compositeH1Error(const gsField<real_t>& computedSol,
                        const gsFunction<real_t>& exactSol,
                        bool exactSolIsParam=false)
{
    gsSeminormH1<real_t> h1DistNorm(computedSol,exactSol,exactSolIsParam);
    h1DistNorm.compute();
    gsSeminormH1<real_t> h1Norm(computedSol);
    h1Norm.compute();
    gsNormL2<real_t> l2DistNorm(computedSol,exactSol,exactSolIsParam);
    l2DistNorm.compute();
    gsNormL2<real_t> l2Norm(computedSol);
    l2Norm.compute();
    return sqrt(l2DistNorm.compute()*l2DistNorm.compute() + h1DistNorm.compute()*h1DistNorm.compute())/
            sqrt(l2Norm.compute()*l2Norm.compute() + h1Norm.compute()*h1Norm.compute());
}

real_t compositeH2Error(const gsField<real_t>& computedSol,
                        const gsFunction<real_t>& exactSol,
                        bool exactSolIsParam=false)
{
    gsSeminormH2<real_t> h2DistNorm(computedSol,exactSol,exactSolIsParam);
    //h2DistNorm.compute();
    gsSeminormH2<real_t> h2Norm(computedSol);
    //h2Norm.compute();
    gsSeminormH1<real_t> h1DistNorm(computedSol,exactSol,exactSolIsParam);
    //h1DistNorm.compute();
    gsSeminormH1<real_t> h1Norm(computedSol);
    //h1Norm.compute();
    gsNormL2<real_t> l2DistNorm(computedSol,exactSol,exactSolIsParam);
    //l2DistNorm.compute();
    gsNormL2<real_t> l2Norm(computedSol);
    //l2Norm.compute();
    return sqrt(l2DistNorm.compute()*l2DistNorm.compute() + h1DistNorm.compute()*h1DistNorm.compute() + h2DistNorm.compute()*h2DistNorm.compute())/
            sqrt(l2Norm.compute()*l2Norm.compute() + h1Norm.compute()*h1Norm.compute() + h2Norm.compute()*h2Norm.compute());
}

//============================================== UTILS ===============================================//

void addAllDirichletBoundaries(const gsMultiPatch<>& mp,gsFunction<real_t>& g,gsBoundaryConditions<real_t> & bcInfo)
{
    std::vector<patchSide> boundaries = mp.boundaries();
    for(size_t i = 0;i<boundaries.size();++i)
        bcInfo.addCondition(boundaries[i].patch, boundaries[i].side(),  condition_type::dirichlet, &g);
}

bool solveSystemCG(const gsSparseMatrix<>& matrix,const gsMatrix<> rhs,gsMatrix<>& solVector,real_t tolerance = 1e-10,size_t iterations=100000000)
{
    gsSparseSolver<>::CGDiagonal solver;
    solver.setMaxIterations(iterations);
    solver.setTolerance(tolerance);
    solver.compute(matrix);
    solVector = solver.solve( rhs );
    if(solver.error()>=tolerance)
        std::cout << "Solver error: " << solver.error() << "\n";
    return solver.error()<tolerance;
}

template<typename T>
bool solveSystemBiCGSTAB(const gsSparseMatrix<>& matrix,const gsMatrix<> rhs,gsMatrix<>& solVector,real_t tolerance = 1e-10,size_t iterations=100000000)
{
    gsSparseSolver<>::BiCGSTABILUT solver(matrix);
    solver.setMaxIterations(iterations);
    solver.setTolerance(tolerance);
    if (solver.preconditioner().info() != Eigen::Success)
    {
        gsWarn << "The preconditioner failed. Aborting...\n";
        return false;
    }
    solVector = solver.solve(rhs);
    if(solver.error()>=tolerance)
        std::cout << "Solver error: " << solver.error() << "\n";
    return solver.error()<tolerance;
}

template <typename T>
gsMatrix<T> reconstructFullCoefs( const gsMatrix<T>& dirichletVal,const gsMatrix<T>& solVal )
{
    const index_t sz  = dirichletVal.rows()+solVal.rows();
    const short_t dim = 1;

    gsMatrix<T> coeffs(sz, dim);
    for (index_t i = 0; i < sz; ++i)
    {
        if ( i<solVal.rows() ) // internal or interface
            coeffs.row(i) = solVal.row(i);
        else // eliminated Dof: fill with Dirichlet data
            coeffs.row(i).array() = dirichletVal.row( i-solVal.rows() );
    }
    return coeffs;
}

template <short_t d,class T>
void gsRefineMarkedElements(gsCompositeHBasis<d,T> & basis,
                            const std::vector<bool> & elMarked,
                            int refExtension = 0)
{
    const short_t dim = basis.dim();

    // numMarked: Number of marked cells on current patch, also currently marked cell
    // poffset  : offset index for the first element on a patch
    // globalCount: counter for the current global element index
    index_t numMarked, poffset = 0, globalCount = 0;

    // refBoxes: contains marked boxes on a given patch
    gsMatrix<T> refBoxes;

    for (size_t pn=0; pn < basis.nPatches(); ++pn )// for all patches
    {
        // Get number of elements to be refined on this patch
        const index_t numEl = basis.getBase(pn).numElements();
        numMarked = std::count_if(elMarked.begin() + poffset,
                                  elMarked.begin() + poffset + numEl,
                                  std::bind2nd(std::equal_to<bool>(), true) );
        poffset += numEl;
        refBoxes.resize(dim, 2*numMarked);
        //gsDebugVar(numMarked);
        numMarked = 0;// counting current patch element to be refined

        // for all elements in patch pn
        typename gsBasis<T>::domainIter domIt = basis.basis(pn).makeDomainIterator();
        for (; domIt->good(); domIt->next())
        {
            if( elMarked[ globalCount++ ] ) // refine this element ?
            {
                const gsVector<T> & ctr = domIt->centerPoint();

                // Construct degenerate box by setting both
                // corners equal to the center
                refBoxes.col(2*numMarked  ) =
                refBoxes.col(2*numMarked+1) = ctr;

                // Advance marked cells counter
                numMarked++;
            }
        }
        // Refine all of the found refBoxes in this patch
        basis.refineWithExtension( pn, refBoxes, refExtension,false );
    }
}

void getConvergenceRatios(const std::vector<real_t>& errors,std::vector<real_t>& rates)
{
    GISMO_ASSERT(errors.size()>1,"size too small for a rate");
    for(size_t i = 0;i<errors.size()-1;++i)
        rates.push_back(errors[i]/errors[i+1]);
}


//========================================= POISSON EQUATION =========================================//

void solvePoisson(gsMappedBasis<2,real_t> * compBasis,const gsBoundaryConditions<real_t> & bcInfo,
                  const gsMultiPatch<real_t> & domain,
                  const gsFunctionExpr<real_t>& f,
                  real_t& time,gsMatrix<real_t>& reconstructedSol)
{
    // recipe based version

    gsMatrix<>        rhs;
    gsMatrix<>        dirichlet_values;

    gsStopwatch stopwatch;

    gsPoissonPde<real_t>poissonPde(domain,bcInfo,f);
    gsRecipeAssemblerPoisson assembler(poissonPde);
    assembler.setDirichletStrategy(dirichlet::elimination);
    std::vector<gsBasis<real_t>*> bases;
    for (size_t i =0;i<compBasis->nPatches();++i)
        bases.push_back(&compBasis->getBase(i));
    std::vector<gsPhysicalSpace*> physSpaces;
    physSpaces.push_back(new gsPhysicalSpaceScalar(bases,domain,INVERSE_COMPOSITION,compBasis->getMapper()));
    assembler.setSpace(physSpaces);
    assembler.setZeroAverage(false);
    assembler.assemble();

    rhs=assembler.getSystemRhs();
    if ( assembler.getElimSize()>0 )
    {
        gsSparseSolver<>::CGDiagonal solver;
        dirichlet_values = solver.compute(assembler.getEliminatedMatrix()).solve(assembler.getEliminatedRhs());
        rhs-=  assembler.getRhsModMatrix()*dirichlet_values;
    }

    gsMatrix<> solVector;
    bool successfulSolve = solveSystemCG(assembler.getSystemMatrix(),rhs,solVector);
    if(!successfulSolve)
        std::cout<<"The solver was not able to reach the desired tolerance.\n";
    GISMO_ASSERT(successfulSolve,"The solver was not able to reach the desired tolerance.");

    reconstructedSol = assembler.reconstructSolution(0,solVector,dirichlet_values);

    time = stopwatch.stop();
}

void solvePoissonCollectResults_oneStep(
        gsMappedBasis<2,real_t> * compBasis,
        const gsBoundaryConditions<real_t> & bcInfo,
        const gsMultiPatch<real_t> & domain,
        const gsFunctionExpr<real_t>& g, const gsFunctionExpr<real_t>& f,
        int iteration,
        real_t& h,real_t& l2error,real_t& h1error,
        size_t & dof,real_t& time,int& plotNr,bool plot=false)
{
    gsMatrix<real_t> reconstructedSol;
    typename gsMappedBasis<2,real_t>::uPtr basis = compBasis->clone();
    solvePoisson(basis.release(),bcInfo,domain,f,time,reconstructedSol);
    std::cout << "goodSol: " << reconstructedSol << std::endl;
    std::vector<gsBasis<real_t> *> vc_refine_bases;
    for (size_t i = 0;i<compBasis->nPatches();i++)
    {
        vc_refine_bases.push_back(new gsMappedSingleBasis<2,real_t>(compBasis->getMappedSingleBasis(i)));
    }
    gsCompositeGeom<2,real_t> geo(*compBasis,reconstructedSol);
    gsMultiPatch<real_t> mp2 = geo.exportToPatches();
    gsField<real_t> sol(domain,mp2);

    if(plot)
    {
        gsFunctionExpr<real_t> g_copy = g;
        gsField<real_t> exact(domain, g_copy, false);
        gsWriteParaview( sol , "Poisson"+util::to_string(plotNr)+"solutionComposite"+util::to_string(iteration), compBasis->size()*2, true );
        gsWriteParaview( exact , "Poisson"+util::to_string(plotNr)+"solutionExact"+util::to_string(iteration), compBasis->size()*2, false );

        gsField<real_t> errorPtr = gsFieldCreator<real_t>::absError(sol,g);
        gsWriteParaview( errorPtr,"Poisson"+util::to_string(plotNr)+"solutionAbsErrComposite"+util::to_string(iteration), 10000, false);
        for(index_t i = 0;i<domain.size();i++)
        {
            gsMultiPatch<real_t> pt;
            pt.addPatch(domain.patch(i).clone());
            gsGradientField<real_t> grad(domain.patch(i),sol.function(i));
            gsField<real_t> gradField(pt,grad,true);
            gsWriteParaview(gradField , "Poisson"+util::to_string(plotNr)+"solutionGradFieldIter"+util::to_string(iteration)+"Patch"+util::to_string(i), 10000);
        }
        std::cout << "prints done - ";
    }
    // Find the element size
    h = math::pow( (real_t) vc_refine_bases[0]->size(), -1.0 / vc_refine_bases[0]->dim() );

    l2error = compositeL2Error(sol,g);
    h1error = compositeH1Error(sol,g);

    dof=compBasis->size();

    //cout<<"iteration " << iteration << " done!\n" << std::flush;
}

////======================================== BIHARMONIC EQUATION =======================================//

void solveBiharmonic(gsMappedBasis<2,real_t>*    compBasis,
               const gsBoundaryConditions<real_t> & /*bcInfo*/,
                     gsMultiPatch<real_t>*          domain,
                     gsFunctionExpr<real_t>       & f,
                     real_t                       & time,
                     real_t                       & conditionNr,
                     gsMatrix<real_t>             & reconstructedSol,
                     std::string                    desc,
                     int                            step)
{
    // recipe based version

    gsMatrix<>        dirichlet_values;

    gsStopwatch stopwatch;

    std::vector<gsBasis<real_t>*> bases;
    for (size_t i =0;i<compBasis->nPatches();++i)
        bases.push_back(&compBasis->getBase(i));

    gsConstantFunction<>  u_dir(0,domain->parDim());
    gsConstantFunction<>  u_neu(0,domain->parDim());
    gsPiecewiseFunction<real_t> f_picewise(f);
    gsBoundaryConditions<> bc1;
    gsBoundaryConditions<> bc2;

    std::vector<patchSide> boundaries = domain->boundaries();
    for(size_t i = 0;i<boundaries.size();++i)
    {
        bc1.addCondition(boundaries[i].patch, boundaries[i].side(),  condition_type::dirichlet, u_dir);
        bc2.addCondition(boundaries[i].patch, boundaries[i].side(),  condition_type::neumann,   u_neu);
    }


    gsBiharmonicPde<real_t> pde(*domain, bc1, bc2, f_picewise);
    gsRecipeAssemblerBiharmonicSogn assembler(pde);
    gsRecipeAssemblerBiharmonicSogn::Discretization spaces;
    spaces.push_back(new gsPhysicalSpaceScalar(bases,*domain,INVERSE_COMPOSITION,compBasis->getMapper()));
    assembler.setSpace(spaces);
    assembler.assemble();

    gsSparseMatrix<real_t> sysMat = assembler.getSystemMatrix();

    gsVector<real_t> diagVec = sysMat.diagonal();
    for(index_t i = 0;i<diagVec.rows();++i)
    {
        diagVec(i)=1/sqrt(diagVec(i));
    }
    sysMat = diagVec.asDiagonal()*sysMat*diagVec.asDiagonal();
    std::stringstream ss;
    ss << desc << "/mat" << step;
    exportMatrixToASCII(ss.str(),sysMat);
    conditionNr = 0;//gsSolverUtils<real_t>::conditionNumber(sysMat);


    gsMatrix<> solVector;
    bool successfulSolve = solveSystemCG(assembler.getSystemMatrix(),assembler.getSystemRhs(),solVector);
    if(!successfulSolve)
        std::cout<<"The solver was not able to reach the desired tolerance.\n";
    GISMO_ASSERT(successfulSolve,"The solver was not able to reach the desired tolerance.");

    dirichlet_values.setZero(assembler.getElimSize(),1);
    reconstructedSol=assembler.reconstructSolution(0,solVector,dirichlet_values);

    time = stopwatch.stop();
}

void solveBiharmonicCollectResults_oneStep(
              gsMappedBasis<2,real_t>*    compBasis,
        const gsBoundaryConditions<real_t> & bcInfo,
              gsMultiPatch<real_t>*          domain,
              gsFunctionExpr<real_t>       & g,
              gsFunctionExpr<real_t>       & /*dg*/,
              gsFunctionExpr<real_t>       & /*ddg*/,
              gsFunctionExpr<real_t>       & f,
              int                            iteration,
              real_t                       & h,
              real_t                       & l2error,
              real_t                       & h1error,
              real_t                       & h2error,
              size_t                       & dof,
              real_t                       & conditionNr,
              real_t                       & time,
              int                          & plotNr,
              std::string                    desc,
              bool                           plot=false)
{
    gsMatrix<real_t> reconstructedSol;
    typename gsMappedBasis<2,real_t>::uPtr basis = compBasis->clone();
    solveBiharmonic(basis.release(),bcInfo,domain,f,time,conditionNr,reconstructedSol,desc,iteration);
    std::vector<gsBasis<real_t> *> vc_refine_bases;
    for (size_t i = 0;i<compBasis->nPatches();i++)
    {
        vc_refine_bases.push_back(new gsMappedSingleBasis<2,real_t>(compBasis->getMappedSingleBasis(i)));
    }
    gsCompositeGeom<2,real_t> geo(*compBasis,reconstructedSol);
    gsMultiPatch<real_t> mp2 = geo.exportToPatches();
    gsField<real_t> sol(*domain,mp2);

    if(plot)
    {
        gsField<real_t> exact(*domain, g, false);
        gsWriteParaview( sol , "Biharmonic"+util::to_string(plotNr)+"solutionComposite"+util::to_string(iteration), 1000, true );
        gsWriteParaview( exact , "Biharmonic"+util::to_string(plotNr)+"solutionExact"+util::to_string(iteration), 1000, false );

        gsField<real_t> errorPtr = gsFieldCreator<real_t>::absError(sol,g);
        gsWriteParaview(errorPtr,"Biharmonic"+util::to_string(plotNr)+"solutionAbsErrComposite"+util::to_string(iteration), 1000//compGeom->getCompBasis().size()*2
                         , false);
        for(index_t i = 0;i<domain->size();i++)
        {
            gsMultiPatch<real_t> pt;
            pt.addPatch(domain->patch(i).clone());
            gsGradientField<real_t> grad(domain->patch(i),sol.function(i));
            gsField<real_t> gradField(pt,grad,true);
            gsWriteParaview(gradField , "Biharmonic"+util::to_string(plotNr)+"solutionGradFieldIter"+util::to_string(iteration)+"Patch"+util::to_string(i), 1000);
        }
        std::cout << "prints done - ";
    }
    // Find the element size
    h = math::pow( (real_t) vc_refine_bases[0]->size(), -1.0 / vc_refine_bases[0]->dim() );

    l2error = compositeL2Error(sol,g);
    h1error = compositeH1Error(sol,g);
    h2error = compositeH2Error(sol,g);

    dof=compBasis->size();

    //std::cout<<"iteration " << iteration << " done!\n" << std::flush;
}

////========================================== QUALITY MEASURES ========================================//

void solveMeasure(gsMappedBasis<2,real_t> * compBasis,
                  gsMultiPatch<real_t>& domain,
                  real_t& time,gsMatrix<real_t>& reconstructedSol,
                  gsMatrix<real_t>* qualityMeasureVals,
                  bool fixBoundaries,
                  const gsQualityMeasureWeights& weights)
{
    // recipe based version
    gsStopwatch stopwatch;

    std::vector<gsBasis<real_t>*> bases;
    for (size_t i =0;i<compBasis->nPatches();++i)
        bases.push_back(&compBasis->getBase(i));
    gsRecipeAssemblerQMOpt2D assembler(domain,bases,compBasis->getMapper(),weights,
                                           true,true,qualityMeasureVals==NULL,false,false,weights.m_approxPoints!=0);
    if(fixBoundaries)
        assembler.fixBoundaries();
    assembler.assemble();

    gsSparseMatrix<real_t>  stiff = assembler.getSystemMatrix();
    gsMatrix<real_t>        rhs = assembler.getSystemRhs();

    gsMatrix<> solVector;

    bool successfulSolve = solveSystemBiCGSTAB<real_t>(stiff,rhs,solVector,1e-10,10000);
    if(!successfulSolve)
        std::cout<<"The solver was not able to reach the desired tolerance.\n"<<std::endl;
    if(!successfulSolve)
        successfulSolve = solveSystemCG(stiff,rhs,solVector,1e-2,10000);
    GISMO_ASSERT(successfulSolve,"The solver was not able to reach the desired tolerance.");

    gsMatrix<>        dirichlet_values;
    dirichlet_values.setZero(assembler.getElimSize(),1);
    reconstructedSol = assembler.reconstructSolution(0,solVector,dirichlet_values);

    time = stopwatch.stop();
}

void calculateAreas(gsMappedBasis<2,real_t> * compBasis,
                    gsMultiPatch<real_t>& domain,
                    gsMatrix<real_t>& areas)
{
    // recipe based version
    std::vector<gsBasis<real_t>*> bases;
    for (size_t i =0;i<compBasis->nPatches();++i)
        bases.push_back(&compBasis->getBase(i));
    gsQualityMeasureWeights weights(compBasis->nPatches());
    gsRecipeAssemblerQMOpt2D assembler(domain,bases,compBasis->getMapper(),weights,
                                           false,false,false,true,false,false);
    assembler.assemble();
    assembler.getAreas(areas);
}

void calculateMeasures(gsMappedBasis<2,real_t> * compBasis,
                        gsMultiPatch<real_t>& domain,
                        gsMatrix<real_t>& qualityMeasures,
                       const gsQualityMeasureWeights& weights)
{
    // recipe based version
    std::vector<gsBasis<real_t>*> bases;
    for (size_t i =0;i<compBasis->nPatches();++i)
        bases.push_back(&compBasis->getBase(i));
    gsRecipeAssemblerQMOpt2D assembler(domain,bases,compBasis->getMapper(),weights,
                                           false,false,true,false,false,false);
    assembler.assemble();
    assembler.getQualityMeasureVals(qualityMeasures);
}

} // end namespace gismo
