/** @file gsProjection.hpp

    @brief Implements projection methods

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst
*/

#include <gsPde/gsBoundaryConditions.h>
#include <gsAssembler/gsExprAssembler.h>
#include <gsAssembler/gsExprEvaluator.h>
#include <gsMatrix/gsSparseMatrix.h>
#include <gsMatrix/gsSparseSolver.h>
#include <gsCore/gsBoundary.h>

namespace gismo {

template<enum ProjectionNorm Norm, typename T>
void gsProjection<Norm,T>::_matrix(const gsMultiBasis<T>         & integrationBasis,
                                   const gsFunctionSet<T>        & projectionBasis,
                                   const gsFunctionSet<T>        & geometryMap,
                                         gsSparseMatrix<T>       & systemMatrix,
                                         short_t                   targetDim,
                                   const gsBoundaryConditions<T> & bc,
                                   const gsOptionList            & options,
                                   T alpha,
                                   T beta,
                                   T gamma)
{
    // Clear the result
    systemMatrix.clear();

    const short_t actualTargetDim = (targetDim < 0) ? geometryMap.targetDim() : targetDim;

    // Create an assembler
    gsExprAssembler<T> A(1,1);

    // Set the integration elements
    A.setIntegrationDomain(integrationBasis.domain());

    // Assign the space
    space u = A.getSpace(projectionBasis, actualTargetDim);

    // Assign the geometry map
    typename gsExprAssembler<T>::geometryMap G = A.getMap(geometryMap);

    // Set up the space
    u.setup(bc,dirichlet::interpolation,options.askInt("Continuity",-1));

    // Initialize the system
    A.initSystem();

    gsProjection<Norm,T>::template _assembleMatrix<Norm>(A,u,G,alpha,beta,gamma);
    const gsSparseMatrix<T> consistentMatrix = A.matrix();
    if (options.askSwitch("Lumped",false))
    {
        gsMatrix<T> ones = gsMatrix<T>::Ones(consistentMatrix.cols(), 1);
        gsMatrix<T> rowSums = consistentMatrix * ones;
        gsSparseEntries<T> entries;
        for (index_t i = 0; i < rowSums.rows(); ++i)
            entries.add(i, i, rowSums(i,0));
        systemMatrix.resize(consistentMatrix.rows(), consistentMatrix.cols());
        systemMatrix.setFrom(entries);
    }
    else
        A.matrix_into(systemMatrix);
}

template<enum ProjectionNorm Norm, typename T>
void gsProjection<Norm,T>::_rhs(const gsMultiBasis<T>         & integrationBasis,
                                const gsFunctionSet<T>        & projectionBasis,
                                const gsFunctionSet<T>        & geometryMap,
                                const gsFunctionSet<T>        & sourceFunction,
                                      gsMatrix<T>             & rhs,
                                const gsBoundaryConditions<T> & bc,
                                const gsOptionList            & options,
                                T alpha,
                                T beta,
                                T gamma)
{
    // Clear the result
    rhs.clear();

    // Create an assembler
    gsExprAssembler<T> A(1,1);

    // Set the integration elements
    A.setIntegrationDomain(integrationBasis.domain());

    // Assign the space
    space u = A.getSpace(projectionBasis,sourceFunction.targetDim());

    // Assign the source function
    auto f = A.getCoeff(sourceFunction);

    // Assign the geometry map
    typename gsExprAssembler<T>::geometryMap G = A.getMap(geometryMap);

    // Set up the space
    u.setup(bc,dirichlet::interpolation,options.askInt("Continuity",-1));

    // Initialize the system
    A.initSystem();

    gsProjection<Norm,T>::template _assembleRhs<Norm>(A,u,f,G,alpha,beta,gamma);
    A.rhs_into(rhs);
}

template<enum ProjectionNorm Norm, typename T>
void gsProjection<Norm,T>::_system(const gsMultiBasis<T>         & integrationBasis,
                                   const gsFunctionSet<T>        & projectionBasis,
                                   const gsFunctionSet<T>        & geometryMap,
                                   const gsFunctionSet<T>        & sourceFunction,
                                         gsSparseMatrix<T>       & systemMatrix,
                                         gsMatrix<T>             & rhs,
                                   const gsBoundaryConditions<T> & bc,
                                   const gsOptionList            & options,
                                   T alpha,
                                   T beta,
                                   T gamma)
{
    // Clear the results
    systemMatrix.clear();
    rhs.clear();

    // Create an assembler
    gsExprAssembler<T> A(1,1);

    // Set the integration elements
    A.setIntegrationDomain(integrationBasis.domain());

    // Assign the space
    space u = A.getSpace(projectionBasis,sourceFunction.targetDim());

    // Assign the source function
    auto f = A.getCoeff(sourceFunction);

    // Assign the geometry map
    typename gsExprAssembler<T>::geometryMap G = A.getMap(geometryMap);

    // Set up the space
    u.setup(bc,dirichlet::interpolation,options.askInt("Continuity",-1));

    // Initialize the system
    A.initSystem();

    gsProjection<Norm,T>::template _assembleMatrix<Norm>(A,u,G,alpha,beta,gamma);
    gsProjection<Norm,T>::template _assembleRhs<Norm>(A,u,f,G,alpha,beta,gamma);
    if (options.askSwitch("Lumped",false))
    {
        gsMatrix<T> ones = gsMatrix<T>::Ones(A.matrix().cols(), 1);
        gsMatrix<T> rowSums = A.matrix() * ones;
        gsSparseEntries<T> entries;
        for (index_t i = 0; i < rowSums.rows(); ++i)
            entries.add(i, i, rowSums(i,0));
        systemMatrix.resize(A.matrix().rows(), A.matrix().cols());
        systemMatrix.setFrom(entries);
        A.rhs_into(rhs);
    }
    else
    {
        A.matrix_into(systemMatrix);
        A.rhs_into(rhs);
    }
}

template<enum ProjectionNorm Norm, typename T>
T gsProjection<Norm,T>::_project(const gsMultiBasis<T>         & integrationBasis,
                                 const gsFunctionSet<T>        & projectionBasis,
                                 const gsFunctionSet<T>        & geometryMap,
                                 const gsFunctionSet<T>        & sourceFunction,
                                       gsMatrix<T>             & coefs,
                                 const gsBoundaryConditions<T> & bc,
                                 const gsOptionList            & options,
                                 T alpha,
                                 T beta,
                                 T gamma)
{
    // Clear the result
    coefs.clear();

    // Create an assembler
    gsExprAssembler<T> A(1,1);

    // Set the integration elements
    A.setIntegrationDomain(integrationBasis.domain());

    // Assign the space
    space u = A.getSpace(projectionBasis,sourceFunction.targetDim());

    // Assign the source function
    auto f = A.getCoeff(sourceFunction);

    // Assign the geometry map
    typename gsExprAssembler<T>::geometryMap G = A.getMap(geometryMap);

    // Set up the space
    u.setup(bc,dirichlet::interpolation,options.askInt("Continuity",-1));

    // Initialize the system
    A.initSystem();

    gsProjection<Norm,T>::template _assembleMatrix<Norm>(A,u,G,alpha,beta,gamma);
    gsProjection<Norm,T>::template _assembleRhs<Norm>(A,u,f,G,alpha,beta,gamma);
    if (options.askSwitch("Lumped",false))
    {
        gsMatrix<T> ones = gsMatrix<T>::Ones(A.matrix().cols(), 1);
        gsMatrix<T> diagVals = A.matrix() * ones;
        gsMatrix<T> rhsVals = A.rhs();
        coefs.resize(diagVals.rows(), 1);
        for (index_t i = 0; i < diagVals.rows(); ++i)
            coefs(i,0) = (diagVals(i,0) != T(0)) ? rhsVals(i,0) / diagVals(i,0) : T(0);
    }
    else
    {
        // Solve the system
        typename gsSparseSolver<T>::uPtr solver = gsSparseSolver<T>::get( options.askString("LinearSolver","SimplicialLDLT") );
        solver->compute(A.matrix());
        coefs = solver->solve(A.rhs());
    }

    if (options.askSwitch("ComputeError",true))
    {
        // Extract the solution and compute the error
        solution sol = A.getSolution(u, coefs);
        gsExprEvaluator<T> ev(A);
        return gsProjection<Norm,T>::template _computeError<Norm>(ev,sol,f,G,alpha,beta,gamma);
    }
    else
        return -1;
}

} // gismo
