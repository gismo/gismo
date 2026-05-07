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
                                   const gsOptionList            & options)
{
    // Clear the result
    systemMatrix.clear();

    // Create an assembler
    gsExprAssembler<T> A(1,1);

    // Set the integration elements
    A.setIntegrationDomain(integrationBasis.domain());

    // Assign the space
    space u = A.getSpace(projectionBasis,targetDim);

    // Assign the geometry map
    typename gsExprAssembler<T>::geometryMap G = A.getMap(geometryMap);

    // Set up the space
    u.setup(bc,dirichlet::interpolation,options.askInt("Continuity",-1));

    // Initialize the system
    A.initSystem();

    // assemble system
    if (!options.askSwitch("Lumped",false))
    {
        gsProjection<Norm,T>::template _assembleMatrix<Norm>(A,u,G);
        A.matrix_into(systemMatrix);
    }
    else
    {
        gsInfo<<"Warning: Lumped mass matrix is not implemented for the matrix-only assembly. Falling back to consistent mass matrix."<<std::endl;
        gsProjection<Norm,T>::template _assembleMatrix<Norm>(A,u,G);
        A.matrix_into(systemMatrix);
    }
}

template<enum ProjectionNorm Norm, typename T>
void gsProjection<Norm,T>::_rhs(const gsMultiBasis<T>         & integrationBasis,
                                const gsFunctionSet<T>        & projectionBasis,
                                const gsFunctionSet<T>        & geometryMap,
                                const gsFunctionSet<T>        & sourceFunction,
                                      gsMatrix<T>             & rhs,
                                const gsBoundaryConditions<T> & bc,
                                const gsOptionList            & options)
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

    gsProjection<Norm,T>::template _assembleRhs<Norm>(A,u,f,G);
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
                                   const gsOptionList            & options)
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

    // assemble system
    if (!options.askSwitch("Lumped",false))
    {
        gsProjection<Norm,T>::template _assembleMatrix<Norm>(A,u,G);
        gsProjection<Norm,T>::template _assembleRhs<Norm>(A,u,f,G);
        A.matrix_into(systemMatrix);
        A.rhs_into(rhs);
    }
    else
    {
        gsInfo<<"Warning: Lumped mass matrix is not implemented for the system assembly. Falling back to consistent mass matrix."<<std::endl;
    }
}

template<enum ProjectionNorm Norm, typename T>
T gsProjection<Norm,T>::_project(const gsMultiBasis<T>         & integrationBasis,
                                 const gsFunctionSet<T>        & projectionBasis,
                                 const gsFunctionSet<T>        & geometryMap,
                                 const gsFunctionSet<T>        & sourceFunction,
                                       gsMatrix<T>             & coefs,
                                 const gsBoundaryConditions<T> & bc,
                                 const gsOptionList            & options)
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

    // assemble system
    if (!options.askSwitch("Lumped",false))
    {
        gsProjection<Norm,T>::template _assembleMatrix<Norm>(A,u,G);
        gsProjection<Norm,T>::template _assembleRhs<Norm>(A,u,f,G);
        // Solve the system
        typename gsSparseSolver<T>::uPtr solver = gsSparseSolver<T>::get( options.askString("LinearSolver","SimplicialLDLT") );
        solver->compute(A.matrix());
        coefs = solver->solve(A.rhs());
    }
    else
    {
        gsProjection<Norm,T>::template _assembleMatrix<Norm>(A,u,G);
        gsProjection<Norm,T>::template _assembleRhs<Norm>(A,u,f,G);
        gsMatrix<> LHS = A.matrix() * gsMatrix<>::Ones(A.matrix().rows(),1);
        A.clearRhs();
        gsProjection<Norm,T>::template _assembleRhs<Norm>(A,u,f,G);
        gsMatrix<> RHS = A.rhs();
        coefs = LHS.cwiseInverse().cwiseProduct(RHS);
    }

    if (options.askSwitch("ComputeError",true))
    {
        // Extract the solution and compute the error
        solution sol = A.getSolution(u, coefs);
        gsExprEvaluator<T> ev(A);
        return gsProjection<Norm,T>::template _computeError<Norm>(ev,sol,f,G);
    }
    else
        return -1;
}

} // gismo
