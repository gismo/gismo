/** @file gsL2Projection.hpp

    @brief Performs L2 projections

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
#include <gsCore/gsBoxTopology.h>
#include <gsCore/gsBoundary.h>
#include <gsCore/gsAffineFunction.h>

namespace gismo {

template<typename T>
T gsL2Projection<T>::projectGeometry(    const gsMultiBasis<T> & basis,
                                            const gsMultiPatch<T> & geometry,
                                            gsMultiPatch<T> & result)
{
    result.clear();

    gsExprAssembler<T> A(1,1);
    gsMatrix<T> solVector;

    A.setIntegrationElements(basis);
    space u = A.getSpace(basis,geometry.targetDim());
    solution sol = A.getSolution(u, solVector);
    auto f = A.getCoeff(geometry);
    geometryMap G = A.getMap(geometry);

    u.setup(0);
    A.initSystem();

    // assemble system
    A.assemble(u*u.tr() * meas(G),u * f * meas(G));

    typename gsSparseSolver<T>::uPtr solver = gsSparseSolver<real_t>::get( "SimplicialLDLT" );
    solver->compute(A.matrix());
    solVector = solver->solve(A.rhs());

    gsExprEvaluator<> ev(A);
    gsDebugVar(ev.integral((sol-f).sqNorm()));

    sol.extract(result);
    result.computeTopology();
    result.closeGaps();

    return ev.integral((sol-f).sqNorm() * meas(G));
}

template<typename T>
T gsL2Projection<T>::projectGeometry(    const gsMultiBasis<T> & basis,
                                            const gsMultiPatch<T> & geometry,
                                            gsMatrix<T> & result)
{
    gsExprAssembler<T> A(1,1);

    A.setIntegrationElements(basis);
    space u = A.getSpace(basis,geometry.targetDim());
    auto f = A.getCoeff(geometry);
    geometryMap G = A.getMap(geometry);

    u.setup(-1);
    A.initSystem();

    // assemble system
    A.assemble(u*u.tr() * meas(G),u * f * meas(G));

    typename gsSparseSolver<T>::uPtr solver = gsSparseSolver<T>::get( "SimplicialLDLT" );
    solver->compute(A.matrix());
    result = solver->solve(A.rhs());

    solution sol = A.getSolution(u, result);
    gsExprEvaluator<> ev(A);
    gsDebugVar(ev.integral((sol-f).sqNorm()));
    return ev.integral((sol-f).sqNorm() * meas(G));
}

template<typename T>
T gsL2Projection<T>::projectGeometry(    const gsMultiBasis<T> & intbasis,
                                            const gsMappedBasis<2,T> & basis,
                                            const gsMultiPatch<T> & geometry,
                                            gsMatrix<T> & result)
{
    gsExprAssembler<T> A(1,1);

    A.setIntegrationElements(intbasis);
    space u = A.getSpace(basis,geometry.targetDim());
    auto f = A.getCoeff(geometry);
    geometryMap G = A.getMap(geometry);

    u.setup(-1);
    A.initSystem();

    // assemble system
    A.assemble(u*u.tr()*meas(G),u * f*meas(G));

    typename gsSparseSolver<T>::uPtr solver = gsSparseSolver<real_t>::get( "SimplicialLDLT" );
    solver->compute(A.matrix());
    result = solver->solve(A.rhs());

    solution sol = A.getSolution(u, result);
    gsExprEvaluator<> ev(A);
    gsDebugVar(ev.integral((sol-f).sqNorm()));
    ev.writeParaview((sol-f).sqNorm(),G,"error");
    return ev.integral((sol-f).sqNorm() * meas(G));
}

template<typename T>
T gsL2Projection<T>::projectFunction(    const gsMultiBasis<T> & basis,
                                            const gsFunctionSet<T> & source,
                                            const gsMultiPatch<T>   & geometry,
                                            gsMultiPatch<T> & result)
{
    result.clear();

    gsExprAssembler<T> A(1,1);
    gsMatrix<T> solVector;

    A.setIntegrationElements(basis);
    space u = A.getSpace(basis,source.targetDim());
    auto  f = A.getCoeff(source);
    solution sol = A.getSolution(u, solVector);
    geometryMap G = A.getMap(geometry);

    u.setup(-1);
    A.initSystem();

    // assemble system
    A.assemble(u*u.tr() * meas(G),u * f * meas(G));

    typename gsSparseSolver<T>::uPtr solver = gsSparseSolver<real_t>::get( "SimplicialLDLT" );
    solver->compute(A.matrix());
    solVector = solver->solve(A.rhs());

    sol.extract(result);
    result.computeTopology();
    result.closeGaps();
    gsExprEvaluator<> ev(A);
    return ev.integral((sol-f).sqNorm() * meas(G));
}

template<typename T>
T gsL2Projection<T>::projectFunction(    const gsMultiBasis<T> & basis,
                                            const gsFunctionSet<T> & source,
                                            const gsMultiPatch<T>   & geometry,
                                            gsMatrix<T> & result)
{
    gsExprAssembler<T> A(1,1);
    gsMatrix<T> solVector;

    A.setIntegrationElements(basis);
    space u = A.getSpace(basis,source.targetDim());
    auto  f = A.getCoeff(source);
    solution sol = A.getSolution(u, solVector);
    geometryMap G = A.getMap(geometry);

    u.setup(-1);
    A.initSystem();

    // assemble system
    A.assemble(u*u.tr() * meas(G),u * f * meas(G));

    typename gsSparseSolver<T>::uPtr solver = gsSparseSolver<real_t>::get( "SimplicialLDLT" );
    solver->compute(A.matrix());
    result = solver->solve(A.rhs());
    gsExprEvaluator<> ev(A);
    return ev.integral((sol-f).sqNorm() * meas(G));
}


template<typename T>
T gsL2Projection<T>::projectFunction(    const gsMultiBasis<T>   & intbasis,
                                            const gsMappedBasis<2,T>& basis,
                                            const gsFunctionSet<T>  & source,
                                            const gsMultiPatch<T>   & geometry,
                                            gsMatrix<T> & result)
{
    gsExprAssembler<T> A(1,1);

    A.setIntegrationElements(intbasis);
    space u = A.getSpace(basis,source.targetDim());
    auto  f = A.getCoeff(source);
    geometryMap G = A.getMap(geometry);

    u.setup(-1);
    A.initSystem();

    // assemble system
    A.assemble(u*u.tr()*meas(G),u * f *meas(G));

    typename gsSparseSolver<T>::uPtr solver = gsSparseSolver<real_t>::get( "SimplicialLDLT" );
    solver->compute(A.matrix());
    result = solver->solve(A.rhs());

    solution sol = A.getSolution(u, result);
    gsExprEvaluator<> ev(A);
    gsDebugVar(ev.integral((sol-f).sqNorm()));
    return ev.integral((sol-f).sqNorm() * meas(G));
}

template<typename T>
T gsL2Projection<T>::projectGeometryBoundaries(const gsMultiBasis<T> & basis,
                                            const gsMultiPatch<T> & geometry,
                                            gsMultiPatch<T> & result)
{
    result.clear();

    gsExprAssembler<T> A(1,1);
    gsMatrix<T> solVector;

    A.setIntegrationElements(basis);
    space u = A.getSpace(basis,geometry.geoDim());
    solution sol = A.getSolution(u, solVector);
    geometryMap G = A.getMap(geometry);

    std::vector<gsMultiPatch<T>> coords(geometry.geoDim());
    for (index_t p=0; p!=geometry.geoDim(); p++)
        coords[p] = geometry.coord(p);

    gsBoundaryConditions<T> bc;
    bc.setGeoMap(geometry);
    for (index_t p=0; p!=geometry.geoDim(); p++)
        for (typename gsMultiPatch<T>::const_biterator bit = geometry.bBegin(); bit != geometry.bEnd(); ++bit)
            bc.addCondition(*bit, condition_type::dirichlet,static_cast<gsFunctionSet<T>*>(&coords[p]), 0, true, p);

    u.setup(bc, dirichlet::l2Projection, 0);

    A.initSystem();

    // assemble system
    A.assemble(u*u.tr()*meas(G),u * G*meas(G));

    typename gsSparseSolver<T>::uPtr solver = gsSparseSolver<real_t>::get( "SimplicialLDLT" );
    solver->compute(A.matrix());
    solVector = solver->solve(A.rhs());

    sol.extract(result);
    result.computeTopology();
    result.closeGaps();
    gsExprEvaluator<> ev(A);
    return ev.integral((sol-G).sqNorm() * meas(G));
}


template<typename T>
T gsL2Projection<T>::projectGeometryPenalty(const gsMultiBasis<T> & basis,
                                            const gsMultiPatch<T> & geometry,
                                            gsMultiPatch<T> & result,
                                            T penalty)
{
    result.clear();

    gsExprAssembler<T> A(1,1);
    gsMatrix<T> solVector;

    A.setIntegrationElements(basis);
    space u = A.getSpace(basis,geometry.geoDim());
    solution sol = A.getSolution(u, solVector);
    geometryMap G = A.getMap(geometry);
    auto Gvar = A.getCoeff(geometry,G);
    element el = A.getElement();

    std::vector<gsMultiPatch<T>> coords(geometry.geoDim());
    for (index_t p=0; p!=geometry.geoDim(); p++)
        coords[p] = geometry.coord(p);

    gsBoundaryConditions<T> bc;
    bc.setGeoMap(geometry);
    for (index_t p=0; p!=geometry.geoDim(); p++)
        for (typename gsMultiPatch<T>::const_biterator bit = geometry.bBegin(); bit != geometry.bEnd(); ++bit)
            bc.addCondition(*bit, condition_type::dirichlet,static_cast<gsFunctionSet<T>*>(&coords[p]), 0, true, p);


    u.setup(bc, dirichlet::l2Projection, -1);


    std::vector<boundaryInterface> iFace;
    for (typename gsMultiPatch<T>::const_iiterator iit = geometry.iBegin(); iit != geometry.iEnd(); ++iit)
        iFace.push_back(*iit);

    A.initSystem();

    // assemble system
    A.assemble(u*u.tr()*meas(G),u * Gvar*meas(G));
    A.assembleBdr(
                  bc.get("Weak Dirichlet"),
                    -(penalty / el.area(G) * u * u.tr()) * tv(G).norm()
                );

    A.assembleIfc(iFace,
                    penalty / el.area(G) * u.left() * u.left().tr() * tv(G).norm()
                    ,
                    - penalty / el.area(G) * u.right()* u.left() .tr() * tv(G).norm()
                    ,
                    - penalty / el.area(G) * u.left() * u.right().tr() * tv(G).norm()
                    ,
                    penalty / el.area(G) * u.right()* u.right().tr() * tv(G).norm()
                     );

    typename gsSparseSolver<T>::uPtr solver = gsSparseSolver<real_t>::get( "SimplicialLDLT" );
    solver->compute(A.matrix());
    solVector = solver->solve(A.rhs());

    sol.extract(result);
    result.computeTopology();
    result.closeGaps();
    gsExprEvaluator<> ev(A);
    return ev.integral((Gvar-G).sqNorm() * meas(G));
}



} // gismo
