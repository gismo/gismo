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
#include <gsMatrix/gsSparseSolver.h>
#include <gsCore/gsBoxTopology.h>
#include <gsCore/gsBoundary.h>
#include <gsCore/gsAffineFunction.h>

namespace gismo {

template<typename T>
void gsL2Projection<T>::projectGeometry(const gsMultiBasis<T> & basis,
                                            const gsFunctionSet<T> & source,
                                            gsMultiPatch<T> & result)
{
    result.clear();

    gsExprAssembler<T> A(1,1);
    gsMatrix<T> solVector;

    A.setIntegrationElements(basis);
    space u = A.getSpace(basis,source.targetDim());
    solution sol = A.getSolution(u, solVector);
    auto G = A.getCoeff(source);

    u.setup(-1);
    A.initSystem();

    // assemble system
    A.assemble(u*u.tr(),u * G);

    typename gsSparseSolver<T>::uPtr solver = gsSparseSolver<real_t>::get( "SimplicialLDLT" );
    solver->compute(A.matrix());
    solVector = solver->solve(A.rhs());

    sol.extract(result);
    result.computeTopology();
    result.closeGaps();
}

template<typename T>
void gsL2Projection<T>::projectGeometryBoundaries(const gsMultiBasis<T> & basis,
                                            const gsMultiPatch<T> & source,
                                            gsMultiPatch<T> & result)
{
    result.clear();

    gsExprAssembler<T> A(1,1);
    gsMatrix<T> solVector;

    A.setIntegrationElements(basis);
    space u = A.getSpace(basis,source.geoDim());
    solution sol = A.getSolution(u, solVector);
    geometryMap G = A.getMap(source);

    std::vector<gsMultiPatch<T>> coords(source.geoDim());
    for (index_t p=0; p!=source.geoDim(); p++)
        coords[p] = source.coord(p);

    gsBoundaryConditions<T> bc;
    bc.setGeoMap(source);
    for (index_t p=0; p!=source.geoDim(); p++)
        for (typename gsMultiPatch<T>::const_biterator bit = source.bBegin(); bit != source.bEnd(); ++bit)
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
}


template<typename T>
void gsL2Projection<T>::projectGeometryPenalty(const gsMultiBasis<T> & basis,
                                            const gsMultiPatch<T> & source,
                                            gsMultiPatch<T> & result,
                                            T penalty)
{
    result.clear();

    gsExprAssembler<T> A(1,1);
    gsMatrix<T> solVector;

    A.setIntegrationElements(basis);
    space u = A.getSpace(basis,source.geoDim());
    solution sol = A.getSolution(u, solVector);
    geometryMap G = A.getMap(source);
    auto Gvar = A.getCoeff(source,G);
    element el = A.getElement();

    std::vector<gsMultiPatch<T>> coords(source.geoDim());
    for (index_t p=0; p!=source.geoDim(); p++)
        coords[p] = source.coord(p);

    gsBoundaryConditions<T> bc;
    bc.setGeoMap(source);
    for (index_t p=0; p!=source.geoDim(); p++)
    {
        for (typename gsMultiPatch<T>::const_biterator bit = source.bBegin(); bit != source.bEnd(); ++bit)
            bc.addCondition(*bit, condition_type::dirichlet,static_cast<gsFunctionSet<T>*>(&coords[p]), 0, true, p);
        // for (typename gsMultiPatch<T>::const_iiterator iit = source.iBegin(); iit != source.iEnd(); ++iit)
        // {
        //     bc.addCondition(iit->first(), condition_type::dirichlet,static_cast<gsFunctionSet<T>*>(&coords[p]), 0, true, p);
        //     bc.addCondition(iit->second(),condition_type::dirichlet,static_cast<gsFunctionSet<T>*>(&coords[p]), 0, true, p);
        // }

    }


    u.setup(bc, dirichlet::l2Projection, -1);


    std::vector<boundaryInterface> iFace;
    for (typename gsMultiPatch<T>::const_iiterator iit = source.iBegin(); iit != source.iEnd(); ++iit)
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
}



} // gismo
