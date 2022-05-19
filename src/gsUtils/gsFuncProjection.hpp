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

namespace gismo {

template<typename T>
void gsFuncProjection<T>::L2(   const gsMultiBasis<T> & basis,
                                const gsMultiPatch<T> & source,
                                gsMultiPatch<T> & result,
                                bool fixSides,
                                bool fixInteriorVerts)
{
    std::vector<gsMultiPatch<T>> coords(source.geoDim());
    for (index_t p=0; p!=source.geoDim(); p++)
        coords[p] = source.coord(p);

    gsBoundaryConditions<T> bc;
    bc.setGeoMap(source);
    if (fixSides)
    {
        for (index_t p=0; p!=source.geoDim(); p++)
            for (typename gsMultiPatch<T>::const_biterator bit = source.bBegin(); bit != source.bEnd(); ++bit)
                bc.addCondition(*bit, condition_type::dirichlet,static_cast<gsFunctionSet<T>*>(&coords[p]), 0, false, p);

    }
    if (fixInteriorVerts)
    {
        gsVector<bool> pars;
        gsMatrix<T> values;
        std::vector<std::vector<patchCorner> > Vs, result;
        source.getEVs(Vs);
        result.insert(result.end(), Vs.begin(), Vs.end());
        source.getOVs(Vs);
        result.insert(result.end(), Vs.begin(), Vs.end());
        for (index_t p=0; p!=source.geoDim(); p++)
            for (std::vector<std::vector<patchCorner>>::const_iterator cont=result.begin(); cont!=result.end(); cont++ )
                for (std::vector<patchCorner>::const_iterator it=cont->begin(); it!=cont->end(); it++ )
                {
                    gsDebug<<"Patch "<<it->patch<<"; Corner "<<it->corner()<<"\n";
                    it->corner().parameters_into(source.domainDim(),pars);
                    source.patch(it->patch).eval_into(pars.template cast<T>(),values);
                    bc.addCornerValue(it->corner(), values(p,0), it->patch, 0, p); // corner, value (which is the function value at the corner), patch index, unknown, component
                }
    }

    _L2(basis,source,bc,result);
}

template<typename T>
void gsFuncProjection<T>::L2(   const gsMultiBasis<T> & basis,
                                const gsFunctionSet<T> & source,
                                gsMultiPatch<T> & result)
{
    // std::vector<gsMultiPatch<T>> coords(source.geoDim());
    // for (index_t p=0; p!=source.geoDim(); p++)
    //     coords[p] = source.coord(p);

    gsBoundaryConditions<T> bc;
    // // bc.setGeoMap(source);
    //
    // TO DO:
    // Take topology from multibasis
    // Find a way to apply the component of a function set as boundary condition
    //
    // if (fixSides)
    // {
    //     for (index_t p=0; p!=source.geoDim(); p++)
    //         for (typename gsFunctionSet<T>::const_biterator bit = source.bBegin(); bit != source.bEnd(); ++bit)
    //             bc.addCondition(*bit, condition_type::dirichlet,static_cast<gsFunctionSet<T>*>(&coords[p]), 0, false, p);

    // }
    // if (fixInteriorVerts)
    // {

    //     for (index_t p=0; p!=source.geoDim(); p++)
    //         for (typename gsFunctionSet<T>::const_biterator bit = source.bBegin(); bit != source.bEnd(); ++bit)
    //             bc.addCondition(*bit, condition_type::dirichlet,static_cast<gsFunctionSet<T>*>(&coords[p]), 0, false, p);
    // }

    _L2(basis,source,bc,result);
}

template<typename T>
void gsFuncProjection<T>::_L2(  const gsMultiBasis<T> & basis,
                                const gsFunctionSet<T> & source,
                                const gsBoundaryConditions<T> & bc,
                                gsMultiPatch<T> & result)
{
    result.clear();

    gsExprAssembler<T> A(1,1);
    gsMatrix<T> solVector;

    A.setIntegrationElements(basis);
    space u = A.getSpace(basis,source.targetDim());
    solution sol = A.getSolution(u, solVector);
    auto G = A.getCoeff(source);

    u.setup(bc, dirichlet::l2Projection, 0);

    A.initSystem();
    A.assemble(u*u.tr(),u * G); // Remove the measure

    typename gsSparseSolver<T>::uPtr solver = gsSparseSolver<real_t>::get( "SimplicialLDLT" );
    solver->compute(A.matrix());
    solVector = solver->solve(A.rhs());

    sol.extract(result);
    result.computeTopology();
    result.closeGaps();
}

} // gismo
