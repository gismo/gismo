/** @file embeddedCurve_example.cpp

    @brief Example for embedded curve

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): C. Chianese, H.M. Verhelst
*/

#include <iostream>

#include <gismo.h>
#include <gsAssembler/gsEmbeddingUtils.h>

using namespace gismo;

int main(int argc, char *argv[])
{

    gsCmdLine cmd("TODO");
    // cmd.addSwitch("plot", "Plot result in ParaView format", plot);
    // cmd.addSwitch("trim", "Basic trim/merge operations", trim);
    // cmd.addString("file", "Path to the file to be read", filename);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    // Make a BSpline surface
    gsKnotVector<> kv_u(0,1,2,3);
    gsKnotVector<> kv_v(0,1,2,3);
    gsTensorBSplineBasis<2, real_t> basis_s(kv_u, kv_v);
    gsMatrix<> coefs_s (basis_s.size(), 3);
    coefs_s << 0, 0, 0,
               0.166, 0, 0.130,
               0.5, 0, 0.215,
               0.833, 0, 0.130,
               1, 0, 0,
               0, 0.166, 0.130,
               0.166, 0.166, 0.260,
               0.5, 0.166, 0.347,
               0.833, 0.166, 0.260,
               1, 0.166, 0.130,
               0, 0.5, 0.215,
               0.166, 0.5, 0.347,
               0.5, 0.5, 0.433,
               0.833, 0.5, 0.347,
               1, 0.5, 0.215,
               0, 0.833, 0.130,
               0.166, 0.833, 0.260,
               0.5, 0.833, 0.347,
               0.833, 0.833, 0.260,
               1, 0.833, 0.130,
               0, 1, 0,
               0.166, 1, 0.130,
               0.5, 1, 0.215,
               0.833, 1, 0.130,
               1, 1, 0;

    gsTensorBSpline<2, real_t>  surf(basis_s, coefs_s);

    // Make a BSpline curve within the surface source domain
    gsKnotVector<> kv_c(0, 1, 2, 3); //start,end,interior knots, start/end multiplicity
    gsMatrix<> coefs_c(5, 2); //u,v;..
    coefs_c << 0, 0,
             0.2618, 0.053,
             0.5, 0.5,
             0.738, 0.9465,
             1,1;

    coefs_c.col(0) = coefs_c.col(0) * kv_u.last();
    coefs_c.col(1) = coefs_c.col(1) * kv_v.last();

    gsBSpline<> curve( kv_c, give(coefs_c));

    gsWriteParaview(surf,"surface",1000,true);
    gsWriteParaview(curve,"curve");

    // 2D->2D
    gsMultiPatch<> mp_surf;
    mp_surf.addPatch(surf);

    gsInfo<<"The surface (R^"<<mp_surf.patch(0).domainDim()<<" -> R^"<<mp_surf.patch(0).targetDim()<<") is:\n"<<mp_surf.patch(0)<<"\n";
    // 1D->2D
    gsMultiPatch<> mp_curve;
    mp_curve.addPatch(curve);
    gsInfo<<"The curve (R^"<<mp_curve.patch(0).domainDim()<<" -> R^"<<mp_curve.patch(0).targetDim()<<") is:\n"<<mp_curve.patch(0)<<"\n";

    gsMultiBasis<> basis_curve(mp_curve);
    gsMultiBasis<> basis_surf(mp_surf);

    // Declare the expression assembler
    gsExprAssembler<> exprAssembler(1,1);

    // Set integration elements
    // Register expressions inside assembler
    auto G_curve = exprAssembler.getMap(mp_curve);
    auto G_surf  = exprAssembler.getMap(mp_surf);
    auto u_curve = exprAssembler.getSpace(basis_curve); // needed for validation
    auto u_surf  = exprAssembler.getSpace(basis_surf,mp_surf.geoDim());

    // Expression to assemble
    auto expr    = u_surf * u_surf.tr();

    gsBoundaryConditions<> bc;
    bc.addCondition(boundary::north, condition_type::dirichlet, 0, 0, false, 0 );
    bc.setGeoMap(mp_surf);

    // Computes Dirichlet values and eliminates Dirichlet DoFs
    u_surf.setup(bc, dirichlet::interpolation, 0);

    // exprAssembler.initSystem();
    index_t N = exprAssembler.numDofs();
    gsSparseMatrix<> matrix(N,N);

    // Setup domain iterator
    const gsBasis<> & basis = basis_curve.basis(0);
    gsBasis<>::domainIter domIt = basis.makeDomainIterator();

    // Define quadrature rule
    gsOptionList options;
    options.addInt("quRule", "", 1);
    options.addReal("quA", "", 1.0);
    options.addInt("quB", "", 1);
    typename gsQuadRule<>::uPtr QuRule = gsQuadrature::getPtr(basis,options);

    // Loop over the curve elements
    gsMatrix<> quPointsCurve, quPointsSurface, quPointsPhysical;
    gsVector<> quWeights;

    // Make expression data
    const gsExprHelper<real_t>::Ptr exprData = exprAssembler.exprData();
    // Parse the expression
    exprData->parse(expr);
    for (; domIt->good(); domIt->next() )
    {
        // Map the Quadrature rule to the element
        QuRule->mapTo( domIt->lowerCorner(), domIt->upperCorner(), quPointsCurve, quWeights);

        // MAYBE NOT NEEDED
        // Map the quadrature points to the surface
        mp_curve.patch(0).eval_into(quPointsCurve, quPointsSurface);
        mp_surf.patch(0).eval_into(quPointsSurface, quPointsPhysical);
        // gsWriteParaviewPoints(quPointsPhysical,"quPointsPhysical");

        gsDebug<<"Curve integration points: "<<quPointsCurve<<"\n";
        gsDebug<<"Surface integration points: "<<quPointsSurface<<"\n";


        // Loop over quadrature points
        for (index_t k = 0; k != quWeights.rows(); ++k)
        {
            // Precompute the expression
            // exprData->points() = quPointsCurve;
            exprData->points() = quPointsSurface.col(k);
            // exprData->weights() = quWeights;
            exprData->precompute(0); // updates the actives for every element

            // Create a local element matrix
            gsMatrix<> localMat = quWeights[k] * expr.eval(0);

            // Push the local element matrix inside the big system
            const expr::gsFeSpace<real_t> & v = expr.rowVar();
            const expr::gsFeSpace<real_t> & u = expr.colVar();
            const index_t rd                  = v.dim();//row
            const index_t cd                  = u.dim();//col
            const gsDofMapper  & rowMap       = v.mapper();
            const gsDofMapper  & colMap       = u.mapper();
            gsMatrix<index_t> & rowInd0       = const_cast<gsMatrix<index_t>&>(v.data().actives);
            gsMatrix<index_t> & colInd0       = const_cast<gsMatrix<index_t>&>(u.data().actives);

            // Push
            for (index_t r = 0; r != rd; ++r)
            {
                const index_t rls = r * rowInd0.rows();     //local stride
                for (index_t i = 0; i != rowInd0.rows(); ++i)
                {
                    const index_t ii = rowMap.index(rowInd0.at(i),v.data().patchId,r); //N_i
                    if ( rowMap.is_free_index(ii) )
                    {
                        for (index_t c = 0; c != cd; ++c)
                        {
                            const index_t cls = c * colInd0.rows();     //local stride
                            for (index_t j = 0; j != colInd0.rows(); ++j)
                            {
                                if ( 0 == localMat(rls+i,cls+j) ) continue;

                                const index_t jj = colMap.index(colInd0.at(j),u.data().patchId,c); // N_j
                                if ( colMap.is_free_index(jj) )
                                {
                                    matrix.coeffRef(ii, jj) += localMat(rls+i,cls+j);
                                }
                            }
                        }
                    }
                }
            }

        }



        // // Create a local element matrix
        // gsMatrix<> localMat;
        // // Integrate the local element matrix
        // localMat = quWeights[0] * expr.eval(0);
        // for (index_t k = 1; k != quWeights.rows(); ++k)
        //     localMat += quWeights[k] * expr.eval(k);

        // gsDebugVar(localMat);

    }
    gsInfo<<"Assembly done\n";


    gsDebugVar(matrix.toDense());

    return EXIT_SUCCESS;
}

template <class T>
T curveCoordinate(const gsGeometry<T> & geometry, const T & value, const short_t & dir)
{
    GISMO_ASSERT(geometry.parDim() == 1, "The geometry must be a curve with parameter dimension 1");
    GISMO_ASSERT(dir < geometry.targetDim(), "The direction must be less than the target dimension of the geometry");
    gsMatrix<T> u(1,1);
    u(0,0) = value;
    return geometry.eval(u)(dir,0);
    // return geometry.deriv(u)(dir,0);
};