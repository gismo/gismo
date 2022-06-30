/** @file gsDirichletValues.h

    @brief The functions compute Dirichlet degrees of freedom using various methods.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris, H.M. Verhelst
*/

#include <gsUtils/gsPointGrid.h>

namespace gismo {

namespace expr
{
template<class T> class gsFeSpace;
};

template<class T>
void gsDirichletValues(
    const gsBoundaryConditions<T> & bc,
    const index_t dir_values,
    const expr::gsFeSpace<T> & u)
{
    if ( bc.container("Dirichlet").empty() && bc.cornerValues().empty()) return;

    const gsDofMapper & mapper = u.mapper();
    gsMatrix<T> & fixedDofs = const_cast<expr::gsFeSpace<T>&>(u).fixedPart();
    fixedDofs.setZero(u.mapper().boundarySize(), 1 );

    switch ( dir_values )
    {
    case dirichlet::homogeneous :
    case dirichlet::user :
        // If we have a homogeneous problem then fill with zeros
        break;
    case dirichlet::interpolation:
        gsDirichletValuesByTPInterpolation(u,bc);
        break;
    case dirichlet::l2Projection:
        gsDirichletValuesByL2Projection(u,bc);
        break;
    default:
        GISMO_ERROR("Something went wrong with Dirichlet values: "<< dir_values);
    }

     // Corner values -- todo
    for ( typename gsBoundaryConditions<T>::const_citerator it = bc.cornerBegin(); it != bc.cornerEnd(); ++it )
    {
        if(it->unknown != u.id())
            continue;

        const int k = it->patch;
        const gsBasis<T> & basis = u.source().basis(k);
        const int i  = basis.functionAtCorner(it->corner);
        const index_t com = it->component;

        for (index_t r = 0; r!=u.dim(); ++r)
        {
            if (com!=-1 && r!=com) continue;
            const int ii = mapper.bindex( i , k, r );
            fixedDofs.at(ii) = it->value;
        }
    }
}

template<class T>
void gsDirichletValuesByTPInterpolation(const expr::gsFeSpace<T> & u,
                                        const gsBoundaryConditions<T> & bc)
{
    const index_t parDim = u.source().domainDim();

    std::vector< gsVector<T> > rr;
    gsMatrix<index_t> boundary;
    gsVector<T> b(1);
    gsMatrix<T> fpts, tmp;

    gsMatrix<T> & fixedDofs = const_cast<expr::gsFeSpace<T>&>(u).fixedPart();
    fixedDofs.setZero(u.mapper().boundarySize(), 1 );

    // Iterate over all patch-sides with Boundary conditions
    typedef gsBoundaryConditions<T> bcList;
    for ( typename bcList::const_iterator it =  bc.begin("Dirichlet");
          it != bc.end("Dirichlet") ; ++it )
    {
        if( it->unknown()!=u.id() ) continue;

        const index_t com = it->unkComponent();

        const int k = it->patch();
        const gsBasis<T> & basis = u.source().basis(k);

        // Get dofs on this boundary
        boundary = basis.boundary(it->side());

        // Get the side information
        const int dir = it->side().direction( );
        const index_t param = (it->side().parameter() ? 1 : 0);

        // Get basis on the boundary
        typename gsBasis<T>::uPtr h = basis.boundaryBasis(it->side());

        //
        for (index_t r = 0; r!=u.dim(); ++r)
        {
            if (com!=-1 && r!=com) continue;

            // If the condition is homogeneous then fill with zeros
            if ( it->isHomogeneous() )
            {
                for (index_t i=0; i!= boundary.size(); ++i)
                {
                    const int ii = u.mapper().bindex( boundary.at(i) , k, r );
                    fixedDofs.at(ii) = 0;
                }
                continue;
            }

            // Compute grid of points on the face ("face anchors")
            rr.clear();
            rr.reserve( parDim );

            for ( int i=0; i < parDim; ++i)
            {
                if ( i==dir )
                {
                    b[0] = ( basis.component(i).support() ) (0, param);
                    rr.push_back(b);
                }
                else
                {
                    rr.push_back( basis.component(i).anchors().transpose() );
                }
            }

            // GISMO_ASSERT(it->function()->targetDim() == u.dim(),
            //              "Given Dirichlet boundary function does not match problem dimension."
            //              <<it->function()->targetDim()<<" != "<<u.dim()<<"\n");

            // Compute dirichlet values
            if ( it->parametric() )
                fpts = it->function()->piece(it->patch()).eval( gsPointGrid<T>( rr ) );
            else
            {
                const gsFunctionSet<T> & gmap = bc.geoMap();
                fpts = it->function()->piece(it->patch()).eval(  gmap.piece(it->patch()).eval(  gsPointGrid<T>( rr ) )  );
            }

            // Interpolate dirichlet boundary
            typename gsGeometry<T>::uPtr geo = h->interpolateAtAnchors(fpts);
            const gsMatrix<T> & dVals = geo->coefs();

            // Save corresponding boundary dofs
            for (index_t l=0; l!= boundary.size(); ++l)
            {
                const int ii = u.mapper().bindex( boundary.at(l) , k, r );
                fixedDofs.at(ii) = dVals.at(l);
            }
        }
    }
}

// Not called and used, todo
template<class T> void
gsDirichletValuesInterpolationTP(const expr::gsFeSpace<T> & u,
                                 const boundary_condition<T> & bc,
                                 gsMatrix<index_t> & boundary,
                                 gsMatrix<T> & values)
{
    const index_t parDim = u.source().domainDim();

    const gsFunctionSet<T> & gmap = bc.geoMap();
    std::vector< gsVector<T> > rr;

    gsVector<T> b(1);
    gsMatrix<T> fpts, tmp;

    gsMatrix<T> & fixedDofs = const_cast<expr::gsFeSpace<T>&>(u).fixedPart();

    if( bc.unknown()!=u.id() ) { boundary.clear(); values.clear(); return; }

    const int k = bc.patch();
    const gsBasis<T> & basis = u.source().basis(k);

    // Get dofs on this boundary
    boundary = basis.boundary(bc.side());

    // If the condition is homogeneous then fill with zeros
    if ( bc.isHomogeneous() )
    {
        const index_t com = bc.unkComponent();
        values.setZero(boundary.size(), (-1==com ? u.dim():1) );
        return;
    }

    // Get the side information
    int dir = bc.side().direction( );
    index_t param = (bc.side().parameter() ? 1 : 0);

    // Compute grid of points on the face ("face anchors")
    rr.clear();
    rr.reserve( parDim );

    for ( int i=0; i < parDim; ++i)
    {
        if ( i==dir )
        {
            b[0] = ( basis.component(i).support() ) (0, param);
            rr.push_back(b);
        }
        else
        {
            rr.push_back( basis.component(i).anchors().transpose() );
        }
    }

    // GISMO_ASSERT(bc.function()->targetDim() == u.dim(),
    //              "Given Dirichlet boundary function does not match problem dimension."
    //              <<bc.function()->targetDim()<<" != "<<u.dim()<<"\n");

    // Compute dirichlet values
    if ( bc.parametric() )
        fpts = bc.function()->eval( gsPointGrid<T>( rr ) );
    else
    {
        const gsFunctionSet<T> & gmap = bc.geoMap();
        fpts = bc.function()->eval(  gmap.piece(bc.patch()).eval(  gsPointGrid<T>( rr ) )  );
    }

    /*
      if ( fpts.rows() != u.dim() )
      {
      // assume scalar
      tmp.resize(u.dim(), fpts.cols());
      tmp.setZero();
      gsDebugVar(!dir);
      tmp.row(!dir) = (param ? 1 : -1) * fpts; // normal !
      fpts.swap(tmp);
      }
    */

    // Interpolate dirichlet boundary
    typename gsBasis<T>::uPtr h = basis.boundaryBasis(bc.side());
    typename gsGeometry<T>::uPtr geo = h->interpolateAtAnchors(fpts);
    values = give( geo->coefs() );
}


template<class T>
void gsDirichletValuesByL2Projection( const expr::gsFeSpace<T> & u,
                                      const gsBoundaryConditions<T> & bc)
{
    const gsFunctionSet<T> & gmap = bc.geoMap();

    const gsDofMapper & mapper = u.mapper();
    gsMatrix<T> & fixedDofs = const_cast<expr::gsFeSpace<T>& >(u).fixedPart();

    // Set up matrix, right-hand-side and solution vector/matrix for
    // the L2-projection
    gsSparseEntries<T> projMatEntries;
    gsMatrix<T>        globProjRhs;
    globProjRhs.setZero(u.mapper().boundarySize(), 1 );

    // Temporaries
    gsVector<T> quWeights;
    gsMatrix<T> basisVals, rhsVals;
    gsMatrix<index_t> globIdxAct, globBasisAct;

    gsMapData<T> md(NEED_MEASURE | SAME_ELEMENT);

    // eltBdryFcts stores the row in basisVals/globIdxAct, i.e.,
    // something like a "element-wise index"
    std::vector<index_t> eltBdryFcts;

    //const gsMultiPatch<T> & mp = static_cast<const gsMultiPatch<T> &>(gmap);

    // Iterate over all patch-sides with Dirichlet-boundary conditions
    typedef gsBoundaryConditions<T> bcList;
    for (typename bcList::const_iterator iter = bc.begin("Dirichlet");
         iter != bc.end("Dirichlet"); ++iter)
    {
        const int unk = iter->unknown();
        if(unk != u.id()) continue;

        const index_t com = iter->unkComponent();// == -1 ? 0 : iter->unkComponent(); // TODO should loop

        const int patchIdx   = iter->patch();
        const gsBasis<T> & basis = u.source().basis(patchIdx);
        const gsFunction<T> & patch = gmap.function(patchIdx);

        // Set up quadrature to degree+1 Gauss points per direction,
        // all lying on iter->side() except from the direction which
        // is NOT along the element
        gsGaussRule<T> bdQuRule(basis, 1.0, 1, iter->side().direction());

        // Create the iterator along the given part boundary.
        typename gsBasis<T>::domainIter bdryIter = basis.makeDomainIterator(iter->side());


        for (; bdryIter->good(); bdryIter->next())
        {
            bdQuRule.mapTo(bdryIter->lowerCorner(), bdryIter->upperCorner(),
                           md.points, quWeights);

            patch.computeMap(md);

            // Indices involved here:
            // --- Local index:
            // Index of the basis function/DOF on the patch.
            // Does not take into account any boundary or interface conditions.
            // --- Global Index:
            // Each DOF has a unique global index that runs over all patches.
            // This global index includes a re-ordering such that all eliminated
            // DOFs come at the end.
            // The global index also takes care of glued interface, i.e., corresponding
            // DOFs on different patches will have the same global index, if they are
            // glued together.
            // --- Boundary Index (actually, it's a "Dirichlet Boundary Index"):
            // The eliminated DOFs, which come last in the global indexing,
            // have their own numbering starting from zero.

            // Get the global indices (second line) of the local
            // active basis (first line) functions/DOFs:
            basis.active_into(md.points.col(0), globBasisAct);

            // Compute the basis function values because they don't
            // depend on the component
            basis.eval_into(md.points, basisVals);

            for (index_t r = 0; r!=u.dim(); ++r)
            {
                if (com!=-1 && r!=com) continue;

                mapper.localToGlobal(globBasisAct, patchIdx, globIdxAct,r);


                // Out of the active functions/DOFs on this element, collect all those
                // which correspond to a boundary DOF.
                // This is checked by calling mapper.is_boundary_index( global Index )

                // eltBdryFcts stores the row in basisVals/globIdxAct, i.e.,
                // something like a "element-wise index"
                eltBdryFcts.clear();
                eltBdryFcts.reserve(mapper.boundarySize());
                for (index_t i = 0; i < globIdxAct.rows(); i++)
                {
                    if (mapper.is_boundary_index(globIdxAct.at(i)))
                    {
                        eltBdryFcts.push_back(i);
                    }
                }

                // the values of the boundary condition are stored
                // to rhsVals. Here, "rhs" refers to the right-hand-side
                // of the L2-projection, not of the PDE.

                // If the condition is homogeneous then fill with zeros
                if ( iter->isHomogeneous() )
                {
                    rhsVals.setZero(1,md.points.size());
                }
                else
                {
                    if ( iter->parametric() )
                        rhsVals = iter->function()->piece(patchIdx).eval(md.points);
                    else
                        rhsVals = iter->function()->piece(patchIdx).eval(gmap.piece(patchIdx).eval(md.points));
                }

                // Do the actual assembly:
                for (index_t k = 0; k < md.points.cols(); k++)
                {
                    const T weight_k = quWeights[k] * md.measure(k);

                    // Only run through the active boundary functions on the element:
                    for (size_t i0 = 0; i0 < eltBdryFcts.size(); i0++)
                    {
                        // Each active boundary function/DOF in eltBdryFcts has...
                        // ...the above-mentioned "element-wise index"
                        const index_t i = eltBdryFcts[i0];
                        // ...the boundary index.
                        const index_t ii = mapper.global_to_bindex(globIdxAct.at(i));

                        for (size_t j0 = 0; j0 < eltBdryFcts.size(); j0++)
                        {
                            const index_t j = eltBdryFcts[j0];
                            const index_t jj = mapper.global_to_bindex(globIdxAct.at(j));

                            // Use the "element-wise index" to get the needed
                            // function value.
                            // Use the boundary index to put the value in the proper
                            // place in the global projection matrix.
                            projMatEntries.add(ii, jj, weight_k * basisVals(i, k) * basisVals(j, k));
                        } // for j

                        globProjRhs.at(ii) += weight_k * basisVals(i, k) * rhsVals.at(k);

                    } // for i
                } // for k
            }// for r
        } // bdryIter
    } // boundaryConditions-Iterator

    gsSparseMatrix<T> globProjMat(mapper.boundarySize(), mapper.boundarySize());
    globProjMat.setFrom(projMatEntries);
    globProjMat.makeCompressed();

    // Solve the linear system:
    // The position in the solution vector already corresponds to the
    // numbering by the boundary index. Hence, we can simply take them
    // for the values of the eliminated Dirichlet DOFs.
    typename gsSparseSolver<T>::CGDiagonal solver;
    fixedDofs = solver.compute(globProjMat).solve(globProjRhs);
} // computeDirichletDofsL2Proj



}; // namespace gismo
