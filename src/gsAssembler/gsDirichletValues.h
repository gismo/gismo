/** @file gsDirichletValues.h

    @brief The functions compute Dirichlet degrees of freedom using various methods.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris, H.M. Verhelst
*/

#pragma once

#include <gsUtils/gsPointGrid.h>
#include <gsCore/gsDofMapper.h>
#include <gsAssembler/gsAssemblerOptions.h>
#include <gsPde/gsBoundaryConditions.h>
#include <gsTensor/gsTensorBasis.h>
#include <gsHSplines/gsHTensorBasis.h>
#include <gsHSplines/gsTHBSplineBasis.h>
#include <gsCore/gsGeometrySlice.h>
#include <gsCore/gsComposedFunction.h>
#include <gsUtils/gsQuasiInterpolate.h>

namespace gismo {

namespace expr
{
template<class T> class gsFeSpace;
};

namespace internal {

/// \brief Whether dirichlet::quasiInterpolation applies to \a basis.
///
/// - Tensor-product B-spline and NURBS bases (parametric dimension 1 to 4),
///   through basis.source().
/// - Truncated hierarchical B-splines (THB) and rational THB (dimension 1 to 4),
///   through basis.source(). With Q^l a local interpolant reproducing the
///   level-l tensor space and supported in one level-l element of
///   \f$\Omega^l\setminus\Omega^{l+1}\f$,
///   \f$Q(f)=\sum_l\sum_{i\in I_l}\lambda_{i,l}(f)\,T_{i,l}\f$ reproduces the
///   THB space (preservation of coefficients),
///   see H. Speleers, C. Manni, Numer. Math. 132 (2016) 155-184.
/// - Non-truncated hierarchical B-splines (HB, dimension 1 to 4), tested on the
///   basis itself: the level-by-level residual quasi-interpolation in
///   gsQuasiInterpolate::localIntpl reproduces the HB space.
///
/// Mapped bases, rational HB and any other basis are excluded. For
/// d >= 2 the boundary basis is quasi-interpolated, and every listed type has a
/// boundaryBasis of a listed type. For d = 1 the side is a point (see
/// dirichletSideQuasiInterpolation), so localIntpl is never called; the same
/// end-point rule serves dirichlet::interpolation and dirichlet::automatic on
/// 1-D patches.
template<class T>
bool dirichletQuasiInterpolationSupported(const gsBasis<T> & basis)
{
    const gsBasis<T> & src = basis.source();
    return dynamic_cast<const gsTensorBasis<1,T>*>(&src)             ||
           dynamic_cast<const gsTensorBasis<2,T>*>(&src)             ||
           dynamic_cast<const gsTensorBasis<3,T>*>(&src)             ||
           dynamic_cast<const gsTensorBasis<4,T>*>(&src)             ||
           dynamic_cast<const gsTHBSplineBasis<1,T,true>*>(&src)     ||
           dynamic_cast<const gsTHBSplineBasis<2,T,true>*>(&src)     ||
           dynamic_cast<const gsTHBSplineBasis<3,T,true>*>(&src)     ||
           dynamic_cast<const gsTHBSplineBasis<4,T,true>*>(&src)     ||
           dynamic_cast<const gsTHBSplineBasis<1,T,false>*>(&basis)  ||   // HB: the basis itself, so
           dynamic_cast<const gsTHBSplineBasis<2,T,false>*>(&basis)  ||   // rational HB is excluded
           dynamic_cast<const gsTHBSplineBasis<3,T,false>*>(&basis)  ||
           dynamic_cast<const gsTHBSplineBasis<4,T,false>*>(&basis);
}

/// \brief Whether dirichlet::interpolation applies to \a basis.
///
/// For d >= 2 it is anchor interpolation, the tensor-product branch: tensor-product
/// B-spline bases interpolate per direction at the Greville points, which are
/// unisolvent by the Schoenberg-Whitney conditions. Rational tensor bases qualify
/// through basis.source(), whose indices and anchors they share.
/// Hierarchical bases are excluded for d >= 2: collocation at the anchors of a
/// mixed-level boundary basis can be singular, since active functions of
/// different levels can share a Greville point.
/// Hierarchical bases use dirichlet::quasiInterpolation instead, see
/// dirichletQuasiInterpolationSupported.
/// On a 1-D patch the side is an end point and every basis accepted by
/// dirichletQuasiInterpolationSupported is accepted (end-point rule, see
/// dirichletSideQuasiInterpolation).
template<class T>
bool dirichletInterpolationSupported(const gsBasis<T> & basis)
{
    if ( 1 == basis.domainDim() )
        return dirichletQuasiInterpolationSupported(basis);
    const gsBasis<T> & src = basis.source();
    return dynamic_cast<const gsTensorBasis<2,T>*>(&src)   ||
           dynamic_cast<const gsTensorBasis<3,T>*>(&src)   ||
           dynamic_cast<const gsTensorBasis<4,T>*>(&src);
}

/// \brief Coefficients of the quasi-interpolant of Dirichlet data on side
/// \a side of the basis \a basis (tensor, NURBS, THB, rational THB or HB).
///
/// The data on the side is \a fun (parametric condition, \a geo == nullptr) or
/// \a fun composed with the geometry map \a geo (condition in physical
/// coordinates). For d >= 2 it is quasi-interpolated in
/// basis.boundaryBasis(side) by gsQuasiInterpolate::localIntpl, which reproduces
/// data in the span of that boundary basis exactly, but is not interpolatory at
/// points. For d = 1 the side is an end point and the coefficient is the data
/// value there: with an open knot vector exactly one function is non-zero at
/// the end point and it equals 1. This d = 1 branch is the end-point rule shared
/// by dirichlet::interpolation, dirichlet::quasiInterpolation and
/// dirichlet::automatic.
///
/// \param[out] coefs basis.boundary(side).size() x fun.targetDim(); row l is
/// the coefficient of the patch function basis.boundary(side)(l).
///
/// Cost O(n (p+1)^{3(d-1)}) for n side functions of degree p (one dense
/// (p+1)^{d-1} LU per function), OpenMP-parallel over the functions.
template<class T>
void dirichletSideQuasiInterpolation(const gsBasis<T> & basis,
                                     const boxSide side,
                                     const gsFunction<T> & fun,
                                     const gsFunction<T> * geo,
                                     gsMatrix<T> & coefs)
{
    // localIntpl runs an omp parallel for, where a throw aborts: all checks first.
    GISMO_ENSURE(dirichletQuasiInterpolationSupported(basis),
                 "Dirichlet quasi-interpolation is not implemented for this basis (supported: tensor-product B-spline, NURBS, THB, rational THB, HB). Use `dirichlet::l2Projection` instead.");
    const short_t d = basis.domainDim();
    if ( nullptr == geo )
    {
        GISMO_ENSURE(fun.domainDim() == d, "Parametric Dirichlet function has domain dimension "
                     << fun.domainDim() << ", but the patch has parametric dimension " << d << ".");
    }
    else
    {
        GISMO_ENSURE(geo->domainDim() == d, "Geometry map has domain dimension "
                     << geo->domainDim() << ", but the patch has parametric dimension " << d << ".");
        GISMO_ENSURE(fun.domainDim() == geo->targetDim(), "Dirichlet function has domain dimension "
                     << fun.domainDim() << ", but the geometry map has target dimension "
                     << geo->targetDim() << ".");
    }

    if ( 1 == d )
    {
        // The side of a 1-D patch is an end point. With an open (clamped) knot vector exactly one
        // function, basis.boundary(side)(0), is non-zero there and it equals 1 (also for NURBS and
        // hierarchical bases: the others vanish, partition of unity); the dof is the data value.
        const gsMatrix<index_t> bnd = basis.boundary(side);
        GISMO_ENSURE(1 == bnd.size(), "Boundary of side " << side << " of a 1-D patch has "
                     << bnd.size() << " functions, but exactly one is expected.");
        gsMatrix<T> pt(1,1);
        pt(0,0) = basis.support()(0, side.parameter() ? 1 : 0);
        GISMO_ENSURE(math::abs(basis.evalSingle(bnd(0,0), pt)(0,0) - (T)1) < (T)1e-10,
                     "Dirichlet values on a 1-D patch need an open (clamped) knot vector: the end "
                     "function is not 1 at the end point.");
        coefs = (nullptr == geo ? fun.eval(pt) : fun.eval(geo->eval(pt))).transpose();   // 1 x targetDim
        return;
    }

    const index_t dir = side.direction();
    const T par = basis.support()(dir, side.parameter() ? 1 : 0);
    typename gsBasis<T>::uPtr h = basis.boundaryBasis(side);
    GISMO_ENSURE(basis.boundary(side).size() == h->size(),
                 "Boundary basis of side " << side << " has " << h->size()
                 << " functions, but the patch has " << basis.boundary(side).size() << " on it.");

    // Non-owning: onSide and data only point to fun, geo and each other, all alive here.
    const gsGeometrySlice<T> onSide(nullptr == geo ? &fun : geo, dir, par);
    if ( nullptr == geo )
        gsQuasiInterpolate<T>::localIntpl(*h, onSide, coefs);
    else
    {
        const gsComposedFunction<T> data(&onSide, &fun);
        gsQuasiInterpolate<T>::localIntpl(*h, data, coefs);
    }
}

/// \brief Anchor interpolation of the Dirichlet condition \a cond of \a u on a
/// tensor-product (or rational tensor) patch; writes the entries of
/// \a fixedDofs belonging to the dofs on cond.side().
/// Parametric dimension d >= 2 only; 1-D patches use dirichletSideInterpolationDofs.
template<class T>
void dirichletSideAnchorInterpolation(const expr::gsFeSpace<T> & u,
                                      const gsBoundaryConditions<T> & bc,
                                      const boundary_condition<T> & cond,
                                      gsMatrix<T> & fixedDofs)
{
    const index_t parDim = u.source().domainDim();
    gsMatrix<T> fpts, pts;
    const index_t com = cond.unkComponent();
    const int k = cond.patch();
    const gsBasis<T> & basis = u.source().basis(k);

    // Get dofs on this boundary
    const gsMatrix<index_t> boundary = basis.boundary(cond.side());

    // Get the side information
    const int dir = cond.side().direction( );
    const index_t param = (cond.side().parameter() ? 1 : 0);

    // Get basis on the boundary
    typename gsBasis<T>::uPtr h = basis.boundaryBasis(cond.side());

    for (index_t r = 0; r!=u.dim(); ++r)
    {
        if (com!=-1 && r!=com) continue;

        // If the condition is homogeneous then fill with zeros
        if ( cond.isHomogeneous() )
        {
            for (index_t i=0; i!= boundary.size(); ++i)
            {
                const int ii = u.mapper().bindex( boundary.at(i) , k, r );
                fixedDofs.at(ii) = 0;
            }
            continue;
        }

        // Compute grid of points on the face ("face anchors")
        const gsMatrix<T> banchors = h->anchors();
        pts.resize(parDim, banchors.cols());
        for (index_t i = 0, j = 0; i != parDim; ++i)
        {
            if ( i==dir )
                pts.row(i).setConstant( basis.support()(dir, param) );
            else
                pts.row(i) = banchors.row(j++);
        }

        // Compute dirichlet values
        if ( cond.parametric() )
            fpts = cond.function()->piece(cond.patch()).eval( pts );
        else
        {
            GISMO_ENSURE(bc.hasGeoMap(), "gsDirichletValues: a Dirichlet condition given in "
                         "physical coordinates needs the geometry map, but the boundary "
                         "conditions carry none; call bc.setGeoMap(...).");
            const gsFunctionSet<T> & gmap = bc.geoMap();
            fpts = cond.function()->piece(cond.patch()).eval(  gmap.piece(cond.patch()).eval(  pts )  );
        }

        // Interpolate dirichlet boundary
        typename gsGeometry<T>::uPtr geo = h->interpolateAtAnchors(fpts);
        const gsMatrix<T> & dVals = geo->coefs();

        // Save corresponding boundary dofs
        const index_t cc = (-1==com ? r : 0);
        GISMO_ENSURE( cc < dVals.cols(),
                      "Dirichlet function has target dimension "<< dVals.cols()
                      <<", which cannot supply component "<< cc <<".");
        for (index_t l=0; l!= boundary.size(); ++l)
        {
            const int ii = u.mapper().bindex( boundary.at(l) , k, r );
            fixedDofs.at(ii) = dVals(l, cc);
        }
    }
}

/// \brief Quasi-interpolation of the Dirichlet condition \a cond of \a u on a
/// patch whose basis passed dirichletQuasiInterpolationSupported; writes the
/// entries of \a fixedDofs belonging to the dofs on cond.side().
/// The caller checks the support of the basis.
template<class T>
void dirichletSideQuasiInterpolationDofs(const expr::gsFeSpace<T> & u,
                                         const gsBoundaryConditions<T> & bc,
                                         const boundary_condition<T> & cond,
                                         gsMatrix<T> & fixedDofs)
{
    const index_t k = cond.patch();
    const gsBasis<T> & basis = u.source().basis(k);
    const index_t com = cond.unkComponent();
    const gsMatrix<index_t> boundary = basis.boundary(cond.side());

    if ( cond.isHomogeneous() )
    {
        for (index_t r = 0; r!=u.dim(); ++r)
        {
            if (com!=-1 && r!=com) continue;
            for (index_t l=0; l!= boundary.size(); ++l)
                fixedDofs.at( u.mapper().bindex(boundary.at(l), k, r) ) = 0;
        }
        return;
    }

    const gsFunction<T> * fun =
        dynamic_cast<const gsFunction<T>*>( &cond.function()->piece(k) );
    GISMO_ENSURE(nullptr != fun, "Dirichlet data on patch " << k << " is not a gsFunction.");

    const index_t ccMax = (-1==com ? u.dim()-1 : 0);
    GISMO_ENSURE( ccMax < fun->targetDim(),
                  "Dirichlet function has target dimension "<< fun->targetDim()
                  <<", which cannot supply component "<< ccMax <<".");

    const gsFunction<T> * geo = nullptr;
    if ( !cond.parametric() )
    {
        GISMO_ENSURE(bc.hasGeoMap(), "gsDirichletValues: a Dirichlet condition given in "
                     "physical coordinates needs the geometry map, but the boundary "
                     "conditions carry none; call bc.setGeoMap(...).");
        geo = dynamic_cast<const gsFunction<T>*>( &bc.geoMap().piece(k) );
        GISMO_ENSURE(nullptr != geo, "The geometry map of patch " << k << " is not a gsFunction.");
    }

    gsMatrix<T> coefs;
    dirichletSideQuasiInterpolation(basis, cond.side(), *fun, geo, coefs);

    for (index_t r = 0; r!=u.dim(); ++r)
    {
        if (com!=-1 && r!=com) continue;
        const index_t cc = (-1==com ? r : 0);
        for (index_t l=0; l!= boundary.size(); ++l)
            fixedDofs.at( u.mapper().bindex(boundary.at(l), k, r) ) = coefs(l, cc);
    }
}

/// \brief Dirichlet dofs of the condition \a cond of \a u for dirichlet::interpolation:
/// anchor interpolation (dirichletSideAnchorInterpolation) on a patch of parametric dimension
/// d >= 2; on a 1-D patch the side is an end point and the dof is the data value there
/// (the end-point rule of dirichletSideQuasiInterpolation, open knot vector required).
/// The caller checks the support of the basis (dirichletInterpolationSupported).
template<class T>
void dirichletSideInterpolationDofs(const expr::gsFeSpace<T> & u,
                                    const gsBoundaryConditions<T> & bc,
                                    const boundary_condition<T> & cond,
                                    gsMatrix<T> & fixedDofs)
{
    if ( 1 == u.source().basis(cond.patch()).domainDim() )
        dirichletSideQuasiInterpolationDofs(u, bc, cond, fixedDofs);
    else
        dirichletSideAnchorInterpolation(u, bc, cond, fixedDofs);
}

} // namespace internal

/// \brief Computes the Dirichlet dofs of \a u (its fixedPart()) by the method
/// \a dir_values (a dirichlet::values) and applies the corner values of \a bc.
/// \param sameElement asserts that each boundary quadrature batch lies in a single Bezier element of
/// the geometry map; passing false evaluates the map per point. Read from no option list.
/// Applies to \c dirichlet::l2Projection, and to \c dirichlet::automatic when that resolves to l2
/// projection; interpolation and quasi-interpolation ignore it.
template<class T>
void gsDirichletValues(
    const gsBoundaryConditions<T> & bc,
    const index_t dir_values,
    const expr::gsFeSpace<T> & u,
    const bool sameElement = true)
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
        gsDirichletValuesByTPInterpolation(u, bc);
        break;
    case dirichlet::quasiInterpolation:
        gsDirichletValuesByQuasiInterpolation(u, bc);
        break;
    case dirichlet::automatic:
        gsDirichletValuesAutomatic(u, bc, sameElement);
        break;
    case dirichlet::l2Projection:
        gsDirichletValuesByL2Projection(u, bc, sameElement);
        break;
    default:
        GISMO_ERROR("Something went wrong with Dirichlet values: "<< dir_values);
    }

     // Corner values -- todo
    for ( typename gsBoundaryConditions<T>::const_citerator it = bc.cornerBegin(); it != bc.cornerEnd(); ++it )
    {
        if(it->unknown!=-1 && it->unknown != u.id())
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

/// \brief Computes the Dirichlet dofs of \a u by interpolation at the anchors of the
/// boundary bases (dirichlet::interpolation). For d >= 2 tensor-product (incl. rational) patches
/// only: throws on any other basis there. On a 1-D patch (any basis supported by
/// quasi-interpolation) the end dof is the data value at the end point; the knot vector must be open.
/// The support of all patches is checked before any dof is written.
template<class T>
void gsDirichletValuesByTPInterpolation(const expr::gsFeSpace<T> & u,
                                        const gsBoundaryConditions<T> & bc)
{
    gsMatrix<T> & fixedDofs = const_cast<expr::gsFeSpace<T>&>(u).fixedPart();
    fixedDofs.setZero(u.mapper().boundarySize(), 1 );

    // Iterate over all patch-sides with Boundary conditions
    typedef gsBoundaryConditions<T> bcList;
    for ( typename bcList::const_iterator it =  bc.begin("Dirichlet");
          it != bc.end("Dirichlet") ; ++it )
    {
        if( it->unknown()!=u.id() ) continue;

        // For d >= 2, basis.boundary(side) and basis.boundaryBasis(side) must number the
        // side's functions identically and the anchors of the boundary basis
        // must be unisolvent: true for tensor and rational tensor bases, not
        // for hierarchical ones. A 1-D side is a point and needs neither.
        GISMO_ENSURE(internal::dirichletInterpolationSupported(u.source().basis(it->patch())),
                     "Dirichlet interpolation only implemented for tensor bases. Use `dirichlet::quasiInterpolation` (hierarchical bases) or `dirichlet::l2Projection` instead.");
    }

    for ( typename bcList::const_iterator it =  bc.begin("Dirichlet");
          it != bc.end("Dirichlet") ; ++it )
    {
        if( it->unknown()!=u.id() ) continue;
        internal::dirichletSideInterpolationDofs(u, bc, *it, fixedDofs);
    }
}

/// \brief Computes the Dirichlet dofs of \a u (its fixedPart()) for dirichlet::quasiInterpolation.
///
/// Per Dirichlet side of \a u, the side data is quasi-interpolated in the boundary basis
/// (internal::dirichletSideQuasiInterpolation) of a tensor-product B-spline, NURBS, THB,
/// rational THB or HB patch. Data in the span of the boundary basis is reproduced exactly.
/// The result is NOT interpolatory: a dof shared by two Dirichlet sides of a patch, or by
/// sides of two patches at an interface, keeps the value of the side processed last, on
/// tensor-product patches too; for data outside the span the sides differ at the
/// approximation-error level. On a 1-D patch the end dof is the data value at the end point.
/// Throws if a Dirichlet side of \a u lies on a patch with any other basis (use
/// dirichlet::l2Projection); the support of all patches is checked before any dof is written.
template<class T>
void gsDirichletValuesByQuasiInterpolation(const expr::gsFeSpace<T> & u,
                                           const gsBoundaryConditions<T> & bc)
{
    gsMatrix<T> & fixedDofs = const_cast<expr::gsFeSpace<T>&>(u).fixedPart();
    fixedDofs.setZero(u.mapper().boundarySize(), 1 );

    typedef gsBoundaryConditions<T> bcList;
    for ( typename bcList::const_iterator it =  bc.begin("Dirichlet");
          it != bc.end("Dirichlet") ; ++it )
    {
        if( it->unknown()!=u.id() ) continue;
        GISMO_ENSURE(internal::dirichletQuasiInterpolationSupported(u.source().basis(it->patch())),
                     "Dirichlet quasi-interpolation is not implemented for the basis of patch "
                     << it->patch() << " (supported: tensor-product B-spline, NURBS, THB, rational THB, HB). "
                     "Use `dirichlet::l2Projection` instead.");
    }

    for ( typename bcList::const_iterator it =  bc.begin("Dirichlet");
          it != bc.end("Dirichlet") ; ++it )
    {
        if( it->unknown()!=u.id() ) continue;
        internal::dirichletSideQuasiInterpolationDofs(u, bc, *it, fixedDofs);
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


/// \param sameElement asserts that each boundary quadrature batch lies in a single Bezier element of
/// the geometry map; passing false evaluates the map per point. Read from no option list.
template<class T>
void gsDirichletValuesByL2Projection( const expr::gsFeSpace<T> & u,
                                      const gsBoundaryConditions<T> & bc,
                                      const bool sameElement = true)
{
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

    unsigned mapFlags = NEED_MEASURE;
    if (sameElement) mapFlags |= SAME_ELEMENT;
    gsMapData<T> md(mapFlags);

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

        GISMO_ENSURE(bc.hasGeoMap(), "gsDirichletValuesByL2Projection: the boundary conditions "
                     "carry no geometry map, which the projection needs for the boundary "
                     "measure; call bc.setGeoMap(...).");
        const gsFunctionSet<T> & gmap = bc.geoMap();

        const index_t com = iter->unkComponent();
        const int patchIdx   = iter->patch();
        const gsBasis<T> & basis = u.source().basis(patchIdx);
        const gsFunction<T> & patch = gmap.function(patchIdx);

        // Set up quadrature to degree+1 Gauss points per direction,
        // all lying on iter->side() except from the direction which
        // is NOT along the element
        gsGaussRule<T> bdQuRule(basis, (T)1, 1, iter->side().direction());

        // Create the iterator along the given part boundary.
        typename gsBasis<T>::domainIter bdrIter    = basis.domain()->beginBdr(iter->side());
        typename gsBasis<T>::domainIter bdrIterEnd = basis.domain()->endBdr(iter->side());

        for (; bdrIter<bdrIterEnd; ++bdrIter)
        {
            bdQuRule.mapTo(bdrIter.lowerCorner(), bdrIter.upperCorner(),
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

                // if the component is not specified and the function evaluates
                // for all target dimensions simultaneous, rhsValues does not
                // need to be updated
                if ((com != -1) || (r == 0))
                {
                  // If the condition is homogeneous then fill with zeros
                  if (iter->isHomogeneous())
                  {
                    rhsVals.setZero((com==-1) ? u.dim() : 1, md.points.size());
                  }
                  else
                  {
                    if (iter->parametric())
                      rhsVals =
                          iter->function()->piece(patchIdx).eval(md.points);
                    else
                      rhsVals = iter->function()->piece(patchIdx).eval(
                          gmap.piece(patchIdx).eval(md.points));
                  }
                }

                GISMO_ASSERT((com!=-1) || rhsVals.rows() == u.dim(),
                    "If no component is specified for Dirichlet boundary, "
                    "target dimension must match field dimension.");
                GISMO_ASSERT((com==-1) || rhsVals.rows() == 1,
                    "If the component is specified for Dirichlet boundary, "
                    "then a scalar function is expected.");

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

                        globProjRhs.at(ii) += weight_k * basisVals(i, k) * rhsVals( (-1==com?r:0) ,k);
                        //globProjRhs.at(ii) += weight_k * basisVals(i, k) * rhsVals(r ,k);

                    } // for i
                } // for k
            }// for r
        } // bdrIter
    } // boundaryConditions-Iterator

    gsSparseMatrix<T> globProjMat(mapper.boundarySize(), mapper.boundarySize());
    globProjMat.setFrom(projMatEntries);
    globProjMat.makeCompressed();

    // Solve the linear system:
    // The position in the solution vector already corresponds to the
    // numbering by the boundary index. Hence, we can simply take them
    // for the values of the eliminated Dirichlet DOFs.
#ifdef GISMO_WITH_PARDISO
    typename gsSparseSolver<T>::PardisoLU solver;
#else
    typename gsSparseSolver<T>::CGDiagonal solver;
#endif
    fixedDofs = solver.compute(globProjMat).solve(globProjRhs);
} // computeDirichletDofsL2Proj

/// \brief Computes the Dirichlet dofs of \a u (its fixedPart()) for dirichlet::automatic.
///
/// If every patch carrying a Dirichlet side of \a u is supported by interpolation or by
/// quasi-interpolation, then per side: tensor-product/NURBS patch -> anchor interpolation (on a 1-D
/// patch, of any supported basis, the end-point rule), as
/// gsDirichletValuesByTPInterpolation (bitwise); other supported patch (THB, rational THB, HB)
/// -> quasi-interpolation, as gsDirichletValuesByQuasiInterpolation, with its last-write-wins
/// at shared dofs. Otherwise the whole unknown is computed by gsDirichletValuesByL2Projection
/// with \a sameElement. Prints nothing.
template<class T>
void gsDirichletValuesAutomatic(const expr::gsFeSpace<T> & u,
                                const gsBoundaryConditions<T> & bc,
                                const bool sameElement = true)
{
    typedef gsBoundaryConditions<T> bcList;
    for ( typename bcList::const_iterator it =  bc.begin("Dirichlet");
          it != bc.end("Dirichlet") ; ++it )
    {
        if( it->unknown()!=u.id() ) continue;
        const gsBasis<T> & basis = u.source().basis(it->patch());
        if ( !internal::dirichletInterpolationSupported(basis) &&
             !internal::dirichletQuasiInterpolationSupported(basis) )
        {
            gsDirichletValuesByL2Projection(u, bc, sameElement);
            return;
        }
    }

    gsMatrix<T> & fixedDofs = const_cast<expr::gsFeSpace<T>&>(u).fixedPart();
    fixedDofs.setZero(u.mapper().boundarySize(), 1 );

    for ( typename bcList::const_iterator it =  bc.begin("Dirichlet");
          it != bc.end("Dirichlet") ; ++it )
    {
        if( it->unknown()!=u.id() ) continue;
        if ( internal::dirichletInterpolationSupported(u.source().basis(it->patch())) )
            internal::dirichletSideInterpolationDofs(u, bc, *it, fixedDofs);
        else
            internal::dirichletSideQuasiInterpolationDofs(u, bc, *it, fixedDofs);
    }
}



}; // namespace gismo
