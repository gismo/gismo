/** @file gsBoxTopology_test.cpp

    @brief Tests gsMultiPatch::computeTopology, gsBoxTopology::allComponents
    and gsDofMapper::markBoundary on grids of boxes against the definitions by
    linear search and by single eliminations

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
**/

#include "gismo_unittest.h"
#include <gsAssembler/gsDofMapperCreator.h>

namespace
{

// Grid of nx x ny x nz unit cubes, numbered lexicographically
gsMultiPatch<> boxGrid(index_t nx, index_t ny, index_t nz)
{
    gsMultiPatch<> mp;
    for (index_t k=0; k<nz; ++k)
        for (index_t j=0; j<ny; ++j)
            for (index_t i=0; i<nx; ++i)
                mp.addPatch( *gsNurbsCreator<>::BSplineCube(1.0, i, j, k) );
    return mp;
}

// Grid of nx x ny unit squares; every second patch has its parameter
// directions swapped and reversed, so the interfaces carry nontrivial
// direction maps and orientations
gsMultiPatch<> squareGridMixed(index_t nx, index_t ny)
{
    gsMultiPatch<> mp;
    for (index_t j=0; j<ny; ++j)
        for (index_t i=0; i<nx; ++i)
        {
            gsGeometry<>::uPtr g = gsNurbsCreator<>::BSplineSquare(1.0, i, j);
            if ((i+j)%2)
            {
                // swap the two parameter directions: the control points of
                // a bilinear patch are (0,0),(1,0),(0,1),(1,1)
                gsMatrix<> c = g->coefs();
                g->coefs().row(1) = c.row(2);
                g->coefs().row(2) = c.row(1);
            }
            mp.addPatch(*g);
        }
    return mp;
}

// The sorted canonic corners of a component, by the linear search of
// gsBoxTopology::getCornerList
std::vector<patchCorner> canonicCorners(const gsBoxTopology& bt, const patchComponent& pc)
{
    std::vector<patchCorner> crns = pc.containedCorners(), result, list;
    for (size_t l=0; l<crns.size(); ++l)
    {
        bt.getCornerList(crns[l], list);
        result.push_back(*std::min_element(list.begin(), list.end()));
    }
    std::sort(result.begin(), result.end());
    return result;
}

bool sameInterfaces(const gsBoxTopology& a, const gsBoxTopology& b)
{
    if (a.nInterfaces() != b.nInterfaces() || a.nBoundary() != b.nBoundary())
        return false;
    for (size_t i=0; i<a.nInterfaces(); ++i)
    {
        const boundaryInterface& x = a.interfaces()[i], &y = b.interfaces()[i];
        if (!(x.first() == y.first()) || !(x.second() == y.second())
            || x.dirMap() != y.dirMap() || x.dirOrientation() != y.dirOrientation())
            return false;
    }
    for (size_t i=0; i<a.nBoundary(); ++i)
        if (!(a.boundaries()[i] == b.boundaries()[i]))
            return false;
    return true;
}

// Every component holds the patch components with equal canonic corners, and
// different components have different canonic corners
bool componentsConsistent(const gsBoxTopology& bt, const std::vector< std::vector<patchComponent> >& comps)
{
    std::set< std::vector<index_t> > keys;
    size_t total = 0;
    for (size_t n=0; n<comps.size(); ++n)
    {
        const std::vector<patchCorner> c0 = canonicCorners(bt, comps[n][0]);
        std::vector<index_t> key;
        key.push_back(comps[n][0].dim());
        for (size_t l=0; l<c0.size(); ++l)
        {
            key.push_back(c0[l].patch);
            key.push_back(c0[l].m_index);
        }
        if (!keys.insert(key).second)
            return false;
        for (size_t i=0; i<comps[n].size(); ++i)
        {
            if (comps[n][i].dim() != comps[n][0].dim())
                return false;
            const std::vector<patchCorner> ci = canonicCorners(bt, comps[n][i]);
            if (ci.size() != c0.size() || !std::equal(ci.begin(), ci.end(), c0.begin()))
                return false;
        }
        total += comps[n].size();
    }
    // every patch has 3^d components
    return total == (size_t)(bt.nBoxes()*math::ipow(3, bt.dim()));
}

}

SUITE(gsBoxTopology_test)
{
    TEST(computeTopology_box_grid)
    {
        const index_t nx = 4, ny = 3, nz = 2;
        gsMultiPatch<> mp = boxGrid(nx, ny, nz);
        CHECK( mp.computeTopology() );
        const size_t nIf = (nx-1)*ny*nz + nx*(ny-1)*nz + nx*ny*(nz-1);
        CHECK( mp.nInterfaces() == nIf );
        CHECK( mp.nBoundary() == 6*nx*ny*nz - 2*nIf );

        // The comparison of the side centers prunes candidate pairs; the
        // comparison of the corners only (all pairs of sides) must give the
        // same interfaces in the same order
        gsMultiPatch<> ref = boxGrid(nx, ny, nz);
        CHECK( ref.computeTopology(1e-4, true) );
        CHECK( sameInterfaces(mp, ref) );
    }

    TEST(computeTopology_orientations)
    {
        gsMultiPatch<> mp = squareGridMixed(5, 4), ref = squareGridMixed(5, 4);
        CHECK( mp.computeTopology() );
        CHECK( ref.computeTopology(1e-4, true) );
        CHECK( mp.nInterfaces() == (size_t)(4*4 + 5*3) );
        CHECK( sameInterfaces(mp, ref) );
    }

    TEST(allComponents_box_grid)
    {
        const index_t nx = 4, ny = 3, nz = 2;
        gsMultiPatch<> mp = boxGrid(nx, ny, nz);
        mp.computeTopology();
        const std::vector< std::vector<patchComponent> > comps = mp.allComponents();

        // vertices, edges, faces and patches of the grid
        gsMatrix<index_t> numbers(4,1), expected(4,1);
        numbers.setZero();
        for (size_t n=0; n<comps.size(); ++n)
            numbers(comps[n][0].dim(),0) += 1;
        expected << (nx+1)*(ny+1)*(nz+1),
                    nx*(ny+1)*(nz+1) + (nx+1)*ny*(nz+1) + (nx+1)*(ny+1)*nz,
                    (nx+1)*ny*nz + nx*(ny+1)*nz + nx*ny*(nz+1),
                    nx*ny*nz;
        CHECK( numbers == expected );
        CHECK( componentsConsistent(mp, comps) );

        // an inner vertex is shared by 8 patches, an inner edge by 4
        size_t maxVertex = 0, maxEdge = 0;
        for (size_t n=0; n<comps.size(); ++n)
        {
            if (comps[n][0].dim() == 0) maxVertex = std::max(maxVertex, comps[n].size());
            if (comps[n][0].dim() == 1) maxEdge   = std::max(maxEdge,   comps[n].size());
        }
        CHECK( maxVertex == 8 );
        CHECK( maxEdge == 4 );
    }

    TEST(allComponents_orientations)
    {
        gsMultiPatch<> mp = squareGridMixed(5, 4);
        mp.computeTopology();
        const std::vector< std::vector<patchComponent> > comps = mp.allComponents();
        gsMatrix<index_t> numbers(3,1), expected(3,1);
        numbers.setZero();
        for (size_t n=0; n<comps.size(); ++n)
            numbers(comps[n][0].dim(),0) += 1;
        expected << 6*5, 5*5 + 6*4, 5*4;
        CHECK( numbers == expected );
        CHECK( componentsConsistent(mp, comps) );
    }
    TEST(markBoundary_equals_eliminateDof)
    {
        // Clamp the sides of minimal x of a 3x3x2 grid: the coupled dofs on
        // these sides are shared with other clamped and free patches
        gsMultiPatch<> mp = boxGrid(3, 3, 2);
        mp.computeTopology();
        gsMultiBasis<> mb(mp);
        mb.uniformRefine();
        gsDofMapper batched = createMapper(mb, gsBoundaryConditions<>(), dirichlet::elimination,
                                           iFace::glue, /*nComp=*/1, /*unk=*/0, /*finalize=*/false);
        gsDofMapper single = batched, allComp = batched;
        index_t eliminated = 0;
        for (gsMultiPatch<>::const_biterator it = mp.bBegin(); it != mp.bEnd(); ++it)
        {
            if (it->index() != boxSide(boundary::west).index())
                continue;
            const gsMatrix<index_t> bnd = mb[it->patch].boundary(it->side());
            batched.markBoundary(it->patch, bnd, 0);
            allComp.markBoundary(it->patch, bnd, -1);
            for (index_t i = 0; i < bnd.rows(); ++i)
                single.eliminateDof(bnd(i,0), it->patch, 0);
            eliminated += bnd.rows();
        }
        CHECK( eliminated > 0 );
        batched.finalize();
        single.finalize();
        allComp.finalize();
        CHECK( allComp.asVector() == single.asVector() );
        CHECK( batched.freeSize() == single.freeSize() );
        CHECK( batched.boundarySize() == single.boundarySize() );
        CHECK( batched.coupledSize() == single.coupledSize() );
        CHECK( batched.boundarySize() > 0 );
        CHECK( batched.asVector() == single.asVector() );
    }
}
