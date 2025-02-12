/** @file domains_example.cpp

    @brief Tutorial on gsDomain class.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): J. Speh
*/

#include <iostream>
#include <gismo.h>
#include <gsDomain/gsPointDomain.h>

using namespace gismo;

int main(int argc, char* argv[])
{

    gsCmdLine cmd("Tutorial on gsDomain class.");

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    // ======================================================================
    // 1D domains
    // ======================================================================

    // The knot vector is 1D domain
    gsKnotVector<real_t> kv(0,1,3,2);
    gsInfo<<"The knot vector has "<<kv.numElements()<<" elements.\n";

    // Iteration over a domain's interior is performed by iterating from
    // domain.beginAll() to domain.endAll()
    gsInfo<<"Elements in the knot vector\nID\tlower\tupper\n";
    for (auto it = kv.beginAll(); it != kv.endAll(); ++it)
        gsInfo<<it.id()<<"\t"<<it.lowerCorner().transpose()<<"\t"<<it.upperCorner().transpose()<<"\n";

    // The domain of a B-spline basis is a knot vector
    gsBSplineBasis<real_t> basis(kv);
    gsInfo<<"The domain of a B-spline basis is a knot vector:\n"<<*basis.domain()<<"\n";

    // ======================================================================
    // 2D domains
    // ======================================================================
    gsTensorBSplineBasis<2,real_t> basis2D(kv,kv);
    gsInfo<<"The basis\n"<<basis2D<<"\nhas domain\n"<<*basis2D.domain()<<"\n";

    // const gsDomain<real_t> & domain2D = *basis2D.domain(); // DOES NOT WORK
    gsDomain<real_t>::Ptr domain2D = basis2D.domain();
    // Interior iteration is like before
    gsInfo<<"Interior elements of the 2D tensor domain\nID\tlower\tupper\n";
    for (auto it = domain2D->beginAll(); it != domain2D->endAll(); ++it)
        gsInfo<<it.id()<<"\t"<<it.lowerCorner().transpose()<<"\t"<<it.upperCorner().transpose()<<"\n";

    // In addition, boundary iteration is possible
    std::vector<boxSide> bdrs = {boundary::west,
                                 boundary::east,
                                 boundary::south,
                                 boundary::north};
    gsInfo<<"Boundary elements of the 2D tensor domain\nID\tside\tlower\tupper\n";
    for (auto bdr : bdrs)
        for (auto it = domain2D->beginBdr(bdr); it != domain2D->endBdr(bdr); ++it)
            gsInfo<<it.id()<<"\t"<<bdr<<"\t"<<it.lowerCorner().transpose()<<"\t"<<it.upperCorner().transpose()<<"\n";

    // We can do the same for Hierarchical domains
    gsTHBSplineBasis<2> thbBasis2D(basis2D);
    std::vector<index_t> box = {1,0,0,2,2};
    thbBasis2D.refineElements(box);
    gsInfo<<"The basis\n"<<thbBasis2D<<"\nhas domain\n"<<*thbBasis2D.domain()<<"\n";

    gsDomain<real_t>::Ptr thbDomain2D = thbBasis2D.domain();
    // Interior and boundary iteration

    gsInfo<<"Interior elements of the 2D THB domain\nID\tlower\tupper\n";
    for (auto it = thbDomain2D->beginAll(); it != thbDomain2D->endAll(); ++it)
        gsInfo<<it.id()<<"\t"<<it.lowerCorner().transpose()<<"\t"<<it.upperCorner().transpose()<<"\n";

    gsInfo<<"Boundary elements of the 2D THB domain\nID\tside\tlower\tupper\n";
    for (auto bdr : bdrs)
        for (auto it = thbDomain2D->beginBdr(bdr); it != thbDomain2D->endBdr(bdr); ++it)
            gsInfo<<it.id()<<"\t"<<bdr<<"\t"<<it.lowerCorner().transpose()<<"\t"<<it.upperCorner().transpose()<<"\n";

    // ======================================================================
    // Multi-domains
    // ======================================================================
    // We can also iterate over domains with multiple sub-domains
    gsMultiBasis<> mb;
    mb.addBasis(basis2D.clone());
    mb.addBasis(thbBasis2D.clone());
    gsDomain<real_t>::Ptr multiDomain2D = mb.domain();
    gsInfo<<"The multi-basis has "<<mb.nPieces()<<" pieces.\n";
    gsInfo<<"Its domain also has "<<multiDomain2D->nPieces()<<" pieces:   \n";
    for (size_t d = 0; d < multiDomain2D->nPieces(); ++d)
        gsInfo<<"* Piece "<<d<<":\n\t"<<*multiDomain2D->subdomain(d)<<"\n";

    // We can iterate over all elements in one go
    gsInfo<<"Interior elements of the multi-domain\nID\tpiece\tlower\tupper\n";
    for (auto it = multiDomain2D->beginAll(); it != multiDomain2D->endAll(); ++it)
        gsInfo<<it.id()<<"\t"<<it.patch()<<"\t"<<it.lowerCorner().transpose()<<"\t"<<it.upperCorner().transpose()<<"\n";

    // ======================================================================
    // Points as domains
    // ======================================================================
    // Let's make some points
    gsMatrix<real_t> points(2,50);
    points.setRandom();

    // We can use these points as a domain
    gsPointDomain<real_t> pd(points);

    // And iterate over them
    // Since the domain is a point domain, the lower and upper corners are the same
    // Therefore, we use gsDomain::centerPoint() to get the point
    gsInfo<<"Points in the point domain\nID\tcenter\n";
    for (auto it = pd.beginAll(); it != pd.endAll(); ++it)
        gsInfo<<it.id()<<"\t"<<it.centerPoint().transpose()<<"\n";

#ifdef _OPENMP
    // ======================================================================
    // Threading
    // ======================================================================
    // The domain iterators are thread-safe
    // We can use them in parallel loops

#pragma omp parallel
{
    const int nt  = omp_get_num_threads();
    const int tid = omp_get_thread_num();

    if (tid == 0)
        gsInfo<<"Parallel iteration with "<<nt<<" threads over the point domain\nID\tthread\tcenter\n";

#pragma omp for
    for (auto it = pd.beginAll(); it != pd.endAll(); ++it)
#pragma omp critical
        gsInfo<<it.id()<<"\t"<<tid<<"\t"<<it.centerPoint().transpose()<<"\n";
}
#endif


    return 0;
}



