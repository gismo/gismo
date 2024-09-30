/** @file gsBasisComponents_test.cpp

    @brief test boxComponent

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): S. Takacs

**/

#include "gismo_unittest.h"

using namespace gismo;


SUITE(gsBasisComponents_test)
{
    TEST(main_test)
    {
        gsMultiPatch<>::uPtr mp(gsReadFile<>("volumes/fichera.xml"));
        gsMultiBasis<> mb(*mp);

        CHECK(mb.nBases() == 7);

        // Get all components
        std::vector< std::vector<patchComponent> > components = mb.topology().allComponents();

        // Check dimensions
        gsMatrix<index_t> numbers(4,1);
        numbers.setZero(4,1);
        for (size_t i=0; i<components.size(); ++i)
        {
            for (size_t j=0; j<components[i].size(); ++j)
                CHECK(components[i][0].dim() == components[i][j].dim());

            numbers( components[i][0].dim(), 0 ) += 1;
        }
        gsMatrix<index_t> numbers_check(4,1);
        numbers_check << 26,51,33,7;
        CHECK( numbers == numbers_check );

        // Check indices
        mb.uniformRefine();
        gsDofMapper dm;
        mb.getMapper(
            dirichlet::elimination,
            iFace::glue,
            gsBoundaryConditions<>(),
            dm,
            0
        );
        const index_t nTotalDofs = dm.freeSize();
        CHECK( nTotalDofs == 117 );
        gsMatrix<index_t> globalIndices = gsMatrix<index_t>::Constant(nTotalDofs,1,-1);
        index_t sz = components.size();
        for (index_t i=0; i<sz; ++i)
        {
            gsMatrix<index_t> indices;
            std::vector<gsBasis<>::uPtr> bases = mb.componentBasis_withIndices(components[i],dm,indices,true);
            for (index_t j=0; j<indices.rows(); ++j)
            {
                index_t l = indices(j,0);
                CHECK( l<nTotalDofs );
                CHECK( globalIndices( l, 0 ) == -1 ); // Index has not been assigned so far
                globalIndices( l, 0 ) = i; // Assign
            }
        }

        gsMatrix<index_t> globalIndices_check(117,1);
        globalIndices_check << 42,8,46,10,0,14,95,56,96,57,26,58,3,24,29,67,27,62,69,102,72,20,33,5,75,31,105,61,77,28,
            106,70,9,47,1,15,78,108,36,81,4,25,30,68,38,85,86,112,34,6,76,32,88,39,114,87,115,44,74,11,90,13,35,2,116,
            80,89,37,91,40,92,41,7,45,93,48,94,49,16,53,97,59,98,50,17,54,99,63,100,19,64,60,101,103,51,71,21,104,65,43,
            12,52,18,55,107,79,109,82,110,22,66,83,111,73,23,113,84;

        CHECK( globalIndices == globalIndices_check );

    }

}
