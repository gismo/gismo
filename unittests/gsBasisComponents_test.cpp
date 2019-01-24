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
        gsBoundaryConditions<> bc;
        const index_t nTotalDofs = 117;
        gsMatrix<index_t> globalIndices = gsMatrix<index_t>::Constant(nTotalDofs,1,-1);
        index_t sz = components.size();
        for (index_t i=0; i<sz; ++i)
        {
            gsMatrix<unsigned> indices;
            std::vector<gsBasis<>::uPtr> bases = mb.componentBasis_withIndices(components[i],bc,gsOptionList(),indices,true);
            for (index_t j=0; j<indices.rows(); ++j)
            {
                index_t l = indices(j,0);
                CHECK( l<nTotalDofs );
                CHECK( globalIndices( l, 0 ) == -1 ); // Index has not been assigned so far
                globalIndices( l, 0 ) = i; // Assign
            }
        }

        gsMatrix<index_t> globalIndices_check(117,1);
        globalIndices_check << 42,8,45,10,0,12,95,54,96,55,21,56,2,19,23,63,22,59,65,102,68,26,28,4,71,30,105,73,74,32,
            106,75,9,46,1,13,76,108,33,78,3,20,24,64,34,81,82,112,29,5,72,31,84,35,114,85,115,86,87,36,88,37,38,6,116,
            89,90,39,91,40,92,41,7,44,93,47,94,48,14,51,97,57,98,49,16,52,99,60,100,17,61,58,101,103,66,67,25,104,70,
            43,11,50,15,53,107,77,109,79,110,18,62,80,111,69,27,113,83;
        CHECK( globalIndices == globalIndices_check );

    }

}
