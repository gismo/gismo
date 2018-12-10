/** @file gsGlobs_test.cpp

    @brief This tests gsMultiBasis::getGlobs_withTransferMatrices

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): S. Takacs
 **/

#include "gismo_unittest.h"

SUITE(gsGlobs_test)
{

    TEST(test)
    {
        gsMultiPatch<> mp = gsNurbsCreator<>::BSplineLShapeMultiPatch_p2();
        gsMultiBasis<> mb(mp);

        mb.uniformRefine();
        
        gsBoundaryConditions<> bc;

        std::vector< std::vector< gsMultiBasis<>::glob > > result = mb.getGlobs_withTransferMatrices(bc);
        
        /*
        for (unsigned i=0; i<result.size(); ++i)
            for (unsigned j=0; j<result[i].size(); ++j)
            {
                gsInfo << "Consider glob [" << i << "][" << j << "]:\n";
                gsInfo << "Basis: " << result[i][j].basis.get() << "\n";
                gsInfo << "Corners: ";
                for (unsigned n=0; n<result[i][j].corners.size(); ++n )
                    gsInfo << "  " << result[i][j].corners[n];
                gsInfo << "\n";
                gsInfo << "Transfer: " << result[i][j].transfer.transpose() << "\n";
            }
        */

        CHECK(result.size() == 3);

        // Have 7 vertices with 1 dof each
        CHECK(result[0].size() == 8);
        for (unsigned i=0; i<8; ++i )
            CHECK( result[0][i].transfer.cols() == 1 && result[0][i].transfer.rows() == 40 );

        // Have 10 edges with 2 dofs each
        CHECK(result[1].size() == 10);
        for (unsigned i=0; i<10; ++i )
            CHECK( result[1][i].transfer.cols() == 2 && result[1][i].transfer.rows() == 40 );

        // Have 3 patches with 4 interior dofs each
        CHECK(result[2].size() == 3);
        for (unsigned i=0; i<3; ++i )
            CHECK( result[2][i].transfer.cols() == 4 && result[2][i].transfer.rows() == 40 );

    }

}
