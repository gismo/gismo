/** @file gsExprAssembler_test.cpp

    @brief Tests for gsExprAssembler

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): R. Schneckenleitner
*/

#include "gismo_unittest.h"


SUITE(gsExprAssembler_test)
{
         TEST(InterfaceExpression)
                {
                    const index_t numRef = 2;
                    gsVector<> translation(2);
                    translation << -0.5, -1;
                    gsFunctionExpr<> ff("if(x>0,1,-1)", 2); // function with jump
                    gsMultiPatch<> patches = gsNurbsCreator<>::BSplineSquareGrid(1,2,1);
                    patches.patch(0).translate(translation);
                    patches.patch(1).translate(translation);

                    gsMultiBasis<> mb(patches);

                    for(index_t i = 0; i < numRef; i++)
                        mb.uniformRefine();

                    gsExprEvaluator<> ev;
                    ev.setIntegrationElements(mb);
                    gsExprEvaluator<>::geometryMap G = ev.getMap(patches);
                    gsExprEvaluator<>::variable f = ev.getVariable(ff, G);

                    ev.integral(f);

                    const real_t v = ev.value();
                    CHECK( v*v < 1e-10 );
                }
        }
