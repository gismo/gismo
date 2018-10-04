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
                    gsFunctionExpr<> ff("if(x>0,1,-1)", 2); // function with jump 2
                    //gsMultiBasis<T> mb;
                    //gsBSplineBasis<> b1( );
                    gsExprEvaluator<> ev;
                    //ev.setIntegrationElements(mb);
                    gsExprEvaluator<>::variable f = ev.getVariable(ff);

                    //ev.integralInterface( avg(f) );

                    const real_t v = ev.value();
                    CHECK_CLOSE( 2.0, v, 1e-9);
                }
        }
