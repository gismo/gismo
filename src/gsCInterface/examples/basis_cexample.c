/** @file basis_example.c

    @brief Tests C interface

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): 
*/

#include <gsCInterface/Cgismo.h>

#include <stdio.h>
#include <string.h>

int main(int argc, char* argv[])
{
    double knots[6];
    knots[0]=knots[1]=knots[2]=0;
    knots[3]=knots[4]=knots[5]=1;

    gsCKnotVector * kv = gsKnotVector_create(knots,6);
    print(kv);

    gsCBasis * b = gsBSplineBasis_create(kv);    
    destroy(kv);
    printf("\n\n");
    print(b);
    printf("\n");

    double udata[4];
    udata[0]=0.00;
    udata[1]=0.33;
    udata[2]=0.66;
    udata[3]=1.00;
    gsCMatrix * u = gsMatrix_create_rcd(1,4,udata);
    gsCMatrix * result = gsMatrix_create();
    gsFunctionSet_eval_into(b, u, result);
    printf("Matrix with  %d rows and %d columns:\n", rows(result), cols(result) );
    print(result);
    printf("\n");
    double * data = data(result);
    printf("First row of that matrix: %f , %f, %f, %f\n", data[0], data[3], data[6], data[9]);

    gsCGeometry * g = gsBSpline_create(b,result);
    print(g);
    printf("\n");
        
    destroy(u);
    destroy(result);
    destroy(b);
    destroy(g);
    
    return 0;
}
