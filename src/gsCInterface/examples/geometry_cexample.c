/** @file geometry_cexample.c

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
    char fname[30];
    strcpy( fname, "surfaces/simple.xml" );
    printf("Reading %s\n", fname);

    gsCGeometry * G = (gsCGeometry *)gsCReadFile(fname);

    print(G);

    // define some interesting parameter pairs (u,v)

    double uv[] = { 0.0, 0.0,
                    0.1, 0.0,
                    0.5, 0.0,
                    0.5, 0.2,
                    0.5, 0.5,
                    0.5, 0.9,
                    1.0, 1.0 };
    int nRows = 2, nCols = 7;
    gsCMatrix * uvm = gsMatrix_create_rcd(nRows,nCols,uv);
    
    printf("Input #rows = %d, #cols = %d\n", rows(uvm), cols(uvm) );
    for (int irow=0; irow<rows(uvm); irow++)
    {
        if (irow==0) { printf("      u: "); }
        if (irow==1) { printf("      v: "); }
        for (int icol=0; icol<cols(uvm); icol++)
            printf(" %9.3f", uv[icol*nRows + irow]);
        printf("\n");
    }

    printf("Input #rows = %d, #cols = %d\n", nRows, nCols);

    // evaluate positions (x,y,z) at given parameter values

    gsCMatrix * out_p = gsMatrix_create();
    gsFunctionSet_eval_into(G, uvm, out_p);
    double* out_data = data(out_p);
    int out_rows = rows(out_p), out_cols = cols(out_p);
    
    printf("Positions : got #rows = %d, #cols = %d\n", out_rows, out_cols);
    for (int irow=0; irow<out_rows; irow++) {
        if (irow==0) { printf("      x: "); }
        if (irow==1) { printf("      y: "); }
        if (irow==2) { printf("      z: "); }
        for (int icol=0; icol<out_cols; icol++) {
            printf(" %9.3f", out_data[icol*out_rows + irow]);
        }
        printf("\n");
    }
    destroy(out_p);

    // evaluate first derivatives d(x,y,z)/du and d(x,y,z)/dv at given parameter values
    gsCMatrix * out_d = gsMatrix_create();
    gsFunctionSet_deriv_into(G, uvm, out_d);
    out_data = data(out_d);
    out_rows = rows(out_d);
    out_cols = cols(out_d);

    printf("Derivatives : got #rows = %d, #cols = %d\n", out_rows, out_cols);
    for (int irow=0; irow<out_rows; irow++) {
        if (irow==0) { printf("  dx/du: "); }
        if (irow==1) { printf("  dx/dv: "); }
        if (irow==2) { printf("  dy/du: "); }
        if (irow==3) { printf("  dy/dv: "); }
        if (irow==4) { printf("  dz/du: "); }
        if (irow==5) { printf("  dz/dv: "); }
        for (int icol=0; icol<out_cols; icol++) {
            printf(" %9.3f", out_data[icol*out_rows + irow]);
        }
        printf("\n");
        if (irow==1 | irow==3) { printf("\n"); }
    }


    gsCMatrix * tMat = gsMatrix_create();
    gsCVector * tVec = gsVector_create();
    gsCGeometryTransform * tr = gsGeometryTransform_create(G, tMat, tVec);
    destroy(tr);
    destroy(tMat);
    destroy(tVec);   

    gsCGeometry * isoparam = gsTensorBSpline2_slice(G, 0, 0.5);
    print(isoparam);
    destroy(isoparam);

        
    // tb.slice(dir, param, result);
    //equidistant sampling. arc length curve (+surface..)
    //
    
    
    destroy(uvm);
    destroy(out_d);
    destroy(G);

    return 0;
}
