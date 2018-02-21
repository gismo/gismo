/** @file gsKronecker_test.cpp

    @brief Tests the Kronecker operators and the Kronecker products

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): C. Hofreither, A. Manzaflaris, S. Takacs
*/

#include "gismo_unittest.h"

SUITE(gsKronecker_test)
{

    TEST(testKroneckerOp)
    {

        gsMatrix<> A(3,3), B(3,3);

        A << 3, 7, 3,
            2, -2, 1,
            9, 5, 2;

        B << 6, 1, 4,
            -3, 4, 2,
            1, 2, 5;

        gsMatrix<> C(9,9);
        C <<    // kron(A,B) computed in Matlab
            18,    3,   12,   42,    7,   28,   18,    3,   12,
            -9,    12,   6,  -21,   28,   14,   -9,   12,    6,
            3,     6,   15,    7,   14,   35,    3,    6,   15,
            12,    2,    8,  -12,   -2,   -8,    6,    1,    4,
            -6,    8,    4,    6,   -8,   -4,   -3,    4,    2,
            2,     4,   10,   -2,   -4,  -10,    1,    2,    5,
            54,    9,   36,   30,    5,   20,   12,    2,    8,
            -27,   36,  18,  -15,   20,   10,   -6,    8,    4,
            9,     18,  45,    5,   10,   25,    2,    4,   10;

        gsMatrix<> x(9,1), y, y2;
        x << 1,2,3,4,5,6,7,8,9;

        gsKroneckerOp<> kron( makeMatrixOp(A), makeMatrixOp(B) );
        kron.apply(x, y);

        // compute Kronecker product directly and compare
        y2 = C * x;

        const real_t err = (y-y2).norm();
        
        CHECK (err <=  pow(10.0, - REAL_DIG ) );
    }

    TEST(testKroneckerProducts)
    {
        gsMatrix<> A(3,4), B(2,3), C(5,9);

        A.setRandom();
        B.setRandom();
        C.setRandom();

        gsKroneckerOp<> kron( makeMatrixOp(A), makeMatrixOp(B), makeMatrixOp(C) );
            
        gsMatrix<> x(4*3*9,3), y;
        x.setRandom();
        
        kron.apply(x, y);   

        // compute Kronecker product directly and compare (ColMajor)
        {
            gsSparseMatrix<> Asp = A.sparseView(), Bsp = B.sparseView(), Csp = C.sparseView();
            gsSparseMatrix<> K = getKroneckerProduct( getKroneckerProduct( Asp, Bsp ), Csp );
            gsMatrix<> y2 = K * x;

            const real_t err = (y-y2).norm();

            CHECK (err <=  pow(10.0, - REAL_DIG+2 ) );
        }

        // compute Kronecker product directly and compare (RowMajor)   
        {
            gsSparseMatrix<real_t, RowMajor> Asp = A.sparseView(), Bsp = B.sparseView(), Csp = C.sparseView();
            gsSparseMatrix<real_t, RowMajor> K = getKroneckerProduct( getKroneckerProduct( Asp, Bsp ), Csp );
            gsMatrix<> y2 = K * x;

            const real_t err = (y-y2).norm();
            CHECK (err <=  pow(10.0, - REAL_DIG+2 ) );
        }
      
        // compute Kronecker for directly and compare (dense)
        {
            gsMatrix<> K = getKroneckerProduct( getKroneckerProduct( A, B ), C );
            gsMatrix<> y2 = K * x;

            const real_t err = (y-y2).norm();

            CHECK (err <=  pow(10.0, - REAL_DIG+2 ) );
        }
    }

}


