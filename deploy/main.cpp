/** @file main.cpp

    @brief Hello G+Smo test program

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): 
*/

#include <gismo.h>

using namespace gismo;

//int main(int argc, char* argv[])
int main()
{
    gsInfo <<  "Hello G+Smo.\n";

    real_t a = 2.0; // a real number, ie. double

    index_t b = 3; // an integer, ie. int

    GISMO_ASSERT( a*b == 6, "This is an error, 2*3 should be 6.");

    // creating the knot-vector: [0, 0, 0, 1, 2, 2, 2]
    gsKnotVector<> kv((real_t)0, a, 1, b);

    return 0;
}
