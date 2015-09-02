/** @file gsPardiso.h

    @brief Header file for using the PARDISO library (http://www.pardiso-project.org)

    Free academic license for PARDISO can be obtained at http://www.pardiso-project.org/branch

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): A. Mantzaflaris
*/

#pragma once

extern "C" 
{

// For Eigen
typedef void *        _MKL_DSS_HANDLE_t;
#ifdef MKL_ILP64
typedef long long int _INTEGER_t;//64-bit integer type
#else
typedef int           _INTEGER_t;// standard int
#endif

// Declaration of pardiso
void pardiso( void *, int *, int *, int *, int *, int * , void *, int *,
              int * , int *, int *, int *, int *, void *, void *, int *);
 
// Declaration of pardiso_64
void pardiso_64( void          *, long long int *, long long int *, long long int *,
                 long long int *, long long int *,          void *, long long int *,
                 long long int *, long long int *, long long int *, long long int *,
                 long long int *,          void *,          void *, long long int *);

} // extern "C"
