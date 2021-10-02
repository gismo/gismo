/** @file gsApproxC1Spline.h

    @brief Provides declaration of Basis abstract interface.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): F. Buchegger
*/

/*
    TO DO
 */

#pragma once

#include<gsIO/gsOptionList.h>

namespace gismo
{

template<short_t d,class T>
class gsApproxC1Spline : public gsC1SplineBase<d,T>
{

public:

    gsApproxC1Spline(const gsMultiPatch<T> & patches)
    :
    {

    };

public:
    // To be overwritten in inheriting classes

    void compute();

    /**
     * @brief      Returns the smoothing matrix into \a matrix
     *
     * @param      matrix  The matrix
     */
    void matrix_into(gsSparseMatrix<T> & matrix) const;

    /**
     * @brief      Returns the basis that is used.
     *
     * @return     { description_of_the_return_value }
     */
    gsMultiBasis<T> localBasis() const;

private:
    // Helper functions

protected:
    // Data members

    // Put here the members of the shared functions


};

}

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsApproxC1Spline.hpp)
#endif
