/** @file gsC1SplineBase.h

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
class gsC1SplineBase
{

public:

    virtual ~gsC1SplineBase() {};

public:
    // To be overwritten in inheriting classes

    /**
     * @brief      Returns the smoothing matrix into \a matrix
     *
     * @param      matrix  The matrix
     */
    virtual void matrix_into(gsSparseMatrix<T> & matrix) const = 0;

    /**
     * @brief      Returns the basis that is used.
     *
     * @return     { description_of_the_return_value }
     */
    virtual gsMultiBasis<T> localBasis() const =0;


// implementations of gsBasis
public:
    // Functions that do not depend on implementation

    /**
     * @brief      Returns the smoothing matrix
     *
     * The number of columns of the matrix corresponds to the number of basis functions in the local basis; this is the sum of all the basis functions over all the patches.
     * The number of rows of the matrix corresponds to the number of global basis functions, i.e. the number of basis functions corresponding to the D-Patch.
     * Multiplying the basis with the local basis function values gives the values of the basis functions in the global basis.
     *
     * @return     A matrix \a result to transfer local to global coefficients
     */
    gsSparseMatrix<T> matrix() const
    {
        gsSparseMatrix<T> result; matrix_into(result);
        return result;
    }

    gsOptionList options() {return m_options;}
    void defaultOptions();
    void getOptions();

private:
    // Helper functions

public:
    gsSparseMatrix<T> m_matrix;

protected:
    // Data members

    // Put here the members of the shared functions


};

}

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsC1SplineBase.hpp)
#endif
