/** @file gsBoehm.h

    @brief Boehm's algorithm for knot insertion

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): J. Speh
*/

#pragma once

#include <gsCore/gsLinearAlgebra.h> // just for index_t

namespace gismo
{


// =============================================================================

// Inserting knots in B-splines.

// =============================================================================


/// Performs insertion of multiple knot on "knots" and coefficients "coefs".
template<class T, class KnotVectorType, class Mat> 
void gsBoehm( 
    KnotVectorType & knots,
    Mat & coefs,
    T val,
    int r = 1,
    bool update_knots = true
    );
    

/// Performs knot insertion once on "knots" and coefficients "coefs".
template<class T, class KnotVectorType, class Mat>
void gsBoehmSingle( 
    KnotVectorType & knots,
    Mat & coefs,
    T val,
    bool update_knots = true
    );



/// Performs knot insertion once on "knots" and coefficients "coefs".
/// Iter points to the first affected knot (ie. knot is the result of findspan)
template<class T, class iter, class Mat>
void gsBoehmSingle( iter knot,   // Knot iterator
                    Mat & coefs, // Coefficients (p+1)
                    int p,       // degree
                    T val        // knot value to insert
    );


/// Performs knot refinement on "knots" and coefficients "coefs", by
/// adding new knots vals
/// NURBS book p.165 modified
template<class KnotVectorType, class Mat, class ValIt>
void gsBoehmRefine( KnotVectorType & knots,
                    Mat & coefs, //all Coefficients
                    int p,       // degree
                    ValIt valBegin, ValIt valEnd, // knot values to insert
                    bool update_knots = true);



// =============================================================================

// Inserting knots in tensor-product B-splines.

// =============================================================================



/// Performs a knot insertion on knots and recomputes coefficients.
///
/// \param knots - vector of knots in direction "direction"
/// \param coefs - coefficients (control points)
/// \param val - value of a knot we will insert
/// \param direction - in which direction (knot vector) we will insert knot
/// \param str - vector of strides
/// \param r - how many times we will insert knot
/// \param update_knots - if we update knots or not
///
/// \ingroup Nurbs
//  Algorithm is based on gsBoehm and THE NURBS BOOK.
//
//  Some tests (for 2D and 3D) are written in gsTensorBoehm_test.
//
//  Possible improvement: function builds new matrix for coefficient (11), so it
//  uses at least [2 * memory(coefs) + epsilon] memory.
template<typename T, typename KnotVectorType, typename Mat>
void gsTensorBoehm(
        KnotVectorType& knots,
        Mat& coefs,
        T val,
        int direction,
        gsVector<unsigned> str,
        int r = 1,
        bool update_knots = true);



/// Performs a knot refinement and recomputes coefficients.
///
/// \param knots - knot vector in direction "direction"
/// \param coefs - coefficients (control points)
/// \param direction - in which direction we will refine knots
/// \param str - vector of strides
/// \param valBegin - iterator pointing to the begining of the vector of the
///                       knots we want to insert
/// \param valEnd - iterator pointing to the end of the vector of the knots
///                     we want to insert
/// \param update_knots - if we should update "knots" or not
///
/// \ingroup Nurbs
//
// Algorithm is based on gsBoehmRefine and THE NURBS BOOK
//
// Some tests (for 2D and 3D) are written in gsTensorBoehm_test
//
// Possible improvement: function builds new matrix for coefficients (!!), so it
// uses at least [2 * memory(coefs) + epsilon] memory.
template <typename KnotVectorType, typename Mat, typename ValIt>
void gsTensorBoehmRefine(
        KnotVectorType& knots,
        Mat& coefs,
        int direction,
        gsVector<unsigned> str,
        ValIt valBegin,
        ValIt valEnd,
        bool update_knots = true);


/// @brief Local refinement algorithm.
///
/// We refine given coefficients (coefs) in given direction with corresponding
/// knots. Knots, we want to insert, can be accessed through iterators valBegin
/// and valEnd. Index is index of the first basis function, in given direction,
/// that corresponds to first coefficient in the coefs matrix. Through variable
/// update_knots we can set if we want to add inserted knots into knot vector
/// knots. Variables nmb_of_coefs, act_size_of_coefs and size_of_coefs describe
/// coefs matrix.
///
/// We can use bigger coefs matrix than needed. Act_size_of_coefs describes
/// actual size of the coefficients. (For example sizes of the edges of the cube
/// that coefs matrix represents in 3D.) Size_of_coefs is size of the
/// coefficients we will populate in this function call. Nmb_of_coefs presents
/// number of nonzero coefficients in coefs matrix.
/// 
/// \ingroup Nurbs
template <short_t d, typename KnotVectorType, typename Mat, typename ValIt>
void gsTensorBoehmRefineLocal(
        KnotVectorType& knots,
        const unsigned index,
        Mat& coefs,
        gsVector<index_t, d>& nmb_of_coefs,
        const gsVector<index_t, d>& act_size_of_coefs,
        const gsVector<index_t, d>& size_of_coefs,
        const unsigned direction,
        ValIt valBegin,
        ValIt valEnd,
        const bool update_knots);


/// @brief Inserts knot @a val such that multiplicity of a @a val in knot vector
///       is equal degree.
///
/// Function inserts the knot @a val into knot vector @a knots, such that
/// multiplicity of @a val becomes @a knots.degree(). Given coefficients
/// @a coefs must be local to the knot @a val (coefficients corresponding to
/// the active base functions at @a val). @a size_of_coefs are sizes of the
/// coefficients. @a direction is direction in which we perform a knot
/// insertion. Sometimes we dont need to insert a knots at all rows at
/// @a direction. With variables @a start and @a end, we can set
/// subcube of the @a coefs where we will perform a knot insertion.
/// 
/// This function should just be used for evaluation via knot insertion (not the
/// full coefficient matrix will be computed).
/// \ingroup Nurbs
template <short_t d, typename T, typename KnotVectorType, typename Mat>
void gsTensorInsertKnotDegreeTimes(
        const KnotVectorType& knots,
        Mat& coefs,
        const gsVector<index_t, d>& size_of_coefs,
        T val,
        const unsigned direction,
        gsVector<index_t, d>& start,
        gsVector<index_t, d>& end);



// -----------------------------------------------------------------------------
//
// Helper functions
//
// -----------------------------------------------------------------------------



namespace bspline
{

/// From strides and current position computes index (row) used in global
/// presentation of the coefficients
///
/// \param stride - vector of strides
/// \param position - current position [x0, x1, x2, x3, ...]
///
/// \return global index (row)
inline
int getIndex(const gsVector<unsigned>& stride,
                      const gsVector<int>& position)
{

    GISMO_ASSERT(stride.size() == position.size(),
                 "stide size is not equal to dimension of coefficients size");

    int ind = 0;
    unsigned n = stride.size();

    for (unsigned i = 0; i < n; ++i)
        ind += static_cast<int>(stride[i]) * position[i];

    return ind;
}


template<short_t d>
index_t getIndex(const gsVector<index_t, d>& stride,
                      const gsVector<index_t, d>& position)
{

    index_t ind = 0;

    for (short_t i = 0; i < d; ++i)
        ind += stride[i] * position[i];

    return ind;
}


/// Computes 2D vector of alpha values.
///
/// \param alpha - to be computed
/// \param knots - knot vector
/// \param value - knot we want to insert
/// \param r - how many times we want to insert knot value
/// \param k - value in [knots[k], knots[k + 1])
/// \param p - degree corresponding B-splines
/// \param s - multiplicity of knot value
///
/// \return void (note "returns" value through alpha)
template <typename T, typename KnotVectorType>
void computeAlpha(std::vector< std::vector<T> >& alpha,
                    const KnotVectorType& knots,
                    T value,
                    int r,
                    int k,
                    int p,
                    int s)
{

    for (int j = 1; j <= r; ++j)
    {
        int L = k - p + j;
        for (int i = 0; i <= p - j - s; ++i)
        {
            T knot_i = knots[L + i];
            alpha[j - 1][i] = (value - knot_i) / (knots[i + k + 1] - knot_i);
        }
    }

}


/// Computes last point of the cube, so that can be passed to nextCubePoint
/// function
///
/// \param stride - vector of strides
/// \param number_of_points - total number of coefficients (coef.rows())
/// \param last_point - throght this parameter function sets last point of the
/// cube
///
/// \return void (note: "returns" value throught last_point)
///
inline
void getLastIndex(const gsVector<unsigned>& stride,
                      const unsigned number_of_points,
                      gsVector<int>& last_point
                      )
{
    index_t stride_length = stride.size();
    for (index_t i = 0; i < stride_length; ++i)
    {
        if (i != stride_length - 1)
            // we need to subtrack 1, because we need actual last index
            last_point[i] = static_cast<int>(stride[i + 1] / stride[i]) - 1;
        else
            last_point[i] = static_cast<int>(number_of_points / stride[i]) - 1;
    }

}


/// Computes 2D vector of alpha values, and also computes new knot values and
/// saves them into nknots.
///
/// \param alpha - 2D vector; to be computed
/// \param nknots - new knot vector, to be computed
/// \param knots - knot vector
/// \param valBegin - iterator pointing to the begining of the vector of the
///                       knots we want to insert
/// \param valEnd - iterator pointing to the end of the vector of the knots
///                     we want to insert
/// \param sparse
/// \return void (note "returns" values through alpha and nknots)
template <typename T, typename KnotVectorType, typename ValIt,
          typename newKnotsType>
void computeTensorAlpha(std::vector< std::vector<T> >& alpha,
                        newKnotsType& nknots,
                        const KnotVectorType& knots,
                        ValIt valBegin,
                        ValIt valEnd,
                        bool sparse = false)
{
    const int a =  knots.iFind(*valBegin)     - knots.begin();
    const int b = (knots.iFind(*(valEnd - 1)) - knots.begin()) + 1;
    const int p = knots.degree(); // degree
    const int nik = std::distance(valBegin, valEnd); // number of inserted knots
    const int nk = knots.size();

    int i = b + p - 1;
    int k = b + p + nik - 1; // nik - 1 == r

    if (!sparse)
    {
        for (int j = 0; j <= a; j++)
            nknots[j] = knots[j];

        for (int j = b + p; j < nk; j++)
            nknots[j + nik] = knots[j];
    }


    for (int j = nik - 1; 0 <= j; --j)
    {
        const T newKnot = *(--valEnd); // X[j]

        // possible improvement: binary search instead of sequential
        while ((newKnot <= knots[i]) && (a < i))
        {
            nknots[k] = knots[i];
            k--;
            i--;
        }

        for (int ell = 1; ell <= p; ell++)
        {
            T alfa = nknots[k + ell] - newKnot;
            if (math::abs(alfa) != 0.0)
                alfa /= nknots[k + ell] - knots[i - p + ell];

            alpha[j][ell - 1] = alfa;
        }

        nknots[k] = newKnot;
        k--;
    }
}


/// Corrects values of a new stride according to the old stride with respect to
/// direction of knot refinment and number of knots we want to insert.
///
/// \param new_str - vector of new strides (will be modified)
/// \param str - vector of old strides
/// \param direction - direction of the refinement
/// \param r - number of knots we will insert
///
/// \return void (note it modifies "new_str")
inline
void correctNewStride(gsVector<unsigned>& new_str,
                      const gsVector<unsigned>& str,
                      const int direction,
                      const int r)
{
    if (direction + 1 != str.size() )
    {
        new_str[direction + 1] += new_str[direction] * r;
    }
    for (index_t i = direction + 2; i < str.size(); ++i)
    {
        new_str[i] = new_str[i - 1] * (str[i] / str[i - 1]);
    }
}


/// Computes last point of a cube. We can pass the point to nextCubePoint
/// function. Function computes last point for local algorithm.
///
/// @param size_of_coef size of the coefficients
/// @param last_point we compute last point into this variable
template<short_t d>
void getLastIndexLocal(const gsVector<index_t, d>& size_of_coef,
                  gsVector<index_t, d>& last_point)
{
    for (short_t i = 0; i < d; ++i)
    {
        last_point[i] = size_of_coef[i] - 1;
    }
}


/// Function builds strides from given sizes of the coefficients. We can use
/// computed strides for calculating indices corresponding to the size of the
/// coefficients
///
/// @param size_of_coefs size of the coefficients
/// @param strides we compute strides into this variable
template <short_t d>
void buildCoeffsStrides(const gsVector<index_t, d>& size_of_coefs,
                         gsVector<index_t, d>& strides)
{
    strides(0) = 1;
    for (short_t dim = 1; dim < d; dim++)
        strides(dim) = size_of_coefs(dim - 1) * strides(dim - 1);
}



/// Function computes necessary number of iterations for given number of the
/// coefficients. We will iterate through cube (its size is nmb_of_coefs) in
/// given direction.
///
/// @param nmb_of_coefs number of the coefficients
/// @param direction direction in which we iterate
///
/// @return number of iterations
inline
unsigned numberOfIterations(const gsVector<index_t>& nmb_of_coefs,
                            const unsigned direction)
{
    unsigned nmb_of_iter = 1;
    for (index_t i = 0; i < nmb_of_coefs.size(); ++i)
    {
        if (static_cast<unsigned>(i) != direction)
            nmb_of_iter *= nmb_of_coefs[i];
    }

    return nmb_of_iter;
}

} //end namespace bspline

} //end namespace gismo

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsBoehm.hpp)
#endif
