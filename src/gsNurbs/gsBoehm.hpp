/** @file gsBoehm.hpp

    @brief Boehm's algorithm for knot insertion

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): J. Speh
*/

#pragma once

#include <gsUtils/gsCombinatorics.h>

namespace gismo
{

using gismo::bspline::getIndex;
using gismo::bspline::computeAlpha;
using gismo::bspline::getLastIndex;
using gismo::bspline::computeTensorAlpha;
using gismo::bspline::correctNewStride;


template<class T, class KnotVectorType, class Mat>
void gsBoehm(
    KnotVectorType & knots,
    Mat & coefs,
    T val,
    int r,
    bool update_knots
    )
{

    GISMO_ASSERT( r >= 1, "Must insert at least one knot." );
    if (r==1)
        return gsBoehmSingle(knots, coefs, val, update_knots);

    GISMO_ASSERT( coefs.rows() == (index_t)(knots.size() - knots.degree()-1),
                  "Incompatible coefficients("<<coefs.rows()
                  <<")/knots("<<knots.size()<<")/degree("<<knots.degree()<<")" ) ;

    const int p = knots.degree();
    typename KnotVectorType::uiterator kit = knots.uFind(val);
    const int k = kit.lastAppearance();
    // current multiplicity of val
    const int s = (*kit == val || *(++kit)==val ? kit.multiplicity() : 0 );
    /*
      typename KnotVectorType::uiterator kit = knots.uFind(val);
      int k = knots.iFind(val) - knots.begin();
      int s = knots.multiplicity(val); // current multiplicity of val
    */

    GISMO_ASSERT( s + r < p + 2  , "Multiplicity can be at most deg+1 ("<<p+1<<")" );
    int np= coefs.rows()-1;

    Mat tmp = coefs.middleRows(k-p, p+1);
    // resize coefficient matrix
    coefs.conservativeResize( coefs.rows()+r, coefs.cols() );

    // shift control points that are not affected
    for( index_t i = np; i>=k-s; --i )
        coefs.row(i+r) = coefs.row(i);

    // Compute new control points
    T a;
    int L = 0;
    for( index_t j = 1; j<=r; ++j )
    {
        L = k - p + j;

        for( index_t i = 0; i<=p-j-s; ++i )
        {
            a = (val - knots[L+i]) / (knots[i+k+1] - knots[L+i]);
            tmp.row(i) = a * tmp.row(i+1) + (1.0-a) * tmp.row(i);
        }
        coefs.row(L)= tmp.row(0);
        coefs.row(k+r-j-s)= tmp.row(math::max(p-j-s,(index_t)0));
    }
    for( index_t i = L+1; i<k-s; ++i )
        coefs.row(i) = tmp.row(i-L);
    //coefs.middleRows(L+1, k-s-L-1) = tmp.

    // Update knot vector
    if ( update_knots )
        knots.insert(val,r);
}



template<class T, class KnotVectorType, class Mat>
void gsBoehmSingle(
    KnotVectorType & knots,
    Mat & coefs,
    T val,
    bool update_knots
    )
{

    GISMO_ASSERT( coefs.rows() == (index_t)(knots.size() - knots.degree()-1),
                  "Incompatible coefficients/knots" ) ;

    int k = knots.iFind(val)-knots.begin();
    int p = knots.degree();

    coefs.duplicateRow( k );
    // // resize coefficient matrix
    //coefs.conservativeResize(coefs.rows() + 1, coefs.cols());
    //
    //for (index_t i = coefs.rows() - 1; i >= k+1; --i)
    //  coefs.row(i) = coefs.row(i-1);

    // Compute new control points
    T a;
    for( index_t i = k; i>=k-p+1; --i )
    {
        a = (val - knots[i]) / (knots[i+p] - knots[i]);
        coefs.row(i) = (T(1)-a) * coefs.row(i-1) + a * coefs.row(i);
    }

    // Update knot vector
    if ( update_knots )
        knots.insert(val);
}



template<class T, class iter, class Mat>
void gsBoehmSingle( iter knot,   // Knot iterator
                    Mat & coefs, // Coefficients (p+1)
                    int p,       // degree
                    T val)       // knot value to insert
{
    GISMO_ASSERT( coefs.rows() == p+1,
                  "Incompatible coefficients/knots" ) ;

    knot++;

    // resize coefficient matrix by adding one coefficient
    coefs.conservativeResize( p+2, coefs.cols() );

    // Shift coefficients by one place
    // for (index_t i = p+1; i >= 1; --i)
    //     coefs.row(i) = coefs.row(i-1);
     coefs.row(p+1) = coefs.row(p);

    // Compute new control points
    T a;
    for( index_t i = p ; i>=1; --i )
    {
        a = (val - *knot) / ( *(knot+p) - *knot );
        coefs.row(i) = (T(1)-a) * coefs.row(i-1) + a * coefs.row(i);
        knot++;
    }
}



template<class KnotVectorType, class Mat, class ValIt>
void gsBoehmRefine( KnotVectorType & knots,
                    Mat & coefs, //all Coefficients
                    int p,       // degree --remove
                    ValIt valBegin, ValIt valEnd, // knot values to insert
                    bool update_knots)
{
    if ( valBegin == valEnd ) return;

    typedef typename std::iterator_traits<ValIt>::value_type T;

    GISMO_ASSERT( (*valBegin >= knots[p]), "Value is before first knot: "
                  << *valBegin <<" < " <<knots[p] );
    // && (val[val.size()-1]<=knots[knots.size()-p-1]));
    //assert( knots[knots.size()-p-1]<=val[val.size()-1] );
    const int np= coefs.rows();
    const int nk= std::distance( valBegin, valEnd );
    coefs.conservativeResize(np+nk, coefs.cols());

    const int a =  knots.iFind( *valBegin    ) - knots.begin();
    const int b = (knots.iFind( *(valEnd-1)  ) - knots.begin()) + 1;
    //const int a = knots.uFind(*valBegin).lastAppearance();
    //const int b = knots.uFind(*(valEnd-1)).lastAppearance() + 1;

    //gsKnotVector<T> nknots(p, knots.size()+nk);
    std::vector<T> nknots(knots.size()+nk);

    // shift control points that are not affected
    for(index_t j = np ; j > b-1; j--)
      coefs.row( j+nk-1) = coefs.row(j-1);

    //std::copy( knots.begin(), knots.begin()+a+1,nknots.begin());
    for(int j = 0; j <= a; j++)
        nknots[j] = knots[j];
    //std::copy( knots.begin()+b+p, knots.bend(),nknots.begin()+nk+b+p);
    for(size_t j = b+p; j < knots.size(); j++)
        nknots[j + nk] = knots[j];

    int i = b + p - 1;
    int k = b + p + nk-1;

    for (int j = nk-1; j>=0; j--)
    {
        const T newKnot = *(--valEnd);

        while( (newKnot <= knots[i]) && (i>a) )
        {
            coefs.row(k-p-1) = coefs.row(i-p-1);

            nknots[k] = knots[i];
            k-- ;
            i-- ;
        }

        coefs.row(k-p-1) = coefs.row(k-p);


        for(int l = 1 ; l <=p; l++)
        {
            const int ind = k-p+l;

            T alfa = nknots[k+l] - newKnot;

            if( math::abs(alfa) == T(0.0) )
                coefs.row(ind-1) = coefs.row(ind);
            else
            {
                alfa /= nknots[k+l]-knots[i-p+l];
                coefs.row(ind-1) = alfa * coefs.row(ind-1) +
                (1.0-alfa)*coefs.row(ind);
            }
        }
        nknots[k] = newKnot;
        k--;
    }

    if ( update_knots )
        knots = KnotVectorType(p, nknots.begin(), nknots.end());
        //knots.insert(valBegin, valEnd); // bug ?
}


template<typename T, typename KnotVectorType, typename Mat>
void gsTensorBoehm(
        KnotVectorType& knots,
        Mat& coefs,
        T val,
        const int direction,
        gsVector<unsigned> str,
        int r,
        bool update_knots)
{
    GISMO_ASSERT(1 <= r, "Can not make insertion r < 1 times");
    GISMO_ASSERT(direction < str.size(),
                 "We can not insert a knot in a given direction");
    GISMO_ASSERT(knots.first() < val && val < knots.last(),
                 "We can not insert a knot outside of the knot vector interval");

    int d = str.size(); // dimension
    int k = knots.iFind(val) - knots.begin();
    int s = knots.multiplicity(val);
    int p = knots.degree();


    GISMO_ASSERT(s + r <= p,
                 "Multiplicity of a knot must be lower (or equal) "
                 "than a degree. Otherwise, we get non-continuous function.");

    // we will compute new coefficients and put them new_coef matrix
    int num_of_points = knots.size() - knots.degree() - 1;
    unsigned points_to_add = (coefs.rows() / num_of_points) * r;
    Mat new_coef = Mat(coefs.rows() + points_to_add, coefs.cols());

    // we precompute alpha (variable we use in knot insertion algorithm)
    std::vector< std::vector<T> > alpha(r, std::vector<T> (p - s));
    computeAlpha<T, KnotVectorType>(alpha, knots, val, r, k, p, s);

    // temporary matrix, for computation
    Mat tmp(p + 1, coefs.cols());

    // compute stride for new coefficients
    gsVector<unsigned> new_str(str);
    correctNewStride(new_str, str, direction, r);

    // this two vectors will hold indices of current points
    gsVector<int> position(d);
    position.fill(0);
    gsVector<int> new_position(position);

    // necessary for computation of the indices
    gsVector<int> first_point(position);
    gsVector<int> last_point(d);
    gsVector<int> new_last_point(d);
    getLastIndex(str, coefs.rows(), last_point);
    getLastIndex(new_str, new_coef.rows(), new_last_point);
    last_point[direction] = 0;
    new_last_point[direction] = 0;

    int ind, new_ind; // actual indices (row of a coefs, new_coef)
    unsigned step = str[direction];
    unsigned new_step = new_str[direction];
    bool flag = true; // safety flag

    do
    {
        if (!flag)
            GISMO_ERROR("This should not happened!"
                        "We do not have an index for the new matrix.");

        ind = getIndex(str, position);
        new_ind = getIndex(new_str, new_position);

        // copy untouched points from the begining
        for (int i = 0; i < k - p + 1; ++i)
        {
            new_coef.row(new_ind + i * new_step) = coefs.row(ind + i * step);
        }

        // make a temporary array of points
        int tmp_ind = ind + step * (k - p);
        for (int i = 0; i < p + 1; ++i)
        {
            tmp.row(i) = coefs.row(tmp_ind + step * i);
        }

        // compute new control points
        int L = 0;
        for (int j = 1; j <= r; ++j)
        {
            L = k - p + j;
            for (int i = 0; i <= p - j - s; ++i)
            {
                T a = alpha[j - 1][i];
                tmp.row(i) = a * tmp.row(i + 1) + (1.0 - a) * tmp.row(i);
            }

            new_coef.row(new_ind + new_step * L) = tmp.row(0);
            new_coef.row(new_ind + new_step * (k + r - j - s)) =
                        tmp.row(p - j - s);
        }

        // put new control point into proper position
        for (int i = L + 1; i < k - s; ++i)
        {
            new_coef.row(new_ind + step * i) = tmp.row(i - L);
        }

        // copy untouched points from the end
        for (int i = k - s; i < num_of_points; ++i)
        {
            new_coef.row(new_ind + (i + r) * step) = coefs.row(ind + i * step);
        }

        flag = nextCubePoint<gsVector<int> >(new_position, first_point,
                                             new_last_point);

    } while(nextCubePoint<gsVector<int> >(position, first_point, last_point));


    coefs = give(new_coef);

    if (update_knots)
    {
        knots.insert(val, r);
    }
}



template <typename KnotVectorType, typename Mat, typename ValIt>
void gsTensorBoehmRefine(
        KnotVectorType& knots,
        Mat& coefs,
        int direction,
        gsVector<unsigned> str,
        ValIt valBegin,
        ValIt valEnd,
        bool update_knots)
{

    typedef typename std::iterator_traits<ValIt>::value_type T;

    const int npts = coefs.rows(); // number of points
    const int nik = std::distance(valBegin, valEnd); // number of inserted knots
    const int nk = knots.size(); // number of knots
    const int p = knots.degree(); // degree
    const int d = str.size();    // dimension

    GISMO_ASSERT(knots[p] <= *valBegin && *(valEnd - 1) <= knots[nk - p - 1],
                 "Can not insert knots, they are out of the knot range");
    GISMO_ASSERT(direction < d,
                 "We can not insert a knot in a given direction");

    const int a =  knots.iFind(*valBegin)     - knots.begin();
    const int b = (knots.iFind(*(valEnd - 1)) - knots.begin()) + 1;

    // we compute new coefficients and put them into new_coef
    int npts_in_dir = nk - p - 1; // number of points in direction d
    int pts_to_add = (npts / npts_in_dir) * nik; // number of point we must add
    Mat new_coefs = Mat(npts + pts_to_add, coefs.cols());

    // allocate a memory for new knots and new control points
    std::vector<T> nknots(nk + nik);

    // compute stride for new coefficients
    gsVector<unsigned> new_str(str);
    correctNewStride(new_str, str, direction, nik);

    gsVector<int> position(d); // position old points
    for (int i = 0; i < d; ++i)
        position[i] = 0;

    gsVector<int> new_position(position); // position new points
    gsVector<int> first_point(position); // first point of a cube
    gsVector<int> last_point(d); // last point of a cube (old points)
    gsVector<int> new_last_point(d); // last point of a cube (new points)

    getLastIndex(str, npts, last_point);
    getLastIndex(new_str, npts + pts_to_add, new_last_point);

    last_point[direction] = 0;
    new_last_point[direction] = 0;

    int ind, new_ind;
    const int step = str[direction];
    const int new_step = new_str[direction];
    bool flag = true;

    // precompute alpha
    std::vector< std::vector<T> > alpha(nik, std::vector<T> (p));
    computeTensorAlpha<T, KnotVectorType, ValIt, std::vector<T> >
            (alpha, nknots, knots, valBegin, valEnd);

    do
    {
        ValIt valEndCopy = valEnd;

        if (!flag)
            GISMO_ERROR("This should not happened! "
                        "We do not have an index for the new matrix.");

        ind = getIndex(str, position);
        new_ind = getIndex(new_str, new_position);

        // copy control points that are not affected
        for (int j = 0; j <= a - p; ++j)
            new_coefs.row(new_ind + j * new_step) = coefs.row(ind + j * step);
        for (int j = npts_in_dir; b - 1 < j; --j)
            new_coefs.row(new_ind + (j + nik - 1) * new_step) =
                    coefs.row(ind + (j - 1) * step);

        // algorithm
        int i = b + p - 1;
        int k = b + p + nik - 1; // nik - 1 == r

        for (int j = nik - 1; 0 <= j; j--)
        {
            const T newKnot = *(--valEndCopy);
            while ((newKnot <= knots[i]) && (a < i))
            {
                new_coefs.row(new_ind + (k - p - 1) * new_step) =
                        coefs.row(ind + (i - p - 1) * step);
                k--;
                i--;
            }

            new_coefs.row(new_ind + (k - p - 1) * new_step) =
                    new_coefs.row(new_ind + (k - p) * new_step);

            for (int ell = 1; ell <= p; ell++)
            {
                const T alfa = alpha[j][ell - 1];
                const int index = k - p + ell;

                if (math::abs(alfa) == T(0.0))
                    new_coefs.row(new_ind + (index - 1) * new_step) =
                            new_coefs.row(new_ind + index * new_step);
                else
                    new_coefs.row(new_ind + (index - 1) * new_step) =
                            alfa * new_coefs.row(new_ind + (index - 1) * new_step) +
                            (1.0 - alfa) * new_coefs.row(new_ind + index * new_step);
            }
            k--;
        }

        flag = nextCubePoint<gsVector<int> >(new_position, first_point,
                                             last_point);

    } while(nextCubePoint<gsVector<int> >(position, first_point,
                                          last_point));

    coefs = give(new_coefs);

    if (update_knots)
        knots = KnotVectorType(p, nknots.begin(), nknots.end());
}


template <short_t d, typename KnotVectorType, typename Mat, typename ValIt>
void gsTensorBoehmRefineLocal(KnotVectorType& knots,
        const unsigned index,
        Mat& coefs,
        gsVector<index_t, d> &nmb_of_coefs,
        const gsVector<index_t, d> &act_size_of_coefs,
        const gsVector<index_t, d> &size_of_coefs,
        const unsigned direction,
        ValIt valBegin,
        ValIt valEnd,
        const bool update_knots)
{

    if (valBegin==valEnd)
        return;

    typedef typename std::iterator_traits<ValIt>::value_type T;

    const index_t nik = std::distance(valBegin, valEnd); // number of inserted knots
    const index_t p = knots.degree(); // degree
    // number of original (not local) points
    const index_t nopts = knots.size() - p - 1;
    //const int d = size_of_coefs.size();    // dimension


    const index_t a =  knots.iFind(*valBegin)     - knots.begin();
    const index_t b = (knots.iFind(*(valEnd - 1)) - knots.begin()) + 1;

    // allocate a memory for new knots and new control points
    gsSparseVector<T> nknots(b + p + nik);


    gsVector<index_t, d> position(d); // position old points
    position.fill(0);

    gsVector<index_t, d> first_point(position); // first point of a cube
    gsVector<index_t, d> last_point(d); // last point of a cube (old points)
    bspline::getLastIndexLocal<d>(nmb_of_coefs, last_point);
    last_point[direction] = 0;

    // build strides
    gsVector<index_t, d> act_str(d);
    bspline::buildCoeffsStrides<d>(act_size_of_coefs, act_str);

    const index_t step = act_str[direction];

    gsMatrix<T> zero(1, coefs.cols());
    zero.fill(0.0);


    // precompute alpha
    std::vector< std::vector<T> > alpha(nik, std::vector<T> (p));
    computeTensorAlpha<T, KnotVectorType, ValIt, gsSparseVector<T> >
            (alpha, nknots, knots, valBegin, valEnd, true);

    unsigned iterations = 0;
    unsigned number_of_iterations = bspline::numberOfIterations(nmb_of_coefs,
                                                                direction);

    do
    {
        if (number_of_iterations <= iterations)
            break;
        ++iterations;

        ValIt valEndCopy = valEnd;

        index_t ind = bspline::getIndex<d>(act_str, position);

        for (index_t j = b - 1; j < nopts; j++)
        {
            index_t indx1 = j + nik - index;
            index_t indx2 = j - index;

            if (0 <= indx1 &&
                    indx1 < act_size_of_coefs[direction])
            {
                if (indx2 < 0)
                    coefs.row(ind + indx1 * step) = zero.row(0);
                else
                    coefs.row(ind + indx1 * step) =
                            coefs.row(ind + indx2 * step);
            }
        }


        int i = b + p - 1;
        int k = b + p + nik - 1;

        for (int j = nik - 1; 0 <= j; j--)
        {
            const T newKnot = *(--valEndCopy);

            while ((newKnot <= knots[i]) && (a < i))
            {
                index_t indx1 = k - p - 1 - index;
                index_t indx2 = i - p - 1 - index;

                if (indx1 < 0 || act_size_of_coefs[direction] <= indx1)
                {
                    k--;
                    i--;
                    continue;
                }

                if (indx2 < 0)
                {
                    coefs.row(ind + indx1 * step) = zero.row(0);
                }
                else
                {
                    coefs.row(ind + indx1 * step) =
                            coefs.row(ind + indx2 * step);
                }
                k--;
                i--;
            }

            index_t indx1 = k - p - 1 - index;

            if (0 <= indx1 && (indx1 + 1) < act_size_of_coefs[direction])
                coefs.row(ind + indx1 * step) =
                        coefs.row(ind + (indx1 + 1) * step);

            if (indx1 == act_size_of_coefs[direction] - 1)
                coefs.row(ind + indx1 * step) = zero.row(0);

            for (int ell = 1; ell <= p; ell++)
            {
                const T alfa = alpha[j][ell - 1];
                const index_t mindex = k - p + ell - index;

                if (mindex <= 0)
                    continue;

                if (act_size_of_coefs[direction] < mindex)
                    break;

                if (math::abs(alfa) == T(0.0))
                {
                    if (mindex == act_size_of_coefs[direction])
                        coefs.row(ind + (mindex - 1) * step) = zero.row(0);
                    else
                        coefs.row(ind + (mindex - 1) * step) =
                                coefs.row(ind + mindex * step);
                }
                else
                {
                    if (mindex == act_size_of_coefs[direction])
                        coefs.row(ind + (mindex - 1) * step) =
                                alfa * coefs.row(ind + (mindex - 1) * step);
                    else
                        coefs.row(ind + (mindex - 1) * step) =
                                alfa * coefs.row(ind + (mindex - 1) * step) +
                                (1.0 - alfa) * coefs.row(ind + mindex * step);
                }
            }
            k--;
        }
    } while(nextCubePoint<gsVector<index_t, d> >(position, first_point,
                                                  last_point));

    nmb_of_coefs[direction] = size_of_coefs[direction];

    if (update_knots)
    {
        for (ValIt knot_iter = valBegin; knot_iter != valEnd; ++knot_iter)
             knots.insert(*knot_iter, 1);
    }
}

/*
This function performs local Tensor Boehm algorithm. It works.

template <short_t d, typename T, typename KnotVectorType, typename Mat>
void gsTensorBoehmLocal(
        const KnotVectorType& knots,
        unsigned index,
        Mat& coefs,
        const gsVector<unsigned, d>& size_of_coefs,
        const gsVector<unsigned, d>& act_size_of_coefs,
        T val,
        unsigned direction,
        unsigned multiplicity,
        gsVector<unsigned, d>& start,
        gsVector<unsigned, d>& end)
{
    int k = knots.findspan(val);
    int s = knots.multiplicity(val);
    int p = knots.degree();

    const unsigned r = multiplicity;

    std::vector<std::vector<T> > alpha(r, std::vector<T> (p - s));
    computeAlpha<T, KnotVectorType>(alpha, knots, val, r, k, p, s);

    // temporary matrix, for computation
    gsVector<unsigned, d> act_coefs_str(d);
    bspline::buildCoeffsStrides<d>(act_size_of_coefs, act_coefs_str);

    start(direction) = 0;
    end(direction) = 0;
    gsVector<unsigned, d> position(start);


    unsigned step = act_coefs_str[direction];

    do
    {
        const unsigned flat_ind = bspline::getIndex<d>(act_coefs_str, position);


        int low = k - s - index;
        for (int j = size_of_coefs(direction) - 1; low <= j; --j)
            coefs.row(flat_ind + (j + r) * step) =
                    coefs.row(flat_ind + j * step);

        Mat tmp(p + 1, coefs.cols());
        for (int j = 0; j < p + 1; j++)
        {
            tmp.row(j) = coefs.row(flat_ind + (k - p + j - index) * step);
        }

        int L = 0;
        for (unsigned j = 1; j <= r; j++)
        {
            L = k - p + j;

            for (unsigned i = 0; i <= p - j - s; i++)
            {
                T alfa = alpha[j - 1][i];
                tmp.row(i) = alfa * tmp.row(i + 1) + (1.0 - alfa) * tmp.row(i);

//                if (ind < 0)
//                    continue;
//                else if (ind == static_cast<int>(act_size_of_coefs(direction)) - 1)
//                    coefs.row(flat_ind + (ind * step) =
//                            (1 - alfa) * coefs.row(flat_ind + ind * step);
//                else
//                {
//                coefs.row(flat_ind + (ind + 1) * step) =
//                        alfa * coefs.row(flat_ind + (ind + 1) * step) +
//                        (1.0 - alfa) * coefs.row(flat_ind + ind * step);
                //                }
            }

            coefs.row(flat_ind + (L - index) * step) = tmp.row(0);
            coefs.row(flat_ind + (k + r - j - s - index) * step) =
                    tmp.row(p - j - s);

        }
    } while(nextCubePoint<gsVector<unsigned, d> >(position, start, end));
}
*/


template <short_t d, typename T, typename KnotVectorType, typename Mat>
void gsTensorInsertKnotDegreeTimes(
        const KnotVectorType& knots,
        Mat& coefs,
        const gsVector<index_t, d>& size_of_coefs,
        T val,
        unsigned direction,
        gsVector<index_t, d>& start,
        gsVector<index_t, d>& end)
{
    int k = knots.iFind(val) - knots.begin();
    int s = knots.multiplicity(val);
    int p = knots.degree();

    // if multiplicity is greater than degree, we insert no knots
    if (p <= s)
        return;

    const unsigned r = p - s; // how many times we will insert the knot val

    std::vector<std::vector<T> > alpha(r, std::vector<T> (p - s));
    computeAlpha<T, KnotVectorType>(alpha, knots, val, r, k, p, s);

    gsVector<index_t, d> coefs_str(d);
    bspline::buildCoeffsStrides<d>(size_of_coefs, coefs_str);

    start(direction) = 0;
    end(direction) = 0;
    gsVector<index_t, d> position(start);

    unsigned step = coefs_str[direction];

    do
    {
        const unsigned flat_ind = bspline::getIndex<d>(coefs_str, position);
        for (unsigned j = 1; j <= r; j++)
        {
            for (unsigned i = 0; i <= p - j - s; i++)
            {
                T alfa = alpha[j - 1][i];
                coefs.row(flat_ind + i * step) =
                        alfa * coefs.row(flat_ind + (i + 1) * step) +
                        (1.0 - alfa) * coefs.row(flat_ind + i * step);
            }
        }
    } while(nextCubePoint<gsVector<index_t, d> >(position, start, end));
}


// =============================================================================
// Knot removal algorithms (Tiller, NURBS Book Algorithm A5.8)
// =============================================================================

template<class T, class KnotVectorType, class Mat>
bool gsKnotRemoveSingle(
    KnotVectorType & knots,
    Mat            & coefs,
    T                val,
    bool             update_knots)
{
    // Piegl & Tiller, "The NURBS Book" 2nd ed., Algorithm A5.8
    const short_t p = knots.degree();
    const index_t n = coefs.rows() - 1; // last coef index (0-based)

    // r: last-occurrence index of val in the knot vector (0-based)
    const index_t r = knots.iFind(val) - knots.begin();
    const index_t s = knots.multiplicity(val);

    GISMO_ASSERT(s >= 1, "Knot value " << val << " not found in knot vector.");

    const index_t ord   = p + 1;
    const index_t first = r - p;
    const index_t last  = r - s;

    // temp stores candidate control points offset so that
    //   temp[k]  represents the point at absolute index (first - 1 + k),
    //   i.e. temp[0] = P[first-1], temp[last-first+2] = P[last+1].
    const index_t sz = last - first + 3;
    Mat temp(sz, coefs.cols());
    temp.row(0)      = coefs.row(first - 1);
    temp.row(sz - 1) = coefs.row(last  + 1);

    // ii and jj are indices INTO temp (1-based left, (sz-2)-based right).
    // i and j are indices into the original coef array.
    index_t i  = first,   ii = 1;
    index_t j  = last,    jj = sz - 2;

    while (ii <= jj)
    {
        const T alfi = (val - knots[i])   / (knots[i + ord] - knots[i]);
        const T alfj = (val - knots[j])   / (knots[j + ord] - knots[j]);
        temp.row(ii) = (coefs.row(i) - (T(1) - alfi) * temp.row(ii - 1)) / alfi;
        temp.row(jj) = (coefs.row(j) -        alfj   * temp.row(jj + 1)) / (T(1) - alfj);
        ++i; ++ii;
        --j; --jj;
    }

    // P&T feasibility check (The NURBS Book, Algorithm A5.8 / Eq. 5.28):
    // After the loop, if j < i the chains crossed (no middle point): check
    // that temp[ii-1] == temp[jj+1].  If j >= i a middle point remains: check
    // that coefs[i] lies on the interpolated line between temp[ii-1] and
    // temp[jj+1] (reference: NURBS-Python helpers.py knot_removal, line 698-712).
    if (j < i)  // chains crossed — odd number of points in [first,last]
    {
        if ((temp.row(ii - 1) - temp.row(jj + 1)).norm() > T(1e-10))
            return false;
    }
    else  // middle point remains — check interpolation
    {
        const T alfa = (val - knots[i]) / (knots[i + ord] - knots[i]);
        if ((coefs.row(i) - (alfa * temp.row(ii - 1) + (T(1) - alfa) * temp.row(jj + 1))).norm() > T(1e-10))
            return false;
    }

    // Write back new control points.
    // Mirror the computation: left chain → coefs[first..], right chain → ..coefs[last].
    // Middle point (if sz is odd) is left unchanged.
    {
        index_t wi = first, wj = last, wii = 1, wjj = sz - 2;
        while (wi < wj)
        {
            coefs.row(wi) = temp.row(wii);
            coefs.row(wj) = temp.row(wjj);
            ++wi; ++wii;
            --wj; --wjj;
        }
    }

    // Shift [last+1 .. n] left by one to close the gap
    for (index_t idx = last + 1; idx <= n; ++idx)
        coefs.row(idx - 1) = coefs.row(idx);

    coefs.conservativeResize(n, coefs.cols());

    if (update_knots)
        knots.remove(val, 1);

    return true;
}


template<class T, class KnotVectorType, class Mat>
index_t gsKnotRemove(
    KnotVectorType & knots,
    Mat            & coefs,
    T                val,
    index_t              t,
    bool             update_knots)
{
    index_t removed = 0;
    for (index_t i = 0; i < t; ++i)
    {
        if (!gsKnotRemoveSingle<T>(knots, coefs, val, update_knots))
            break;
        ++removed;
    }
    return removed;
}


template<typename T, typename KnotVectorType, typename Mat>
index_t gsTensorKnotRemove(
        KnotVectorType        & knots,
        Mat                   & coefs,
        T                       val,
        index_t                 direction,
        gsVector<index_t>       str,
        index_t                 t,
        bool                    update_knots)
{
    GISMO_ASSERT(t >= 1, "Must remove at least once.");
    GISMO_ASSERT(direction < str.size(),
                 "Direction out of range.");

    const short_t d           = str.size();
    const index_t num_in_dir  = knots.size() - knots.degree() - 1;
    const index_t num_fibers  = coefs.rows() / num_in_dir;
    const index_t step        = str[direction];

    // Helper: extract a single fiber (column of coefficients along direction)
    // and write it back.
    index_t total_removed = t; // will be reduced to minimum over fibers per pass

    for (index_t pass = 0; pass < t; ++pass)
    {
        // Try one removal on all fibers; if any fails the pass is aborted
        // and coefs/knots are restored from saved copies.
        Mat coefs_backup = coefs;

        // Re-compute strides from current state
        gsVector<index_t> position(d);
        position.fill(0);
        gsVector<index_t> first_point(d);
        first_point.fill(0);
        gsVector<index_t> last_point(d);
        getLastIndex(str, coefs.rows(), last_point);
        last_point[direction] = 0;

        bool all_ok = true;

        do
        {
            const index_t base = getIndex(str, position);

            // Build fiber coefficient matrix (num_in_dir rows)
            Mat fiber(num_in_dir, coefs.cols());
            for (index_t i = 0; i < num_in_dir; ++i)
                fiber.row(i) = coefs.row(base + i * step);

            // Attempt removal on fiber (don't update knots yet)
            KnotVectorType knots_copy = knots;
            if (!gsKnotRemoveSingle<T>(knots_copy, fiber, val, false))
            {
                all_ok = false;
                break;
            }

            // Write fiber back (now has num_in_dir - 1 rows)
            // We write into coefs_backup at the right positions
            // Note: coefs_backup still has the old size; we will
            // rebuild coefs after all fibers succeed.
            for (index_t i = 0; i < num_in_dir - 1; ++i)
                coefs_backup.row(base + i * step) = fiber.row(i);
            // The "last" slot is temporarily stale — will be fixed
            // when we compact below.

        } while (nextCubePoint<gsVector<index_t>>(position, first_point, last_point));

        if (!all_ok)
        {
            // Restore and stop
            total_removed = pass;
            break;
        }

        // All fibers succeeded: compact the coefficient matrix
        // (remove one "layer" per direction stride)
        const index_t new_num_in_dir = num_in_dir - 1;
        const index_t new_npts = num_fibers * new_num_in_dir;
        Mat new_coefs(new_npts, coefs.cols());

        // Recompute stride for the compacted matrix
        gsVector<index_t> new_str(str);
        correctNewStride(new_str, str, direction, -1);  // -1 removal

        gsVector<index_t> pos(d);  pos.fill(0);
        gsVector<index_t> fp(d);   fp.fill(0);
        gsVector<index_t> lp(d);
        getLastIndex(str, coefs.rows(), lp);
        lp[direction] = 0;
        gsVector<index_t> new_pos(d); new_pos.fill(0);
        gsVector<index_t> new_lp(d);
        getLastIndex(new_str, new_npts, new_lp);
        new_lp[direction] = 0;

        const index_t new_step = new_str[direction];

        bool flag2 = true;
        do
        {
            if (!flag2)
                GISMO_ERROR("gsTensorKnotRemove: index error during compaction.");

            const index_t ind     = getIndex(str,     pos);
            const index_t new_ind = getIndex(new_str, new_pos);

            // Rebuild fiber from coefs_backup
            Mat fiber(num_in_dir, coefs.cols());
            for (index_t i = 0; i < num_in_dir; ++i)
                fiber.row(i) = coefs.row(ind + i * step);

            // Re-run the (guaranteed) removal to get the compacted fiber
            KnotVectorType knots_dummy = knots;
            gsKnotRemoveSingle<T>(knots_dummy, fiber, val, false);

            for (index_t i = 0; i < new_num_in_dir; ++i)
                new_coefs.row(new_ind + i * new_step) = fiber.row(i);

            flag2 = nextCubePoint<gsVector<index_t>>(new_pos, fp, new_lp);

        } while (nextCubePoint<gsVector<index_t>>(pos, fp, lp));

        coefs = give(new_coefs);
        str   = new_str;

        if (update_knots)
            knots.remove(val, 1);
    }

    if (total_removed == t) // loop completed without break
        return t;
    return total_removed;
}

} // namespace gismo
