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

    GISMO_ASSERT( coefs.rows() == index_t(knots.size() - knots.degree()-1),
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

    GISMO_ASSERT( coefs.rows() == index_t(knots.size() - knots.degree()-1),
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
        coefs.row(i) = (1-a) * coefs.row(i-1) + a * coefs.row(i);
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
        coefs.row(i) = (1-a) * coefs.row(i-1) + a * coefs.row(i);
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

            if( math::abs(alfa) == 0.0 )
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

                if (math::abs(alfa) == 0.0)
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

                if (math::abs(alfa) == 0.0)
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

} // namespace gismo
