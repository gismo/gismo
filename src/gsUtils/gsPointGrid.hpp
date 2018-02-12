/** @file gsPointGrid.hpp

    @brief Provides implementation of functions which generate
    structured point data

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): C. Hofreither, A. Mantzaflaris
*/

#pragma once


namespace gismo 
{

/* **************** Utility functions **************** */


template<typename T> //cwiseSampleCount
gsVector<unsigned> uniformSampleCount (const gsVector<T>& lower, 
                                       const gsVector<T>& upper, 
                                       int numPoints)
{
    const index_t d = lower.rows();
    assert( d == upper.rows() );

    // TO do : phys. volume criterion
    gsVector<T> span = upper - lower;
    const T volume = span.prod();
    const T h = math::pow(volume / numPoints, 1.0 / d);

    gsVector<unsigned> np(d);

    for (index_t i = 0; i < d; ++i)
    {
        np[i] = cast<T,unsigned>(math::ceil( span[i] / h ) );
        GISMO_ASSERT( np[i] > 0, "Something went wrong, number of points is zero..");
    }

    return np;
}

template<typename T>
void uniformIntervals(const gsVector<T>& lower, 
                      const gsVector<T>& upper, 
                      std::vector< std::vector<T> >& intervals, 
                      int numIntervals)
{
    const int d = lower.rows();
    assert( d == upper.rows() );

    gsVector<unsigned> np = uniformSampleCount( lower, upper, numIntervals );

    // resize to dimension d without copying old contents, if any
    intervals.clear();
    intervals.resize(d);

    for (int i = 0; i < d; ++i)
    {
        int numInt = np[i] - 1;
        if (numInt <= 1)
            numInt = 1;

        const T h = T(1) / numInt;

        intervals[i].resize(numInt + 1);
        for (int j = 0; j <= numInt; ++j)
            intervals[i][j] = j * h;
    }
}


/* **************** Uniform grids described by limits/corners **************** */


template<class T>
gsMatrix<T> gsPointGrid( gsVector<T> const & a, 
                         gsVector<T> const & b, 
                         gsVector<unsigned> const & np )
{
    gsMatrix<T> res(a.size(), np.prod() );
    gsGridIterator<T,CUBE> pt(a, b, np.cast<index_t>());
    for(index_t c = 0; pt; ++pt, ++c)
        res.col(c) = *pt;
    return res;
}

template<typename T>
gsMatrix<T> uniformPointGrid(const gsVector<T>& lower, 
                             const gsVector<T>& upper, 
                             int numPoints)
{
    const gsVector<unsigned> cwisePoints = uniformSampleCount(lower, upper, numPoints);
    return gsPointGrid(lower, upper, cwisePoints); // note: structure lost
}

};// namespace gismo
