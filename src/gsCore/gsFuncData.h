/** @file gsFuncData.h

    @brief
    This object is a cache for computed values from an evaluator.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Bressan, A. Manzaflaris
*/

#pragma once

#include<gsCore/gsLinearAlgebra.h>
#include<gsCore/gsBoundary.h>

namespace gismo
{

template <typename T> class gsFunctionSet;

/**
   @brief the gsFuncData is a cache of pre-computed function sets values.

   Which data is contained in the gsFuncData is specidied by a flag system.
   The user must set the \a flags memeber using a combination of the constants
   \a NEED_VALUE, \a NEED_DERIV, ... then the cache is filled by calling the
   \a gsFunctionSet::compute(points, gsFuncData&) where points can either be
   a gsMatrix<> containing the point coordinates or a gsMapData object.

   The row matrix data are public members. There are also accessor functions
   that provide a per-point view of the data in a different format: each column
   corresponds to a different function object.

   \tparam T numeric type
*/
template <typename T>
class gsFuncData
{
    friend class gsFunctionSet<T>;

public:
    typedef typename gsMatrix<T>::constColumn constColumn;
    // Types for returning quick access to data in matrix format
    typedef gsAsConstMatrix<T, -1, -1>                  matrixView;
    typedef Eigen::Transpose<typename matrixView::Base> matrixTransposeView;

public:
    mutable unsigned flags;
    index_t      patchId; // move to mapdata

    gsMatrix<index_t> actives;

    /// Stores values and derivatives
    /// values[0] for base
    /// values[n] for n-th derivative
    std::vector<gsMatrix<T> > values;

    gsMatrix<T> curls;
    gsMatrix<T> divs;
    gsMatrix<T> laplacians;

    /// \brief Dimension of the (source) domain and the target (image) space.
    /// dim.first refers to ParDim, dim.second refers to GeoDim
    /// @return For \f$f:\mathbb{R}^n\rightarrow\mathbb{R}^m\f$ returns \f$n\f$.
    std::pair<short_t, short_t> dim;

public:
    /**
     * @brief Main constructor
     * @param flags what to compute
     * @param patch in case of multipatch structures, on which patch to compute
     */
    explicit gsFuncData(unsigned flags = 0, int patch = 0)
    : flags(flags), patchId(patch)
    { }

public:

    /**
     * @brief addFlags
     * set the evaluator to compute additional values
     * @param newFlags
     */
    void addFlags (unsigned newFlags)
    { flags = flags|newFlags; }


    int maxDeriv() const
    {
        if (flags & (NEED_LAPLACIAN|NEED_DERIV2) )
            return 2;
        else if (flags & (NEED_DERIV|NEED_CURL|NEED_DIV) )
            return 1;
        else if (flags & (NEED_VALUE) )
            return 0;
        return -1;
    }

    /**
     * @brief Provides memory usage information
     * @return the number of bytes occupied by this object
     */
    unsigned bytesUsed() const
    { /*to do*/
        return 1;
    }

    /// \brief Clear the memory that this object uses
    void clear() { /*to do*/}

    /// \brief Swaps this object with \a other
    void swap(gsFuncData & other)
    {
        std::swap(flags  , other.flags  );
        std::swap(patchId, other.patchId);
        std::swap(dim, other.dim);
        actives   .swap(other.actives   );
        values    .swap(other.values    );
        curls     .swap(other.curls     );
        divs      .swap(other.divs      );
        laplacians.swap(other.laplacians);
    }

public:

    inline const gsMatrix<index_t> & allActives() const
    {
        GISMO_ASSERT(flags & NEED_ACTIVE,
                   "actives are not computed unless the NEED_ACTIVE flag is set.");
        GISMO_ASSERT(0!=actives.size(), "actives were not computed.");
        return actives;
    }

    inline const gsMatrix<T> & allValues() const
    {
        GISMO_ASSERT(flags & NEED_ACTIVE,
                     "values are not computed unless the NEED_ACTIVE flag is set.");
        GISMO_ASSERT(0!=values.size(), "values were not computed.");
        return values.front();
    }

    inline gsMatrix<index_t>::constColumn active(index_t point = 0) const
    {
        GISMO_ASSERT(flags & NEED_ACTIVE,
                   "actives are not computed unless the NEED_ACTIVE flag is set.");
        return actives.col(point);
    }

    inline matrixView eval (index_t point) const
    {
        GISMO_ASSERT(flags & NEED_VALUE,
                   "values are not computed unless the NEED_VALUE flag is set.");
        return values[0].reshapeCol(point, dim.second, values[0].rows()/dim.second);
    }

    inline matrixView deriv (index_t point) const
    {
        GISMO_ASSERT(flags & NEED_DERIV,
                   "derivs are not computed unless the NEED_DERIV flag is set.");
        return values[1].reshapeCol(point, derivSize(), values[1].rows()/derivSize());
    }

    inline matrixView deriv2 (index_t point) const
    {
        GISMO_ASSERT(flags & NEED_DERIV2,
                   "deriv2s are not computed unless the NEED_DERIV2 flag is set.");
        return values[2].reshapeCol(point, deriv2Size(), values[2].rows()/deriv2Size());
    }

    inline matrixView curl (index_t point) const
    {
        GISMO_ASSERT(flags & NEED_CURL,
                   "curls are not computed unless the NEED_CURL flag is set.");
        return curls.reshapeCol(point, dim.second, curls.rows()/dim.second );
    }

    inline matrixView div (index_t point) const
    {
        GISMO_ASSERT(flags & NEED_DIV,
                   "divs are not computed unless the NEED_DIV flag is set.");
        return divs.reshapeCol(point, divSize(), divs.rows()/divSize() );
    }

    inline matrixView laplacian (index_t point) const
    {
        GISMO_ASSERT(flags & NEED_LAPLACIAN,
                   "laplacians are not computed unless the NEED_LAPLACIAN flag is set.");
        return laplacians.reshapeCol(point, dim.second, laplacians.rows()/dim.second );
    }


    inline matrixTransposeView jacobian (index_t point, index_t func = 0) const
    {
       GISMO_ASSERT(flags & NEED_DERIV,
                  "jacobian access needs the computation of derivs: set the NEED_DERIV flag.");
       return gsAsConstMatrix<T, Dynamic, Dynamic>(&values[1].coeffRef(func*derivSize(),point), dim.first,dim.second).transpose();
    }

//protected:

    /// Number of partial derivatives (<em>dim.second*dim.first</em>).
    int  derivSize () const {return dim.first*dim.second;}
    /// Number of 2nd derivatives (<em>dim.second*dim.first*(dim.first+1)/2</em>).
    int  deriv2Size() const {return dim.second*dim.first*(dim.first+1) / 2; }
    /// Size of computed divergence (<em>dim.second/dim.first</em>).
    int  divSize   () const {return dim.second/dim.first;}
};


/**
 @brief the gsMapData is a cache of pre-computed function (map) values.

 \tparam T numeric type
 \sa gsFuncData
 */
template <typename T>
class gsMapData : public gsFuncData<T>
{
public:
    typedef gsFuncData<T> Base;
    typedef typename Base::constColumn         constColumn;
    typedef typename Base::matrixView          matrixView;
    typedef typename Base::matrixTransposeView matrixTransposeView;

public:
    /**
     * @brief Main constructor
     * @param flags what to compute
     */
    explicit gsMapData(unsigned flags = 0)
    : Base(flags), side(boundary::none)
    { }

public:
    using Base::flags;
    using Base::values;
    using Base::dim;

    boxSide     side;

    gsMatrix<T> points;     ///< input (parametric) points

    gsMatrix<T> measures;
    gsMatrix<T> fundForms;  // First fundumental forms
    gsMatrix<T> normals;
    gsMatrix<T> outNormals; // only for the boundary

public:
    inline constColumn point(const index_t point) const { return points.col(point);}

    inline T measure(const index_t point) const
    {
        GISMO_ASSERT(flags & NEED_MEASURE,
                   "measures are not computed unless the NEED_MEASURE flag is set.");
        return measures(0,point);
    }

    inline matrixView fundForm(const index_t point) const
    {
        GISMO_ASSERT(flags & NEED_GRAD_TRANSFORM,
                   "fundForms are not computed unless the NEED_GRAD_TRANSFORM flag is set.");
        return fundForms.reshapeCol(point, dim.second, dim.first);
    }

    inline constColumn normal(const index_t point) const
    {
        GISMO_ASSERT(flags & NEED_NORMAL,
                   "normals are not computed unless the NEED_NORMAL flag is set.");
        return normals.col(point);
    }

    inline constColumn outNormal(const index_t point) const
    {
        GISMO_ASSERT(flags & NEED_OUTER_NORMAL,
                   "normals are not computed unless the NEED_NORMAL flag is set.");
        return outNormals.col(point);
    }

    inline matrixTransposeView jacobians () const
    {
       GISMO_ASSERT(flags & NEED_DERIV,
                  "jacobian access needs the computation of derivs: set the NEED_DERIV flag.");
       return gsAsConstMatrix<T, Dynamic, Dynamic>(&values[1].coeffRef(0,0), dim.first,dim.second*values[1].cols()).transpose();
    }
};

} // namespace gismo

namespace std
{

// Specializations of std::swap
template <typename T>
void swap(gismo::gsFuncData<T> & f1, gismo::gsFuncData<T> & f2)
{
    f1.swap(f2);
}

} // namespace std
