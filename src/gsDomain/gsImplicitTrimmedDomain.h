/** @file gsTrimmedDomain.h

    @brief TODO

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst, A. Mantzaflaris
*/

#pragma once

#include <gsDomain/gsTrimmedDomain.h>

#include <type_traits>


namespace gismo
{

/**
 * @brief TODO
 *
 * \ingroup Domain
 */

/// Class representing an implicit trimmed domain
template<short_t d, class T, class Z=size_t>
class gsImplicitTrimmedDomain : public gsTrimmedDomain<d,T,Z>
{
    static bool constexpr implicit = true;
    typedef gsTrimmedDomain<d,T,Z> Base;
private:
    typename gsFunction<T>::Ptr m_implFunction; // implicit function

public:
    gsImplicitTrimmedDomain(const gsFunction<T> & fnc,
                        const gsTensorBSplineBasis<d,T> & tbasis,
                        index_t samples = 5) :
    m_implFunction(memory::make_shared_not_owned(&fnc))
    {
        this->init(tbasis, samples);
    }

    /// Constructor for hierarchical tensor bases
    gsImplicitTrimmedDomain(const gsFunction<T> & fnc,
                        const gsHTensorBasis<d,T> & htbasis,
                        index_t samples = 5) :
    m_implFunction(memory::make_shared_not_owned(&fnc))
    {
        this->init(htbasis, samples);
    }

    /// Constructor for bounding-box + uniform grid.
    /// No basis is involved here, so the degree used for quadrature sizing
    /// (Base::degree()) is taken from \a deg.
    gsImplicitTrimmedDomain(const gsFunction<T>       & fnc,
                        const gsMatrix<T>            & bbox,
                        const gsVector<index_t,d>    & numCells,
                        index_t samples = 5,
                        short_t deg     = 1) :
    m_implFunction(memory::make_shared_not_owned(&fnc))
    {
        this->init(bbox, numCells, samples, deg);
    }

    /// Constructor for the size-based adaptive grid (no basis, no bounding box given).
    /// There is no basis on this path, so the degree used for quadrature sizing
    /// (Base::degree()) is taken from \a deg.
    /// \note This requires a level set whose support() returns a finite \c d x 2
    /// bounding box. gsFunctionExpr and most analytic level sets have unbounded
    /// support — use gsImplicitTrimmedDomain(fnc, bbox, numCells, samples, deg)
    /// instead which supplies the box explicitly.
    gsImplicitTrimmedDomain(const gsFunction<T> & fnc, short_t deg = 1)
    : m_implFunction(memory::make_shared_not_owned(&fnc))
    {
        this->init(T(1), T(0.1), 10, deg);
    }

    /// A floating-point degree is always a mistake: gsImplicitTrimmedDomain(phi, 2.9)
    /// would otherwise compile and silently truncate to degree 2 (parenthesised
    /// direct-initialisation does not diagnose narrowing, and this build does not
    /// enable -Wconversion). Deleting it turns the mistake into a compile error.
    /// NOTE: this MUST stay a constrained template. A plain non-template
    /// `gsImplicitTrimmedDomain(const gsFunction<T>&, double) = delete;` would make the
    /// legitimate call `d(phi, 3)` AMBIGUOUS: int->double and int->short_t are both
    /// standard conversions of the same rank.
    template<class D,
             typename std::enable_if<std::is_floating_point<D>::value, int>::type = 0>
    gsImplicitTrimmedDomain(const gsFunction<T> & fnc, D deg) = delete;

    const gsFunction<T> & implicitFunction() const { return *m_implFunction; };

    inline gsVector<short_t> sign(const gsMatrix<T> & u)
    {
        gsVector<T> val = m_implFunction->eval(u).row(0);
        return gsVector<short_t>(val.array().sign().template cast<short_t>());
    }

    gsMatrix<T> boundingBox() const override
    {
        //always exists?
        return m_implFunction->support(); 
    }

};

} // namespace gismo
