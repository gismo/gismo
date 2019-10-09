/** @file gsFunctionSet.hpp

    @brief implementation of default functions of the gsFunctionSet

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Bressan
**/

#include <gsCore/gsFuncData.h>
#include <gsCore/gsFunction.h>
#include <gsCore/gsBasis.h>

namespace gismo
{

template <class T>
gsFunctionSet<T>::gsFunctionSet() {}

template <class T>
gsFunctionSet<T>::gsFunctionSet(const gsFunctionSet &) {}

template <class T>
gsFunctionSet<T>::~gsFunctionSet () {}

template <class T>
const gsFunction<T>& gsFunctionSet<T>::function(const index_t k) const
{
    GISMO_ASSERT(dynamic_cast<const gsFunction<T>*>(&piece(k)),
                 "No function found, instead: "<< piece(k));
    return static_cast<const gsFunction<T>&>(piece(k));
}

template <class T>
const gsBasis<T>& gsFunctionSet<T>::basis(const index_t k) const
{
    GISMO_ASSERT(dynamic_cast<const gsBasis<T>*>(&piece(k)),
                 "No basis found, instead: "<< piece(k));
    return static_cast<const gsBasis<T>&>(piece(k));
}

// support (domain of definition)
template <class T>
gsMatrix<T> gsFunctionSet<T>::support() const
{
    return gsMatrix<T>();
}

// actives

template <typename T>
void gsFunctionSet<T>::active_into (const gsMatrix<T> &, gsMatrix<index_t> &) const
{
    GISMO_NO_IMPLEMENTATION
    // Single function 0 globally active:
    // result.setConstant(1,u.cols(),0);
}

// evaluation

template <typename T>
void gsFunctionSet<T>::eval_into (const gsMatrix<T> &, gsMatrix<T> &) const
{GISMO_NO_IMPLEMENTATION}

template <typename T>
void gsFunctionSet<T>::deriv_into (const gsMatrix<T> &, gsMatrix<T> &) const
{GISMO_NO_IMPLEMENTATION}

template <typename T>
void gsFunctionSet<T>::deriv2_into (const gsMatrix<T> &, gsMatrix<T> &) const
{GISMO_NO_IMPLEMENTATION}

template <typename T>
void gsFunctionSet<T>::evalAllDers_into(const gsMatrix<T> & u, const int n,
                                        std::vector<gsMatrix<T> > & result) const
{
    result.resize(n+1);

    switch(n)
    {
    case 0:
        eval_into(u, result[0]);
        break;
    case 1:
        eval_into (u, result[0]);
        deriv_into(u, result[1]);
        break;
    case 2:
        eval_into  (u, result[0]);
        deriv_into (u, result[1]);
        deriv2_into(u, result[2]);
        break;
    default:
        GISMO_ERROR("evalAllDers implemented for order up to 2<"<<n ); //<< " for "<<*this);
        break;
    }
}

template <class T>
gsMatrix<T>
gsFunctionSet<T>::eval(const gsMatrix<T>& u) const
{
    gsMatrix<T> result;
    this->eval_into( u, result );
    return result;
}

template <class T>
gsMatrix<T>
gsFunctionSet<T>::deriv(const gsMatrix<T>& u) const
{
    gsMatrix<T> result;
    this->deriv_into( u, result );
    return result;
}

template <class T>
gsMatrix<T>
gsFunctionSet<T>::deriv2(const gsMatrix<T>& u) const
{
    gsMatrix<T> result;
    this->deriv2_into( u, result );
    return result;
}

/*
template <typename T>
void gsFunctionSet<T>::div_into       (const gsMatrix<T> & u, gsMatrix<T> &result) const
{
    gsMatrix<T> tmp;
    deriv_into(u,tmp);
    convertValue<T>::derivToDiv(tmp, result, info());
}

template <typename T>
void gsFunctionSet<T>::curl_into      (const gsMatrix<T> & u, gsMatrix<T> &result) const
{
    gsMatrix<T> tmp;
    deriv_into(u,tmp);
    convertValue<T>::derivToCurl(tmp, result, info());
}

template <typename T>
void gsFunctionSet<T>::laplacian_into (const gsMatrix<T> & u, gsMatrix<T> &result) const
{
    gsMatrix<T> tmp;
    deriv2_into(u,tmp);
    convertValue<T>::deriv2ToLaplacian(tmp, result, info());
}
*/


// Returns quantities either on the target domain or on the parametric
// domain depending on the representation of the object
template <typename T>
void gsFunctionSet<T>::compute(const gsMatrix<T> & in,
                               gsFuncData<T> & out   ) const
{
    const unsigned flags = out.flags;

    out.dim = this->dimensions();

    const int md = out.maxDeriv();
    if (md != -1)
        evalAllDers_into(in, md, out.values);

    if (flags & NEED_ACTIVE && flags & SAME_ELEMENT)
    {
        GISMO_ASSERT(0!=in.cols(), "The points are empty.");
        active_into(in.col(0), out.actives);
    }
    else if (flags & NEED_ACTIVE)
        active_into(in, out.actives);

    // if ( flags & NEED_DIV )
    //     convertValue<T>::derivToDiv(out.values[1], out.divs, info());
    // if ( flags & NEED_CURL )
    //     convertValue<T>::derivToCurl(out.values[1], out.curls, info());
    if (flags & NEED_LAPLACIAN)
    {
        const index_t dsz    = out.deriv2Size();
        const index_t numact = out.values[2].rows() / dsz;
        out.laplacians.resize(numact, in.cols());
        for (index_t i=0; i!= numact; ++i)
            out.laplacians.row(i) =
                out.values[2].middleRows(dsz*i,out.dim.first).colwise().sum();
    }
}

} // namespace gismo
