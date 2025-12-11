/** @file BaseObject.h
 *
    @brief

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M.Verhelst
*/

#pragma once

#include <gsMatrix/gsMatrix.h>
#include <gsNewExpressions/ExpressionUtils.h>

namespace gismo
{
namespace Expr
{
template <typename E>
class BaseObject : public BaseExpression<E>
{
    using Base = BaseExpression<E>;

protected:
    typedef typename Base::Scalar T;

    const std::array<size_t, Base::Order> sizes_;
    const size_t domainDim_;
    const std::string label_;

public:
    // Expose the static traits publicly so they can be accessed as LhsExpr::Order etc.
    using Base::Order;
    using Base::Space;
    using Base::Deriv;
    using Base::IsConstant;
    using Base::Scalar;

    // Access to sizes is direct
    const std::array<size_t, Order>& sizes() const  { return sizes_; }

    size_t domainDim() const { return domainDim_; }

public:

    ExpressionValue<T> eval(const index_t) const
    {
        GISMO_NO_IMPLEMENTATION;
    }

    void parse(gismo::ExpressionHelper<T> &) const
    {
        GISMO_NO_IMPLEMENTATION;
    }

    ~BaseObject() = default;

    explicit BaseObject(size_t domainDim, const std::array<size_t, Order> & input_sizes, std::string label = "a")
    :
    Base(), // Initialize base class properly
    sizes_(input_sizes),
    domainDim_(domainDim),
    label_(label)
    {}

    BaseObject(const BaseObject& other)
    :
    Base(),
    sizes_(other.sizes_),
    domainDim_(other.domainDim_),
    label_(other.label_)
    {}

    BaseObject& operator=(const BaseObject&)
    {
        return *this;
    }

    std::string label() const { return label_; }

protected:



    // // Helper to calculate total elements from sizes (used internally by derived classes)
    // static size_t tensorSize(std::array<size_t, Order> dims)
    // {
    //     if (Order == 0) return 1; // Scalar has 1 element
    //     size_t total = 1;
    //     for (size_t dim_size : dims)
    //     {
    //         total *= dim_size;
    //     }
    //     return total;
    // }

};

template <typename E, typename SpaceObj>
VariationObject<BaseObject<E>, SpaceObj> variation(const BaseObject<E> & sol,const SpaceObj & space)
{
    return VariationObject<BaseObject<E>, SpaceObj>(sol,space, false);
}

}//namespace Expr
}//namespace gismo