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
    // Traits class for BaseObject to avoid circular dependency
    template <class T, size_t _order, bool _isConstant, size_t _space>
    struct ExpressionTraits<BaseObject<T, _order, _isConstant, _space>>
    {
        typedef T Scalar;
        static constexpr size_t order = _order;  // Use template parameter, not BaseObject::order
        static constexpr size_t space = _space;
        static constexpr size_t deriv = 0;
        static constexpr bool isConstant = _isConstant;
    };

    // IsConstant: Flag that indicates if the expression is constant, e.g., its derivatives are zero
    // Space: Flag that indicates whether the expression is a space
    template <class T, size_t _order, bool _isConstant = true, size_t _space = Space::None>
    class BaseObject : public BaseExpression<BaseObject<T,_order,_isConstant,_space>>
    {
        using Base = BaseExpression<BaseObject<T,_order,_isConstant,_space>>;

    protected:
        const std::array<size_t, _order> sizes_;
        const size_t domainDim_;
        const std::string label_;

    public:
        typedef typename ExpressionTraits<BaseObject<T,_order,_isConstant,_space>>::Scalar Scalar;
        static constexpr size_t order = ExpressionTraits<BaseObject<T,_order,_isConstant,_space>>::order;
        static constexpr size_t space = ExpressionTraits<BaseObject<T,_order,_isConstant,_space>>::space;
        static constexpr size_t deriv = ExpressionTraits<BaseObject<T,_order,_isConstant,_space>>::deriv;
        static constexpr bool isConstant = ExpressionTraits<BaseObject<T,_order,_isConstant,_space>>::isConstant;

        // Access to sizes is direct
        const std::array<size_t, _order>& sizes() const override { return sizes_; }

        size_t domainDim() const { return domainDim_; }


    public:

        gsMatrix<Scalar> eval(const index_t) const
        {
            GISMO_NO_IMPLEMENTATION;
        }

        void parse(gismo::ExpressionHelper<Scalar> &) const
        {
            GISMO_NO_IMPLEMENTATION;
        }

        virtual const SpaceObject<Scalar,Space::None,_order> & rowVar() const {return NullObject<T,_order>::get();}
        virtual const SpaceObject<Scalar,Space::None,_order> & colVar() const {return NullObject<T,_order>::get();}

    protected:

        explicit BaseObject(size_t domainDim, const std::array<size_t, order> & input_sizes, std::string label = "a")
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
        ~BaseObject() = default;

        // // Helper to calculate total elements from sizes (used internally by derived classes)
        // static size_t tensorSize(std::array<size_t, order> dims)
        // {
        //     if (order == 0) return 1; // Scalar has 1 element
        //     size_t total = 1;
        //     for (size_t dim_size : dims)
        //     {
        //         total *= dim_size;
        //     }
        //     return total;
        // }

    };
}//namespace Expr
}//namespace gismo