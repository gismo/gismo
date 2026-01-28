/** @file NullObject.h

    @brief

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M.Verhelst
*/

#pragma once

#include <gsNewExpressions/ExpressionTraits.h>
#include <gsNewExpressions/SpaceObject.h>

namespace gismo
{
namespace Expr
{
    // Need to define ExpressionTraits for NullObject
    template<class T, enum SpaceType _Space, size_t _Order>
    struct ExpressionTraits<NullObject<T, _Space, _Order>>
    {
        typedef T Scalar;
        static constexpr size_t Order = _Order;
        static constexpr SpaceType Space = _Space;
        static constexpr size_t Deriv = 0;
        static constexpr bool IsConstant = true;
    };

    template <class T, enum SpaceType _Space, size_t _Order>
    class NullObject : public SpaceObject<T, _Space, _Order>
    {
        using Base = SpaceObject<T, _Space, _Order>;

    public:
        // Expose the static traits publicly
        using Base::Order;
        using Base::Space;
        using Base::Deriv;
        using Base::IsConstant;
        using Base::sizes_;
        typedef typename Base::Scalar Scalar;

    public:

        operator const SpaceObject<T,SpaceType::None,_Order> & () const
        {
            static SpaceObject<T,SpaceType::None,_Order> vv(0,{}, 0);
            return vv;
        }

        ExpressionResult<Scalar> eval(const index_t) const
        {
            GISMO_ERROR("The null expression should not be evaluated");
        }

        void print(std::ostream &) const
        {
            GISMO_NO_IMPLEMENTATION;
        }

        explicit NullObject()
        :
        Base(0,{}, -1)
        {}

        NullObject(const NullObject&) = default;
        NullObject& operator=(const NullObject&) = default;
        ~NullObject() = default;

        static const NullObject & get()
        {
            static NullObject o;
            return o;
        }
    };
}//namespace Expr
}//namespace gismo