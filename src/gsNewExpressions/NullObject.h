/** @file NullObject.h

    @brief

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M.Verhelst
*/

#pragma once

#include <gsNewExpressions/SpaceObject.h>

namespace gismo
{
namespace Expr
{

    template<class T, size_t _Space, size_t _Order>
    struct ExpressionTraits<NullObject<T, _Space, _Order>>
    {
        typedef T Scalar;
        static constexpr size_t Order = _Order;
        static constexpr size_t Space = _Space;
        static constexpr size_t Deriv = 0;
        static constexpr bool IsConstant = true;
    };

    template <class T, size_t _Space, size_t _Order>
    class NullObject : public BaseObject<NullObject<T, _Space, _Order>>
    {
        using Base = BaseObject<NullObject<T, _Space, _Order>>;

    public:
        // Expose the static traits publicly
        using Base::Order;
        using Base::Space;
        using Base::Deriv;
        using Base::IsConstant;
        using Base::sizes_;
        typedef typename Base::Scalar Scalar;

    public:

        ExpressionValue<Scalar> eval(const index_t) const
        {
            GISMO_ERROR("The null expression should not be evaluated");
        }

        void print(std::ostream &) const
        {
            GISMO_NO_IMPLEMENTATION;
        }

        // operator const NullSpaceObject<T,_Space,_Order> & () const
        // {
        //     static NullSpaceObject<T,_Space,_Order> vv(0,{}, 0);
        //     return vv;
        // }

        explicit NullObject()
        :
        Base(0,{})
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