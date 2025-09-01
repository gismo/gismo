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

    template <class T, size_t _order>
    class NullObject : public BaseObject<T, _order, true, 0>
    {
        using Base = BaseObject<T, _order, true, 0>;
    public:
        typedef T Scalar; // Define Scalar type for this expression
        static constexpr size_t order = Base::order;
        static constexpr size_t space = Base::space;
        static constexpr size_t deriv = Base::deriv;
        static constexpr bool isConstant = Base::isConstant;

    public:

        gsMatrix<Scalar> eval(const index_t) const
        {
            GISMO_ERROR("The null expression should not be evaluated");
        }

        void print(std::ostream &) const
        {
            GISMO_NO_IMPLEMENTATION;
        }

        operator const SpaceObject<T,Space::None,_order> & () const
        {
            static SpaceObject<T,Space::None,_order> vv(0,{}, 0);
            return vv;
        }

        explicit NullObject()
        :
        BaseObject<T, _order, true, 0>(0,{})
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