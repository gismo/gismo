/** @file BaseExpression.h

    @brief

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M.Verhelst
*/

#include <gismo.h>

using namespace gismo;

namespace Expr
{

    template <class T>
    class ScalarExpression : public BaseExpression<ScalarExpression<T>, T, 0, true>
    {
        using Base = BaseExpression<ScalarExpression<T>, T, 0, true>;
    public:
        typedef T Scalar; // Define Scalar type for this expression
        static constexpr int order = Base::order;
    private:
        Scalar m_value; // Store scalar value directly
    public:
        // Explicit constructor for a scalar value
        explicit ScalarExpression(Scalar value)
        :
        BaseExpression<ScalarExpression<Scalar>,Scalar, 0>({}),    // Pass empty array for sizes_ for order 0 (scalar)
        m_value(value)
        {}

        gsMatrix<Scalar> eval(const index_t k) const
        {
            // For a scalar, k might be ignored or used as a parameter if it's a field.
            // Returns a 1x1 gsMatrix containing the scalar value.
            gsMatrix<Scalar> result(1, 1);
            result(0, 0) = m_value;
            return result;
        }

    };
}