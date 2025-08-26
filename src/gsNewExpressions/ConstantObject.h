/** @file BaseObject.h

@brief

This file is part of the G+Smo library.

This Source Code Form is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this
file, You can obtain one at http://mozilla.org/MPL/2.0/.

Author(s): H.M.Verhelst
*/


#pragma once

namespace gismo
{
namespace Expr
{

template <class T, size_t _order>
struct ExpressionTraits<ConstantObject<T, _order>>
{
    using Base = BaseObject<T, _order, true, Space::None>;
    typedef typename ExpressionTraits<Base>::Scalar Scalar;
    static constexpr size_t order = ExpressionTraits<Base>::order;  // Use template parameter, not Base
    static constexpr size_t space = ExpressionTraits<Base>::space;
    static constexpr size_t deriv = ExpressionTraits<Base>::deriv;
    static constexpr bool isConstant = ExpressionTraits<Base>::isConstant;
};

template <class T, size_t _order>
class ConstantObject : public BaseObject<T, _order, true, Space::None>
{
    using Base = BaseObject<T, _order, true, Space::None>;
public:

    typedef typename ExpressionTraits<ConstantObject<T, _order>>::Scalar Scalar;
    static constexpr size_t order = ExpressionTraits<ConstantObject<T, _order>>::order;
    static constexpr size_t space = ExpressionTraits<ConstantObject<T, _order>>::space;
    static constexpr size_t deriv = ExpressionTraits<ConstantObject<T, _order>>::deriv;
    static constexpr bool isConstant = ExpressionTraits<ConstantObject<T, _order>>::isConstant;

private:
    gsMatrix<Scalar> m_value;

public:
    ConstantObject(const std::array<size_t, order> & input_sizes)
    :
    BaseObject<Scalar, _order, true, Space::None>(0, input_sizes)
    {
    }

    gsMatrix<Scalar> eval(const index_t) const
    {
        return m_value;
    }

    void setValue(const gsMatrix<T> & value) { m_value = value;}

    void parse(gismo::ExpressionHelper<Scalar> &) const {}

    void print(std::ostream & os) const
    {
        _print_impl<order>(os);
    }

protected:

    template <size_t _ORDER>
    typename std::enable_if<_ORDER==0,void>::type
    _print_impl(std::ostream & os) const { os<<"α"; }

    template <size_t _ORDER>
    typename std::enable_if<_ORDER==1,void>::type
    _print_impl(std::ostream & os) const { os<<"a"; }

    template <size_t _ORDER>
    typename std::enable_if<_ORDER!=0 && _ORDER!=1,void>::type
    _print_impl(std::ostream & os) const { os<<"A"; }
};
}//namespace Expr
}//namespace gismo