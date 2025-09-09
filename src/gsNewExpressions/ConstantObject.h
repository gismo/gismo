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
    typedef T Scalar;
    static constexpr size_t order = _order;
    static constexpr size_t space = Space::None;
    static constexpr size_t deriv = 0;
    static constexpr bool isConstant = true;
};

template <class T, size_t _order>
class ConstantObject : public BaseObject<ConstantObject<T, _order>>
{
    using Base = BaseObject<ConstantObject<T, _order>>;
    typedef typename Base::Scalar Scalar;
    using Base::order;
    using Base::sizes_;

private:
    gsMatrix<Scalar> m_value;

public:
    ConstantObject(const std::array<size_t, order> & input_sizes, std::string label = "a")
    :
    Base(0, input_sizes, label)
    {
        initSize();
        m_value.setZero();
    }

    explicit ConstantObject(const Scalar & value)
    :
    Base(0, std::array<size_t, order>{}, std::to_string(value))
    {
        m_value.resize(1,1);
        m_value.setConstant(value);
    }

    gsMatrix<Scalar> eval(const index_t) const
    {
        return m_value;
    }

    void setValue(const gsMatrix<T> & value) { m_value = value;}
    void setConstant(const T value) { m_value.setConstant(value);}
    void parse(gismo::ExpressionHelper<Scalar> &) const {}

    void print(std::ostream & os) const
    {
        os<<Base::label_;
    }

    size_t domainDim() const
    {
        GISMO_NO_IMPLEMENTATION;
    }

protected:
    void initSize()
    {
        if (order==0)
            m_value.resize(1,1);
        else if (order==1)
            m_value.resize(sizes_[0],1);
        else if (order==2)
            m_value.resize(sizes_[0],sizes_[1]);
        else
            GISMO_ERROR("ConstantObject only implemented for order 0, 1, or 2");
    }
};



}//namespace Expr
}//namespace gismo