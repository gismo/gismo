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

template <class T, size_t _Order>
struct ExpressionTraits<ConstantObject<T, _Order>>
{
    typedef T Scalar;
    static constexpr size_t Order = _Order;
    static constexpr size_t Space = SpaceType::None;
    static constexpr size_t Deriv = 0;
    static constexpr bool IsConstant = true;
};

template <class T, size_t _Order>
class ConstantObject : public BaseObject<ConstantObject<T, _Order>>
{
    using Base = BaseObject<ConstantObject<T, _Order>>;

public:
    // Expose the static traits publicly
    using Base::Order;
    using Base::Space;
    using Base::Deriv;
    using Base::IsConstant;
    using Base::sizes_;
    typedef typename Base::Scalar Scalar;

private:
    gsMatrix<Scalar> m_value;

public:
    ConstantObject(const std::array<size_t, Order> & input_sizes, std::string label = "a")
    :
    Base(0, input_sizes, label)
    {
        initSize();
        m_value.setZero();
    }

    explicit ConstantObject(const Scalar & value)
    :
    Base(0, std::array<size_t, Order>{}, std::to_string(value))
    {
        m_value.resize(1,1);
        m_value.setConstant(value);
    }

    template <int _Rows>
    explicit ConstantObject(const gsVector<T, _Rows> & value)
    :
    Base(0, std::array<size_t, 1>{value.rows()}, "c")
    {
        m_value = value;
    }

    template <int _Rows, int _Cols>
    explicit ConstantObject(const gsMatrix<T, _Rows, _Cols> & value)
    :
    Base(0, std::array<size_t, 2>{value.rows(), value.cols()}, "C")
    {
        m_value = value;
    }

    void setValue(const gsMatrix<T> & value) { m_value = value;}
    void setConstant(const T value) { m_value.setConstant(value);}
    void parse(gismo::ExpressionHelper<Scalar> &) const {}

    ExpressionValue<Scalar> eval(const index_t) const
    {
        return ExpressionValue<Scalar>(m_value);
    }

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
        if (Order==0)
            m_value.resize(1,1);
        else if (Order==1)
            m_value.resize(sizes_[0],1);
        else if (Order==2)
            m_value.resize(sizes_[0],sizes_[1]);
        else
            GISMO_ERROR("ConstantObject only implemented for order 0, 1, or 2");
    }
};



}//namespace Expr
}//namespace gismo